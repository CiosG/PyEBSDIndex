"""Index an H5OINA file in this folder and write an .ang map.

Edit the settings below before running the script.
"""
from pathlib import Path
import shutil

import h5py
import numpy as np

from pyebsdindex import ebsd_index, rotlib


# ------------------------- ustawienia użytkownika -------------------------
FOLDER = Path(__file__).resolve().parent.parent  # folder z plikiem .h5oina
INPUT_FILE = "Ni Max new Specimen 1 10 kV Map Data 3.h5oina"
OUTPUT_FILE = "results_pyebsdindex/Ni Max new Specimen 1 10 kV Map Data 3_indexed.h5oina"
OUTPUT_FORMAT = "h5oina"                 # "h5oina" albo "ang"

VENDOR = "OXFORD"
PHASELIST = ["FCC"]                       # np. ["FCC", "BCC"]
PC_MODE = "mapsweeper_per_pattern"       # mapsweeper albo mapsweeper_per_pattern
USE_CPU = True                             # False = GPU/OpenCL, jeśli dostępne
SAMPLE_TILT = None                         # None = odczyt z H5OINA; np. 70.0
CAM_ELEV = None                            # None = odczyt z H5OINA; np. 5.3
PATSTART = 0
NPATS = -1                                 # -1 = wszystkie patterny
# ---------------------------------------------------------------------------

folder = Path(FOLDER).expanduser().resolve()
inputs = [folder / INPUT_FILE] if INPUT_FILE else sorted(folder.glob("*.h5oina"))
if not inputs:
    raise FileNotFoundError(f"Nie znaleziono pliku .h5oina w {folder}")

source = inputs[0]
extension = ".h5oina" if OUTPUT_FORMAT == "h5oina" else ".ang"
output = (folder / OUTPUT_FILE if OUTPUT_FILE else
          source.with_name(source.stem + "_indexed" + extension))
if output.resolve() == source.resolve():
    raise ValueError("OUTPUT_FILE musi być inny niż plik wejściowy")

data, _, indexer = ebsd_index.index_pats(
    filename=str(source),
    vendor=VENDOR,
    phaselist=PHASELIST,
    PC=PC_MODE,
    sampleTilt=SAMPLE_TILT,
    camElev=CAM_ELEV,
    patstart=PATSTART,
    npats=NPATS,
    useCPU=USE_CPU,
    return_indexer_obj=True,
)

if OUTPUT_FORMAT == "ang":
    from pyebsdindex import ebsdfile
    ebsdfile.writeang(str(output), indexer, data)
elif OUTPUT_FORMAT == "h5oina":
    # Zachowaj cały oryginalny plik i nadpisz standardowe pola Euler.
    shutil.copy2(source, output)
    best = data[-1]
    n = len(best)
    eulers = np.asarray(rotlib.qu2eu(best["quat"]), dtype=np.float32)
    # H5OINA phase IDs are one-based; PyEBSDIndex uses zero-based IDs.
    phase = np.where(best["phase"] >= 0, best["phase"] + 1, 0).astype(np.uint8)
    # PyEBSDIndex fit is reported in degrees; H5OINA MAD is stored in radians.
    mad = np.deg2rad(np.asarray(best["fit"], dtype=np.float32))
    acquisition = indexer.fID.h5patdatpth.strip('/').split('/')[0]
    with h5py.File(output, "r+") as h5:
        for data_group in (f"/{acquisition}/EBSD/Data",
                           f"/{acquisition}/Data Processing/Data"):
            dataset = data_group + "/Euler"
            if dataset in h5:
                stored = h5[dataset]
                if stored.ndim != 2 or stored.shape[1] != 3:
                    raise ValueError(f"Nieprawidłowy dataset Euler: {dataset}")
                if PATSTART + n > stored.shape[0]:
                    raise ValueError(f"Za mało wierszy w dataset Euler: {dataset}")
                stored[PATSTART:PATSTART + n] = eulers
                stored.attrs["Unit"] = "rad"
            else:
                h5.create_dataset(dataset, data=eulers)
                h5[dataset].attrs["Unit"] = "rad"
            for name, values, unit in (("Phase", phase, None),
                                       ("Mean Angular Deviation", mad, "rad")):
                dataset = data_group + "/" + name
                if dataset in h5:
                    stored = h5[dataset]
                    if stored.ndim != 1 or PATSTART + n > stored.shape[0]:
                        raise ValueError(f"Nieprawidłowy dataset {name}: {dataset}")
                    stored[PATSTART:PATSTART + n] = values
                else:
                    h5.create_dataset(dataset, data=values)
                if unit is not None:
                    h5[dataset].attrs["Unit"] = unit
else:
    raise ValueError('OUTPUT_FORMAT musi być "h5oina" albo "ang"')
print(f"Zapisano: {output}")
