Oxford detector elevation and sample tilt
=========================================

When constructing an indexer from an Oxford ``.h5oina`` file, leaving
``camElev=None`` (the default) reads the selected acquisition's optional
``EBSD/Header/Detector Orientation Euler`` dataset. Oxford stores Bunge Euler
angles in radians. The scalar camera elevation is obtained as
``degrees(Phi) - 90``, where ``Phi`` is the second angle::

    from pyebsdindex import ebsd_index

    indexer = ebsd_index.EBSDIndexer(filename="scan.h5oina")
    print(indexer.camElev)  # Degrees

An explicit ``camElev`` always takes precedence, including zero. To retain
the previous default even for files containing detector metadata, pass
``camElev=5.3``. Missing metadata, other file formats, and array-only inputs
retain the fallback of 5.3 degrees. Malformed or non-finite optional Oxford
detector metadata emits a warning and uses that fallback.

Similarly, ``sampleTilt=None`` reads the selected acquisition's optional
``EBSD/Header/Tilt Angle`` and converts radians to degrees. An explicit
``sampleTilt`` overrides metadata, including zero. If the field is absent,
or for array-only inputs and formats without tilt metadata, the previous
70 degree fallback is retained. Invalid or non-finite metadata emits a
warning and uses that fallback. For example::

    indexer = ebsd_index.EBSDIndexer(filename="scan.h5oina")
    print(indexer.sampleTilt)  # Degrees, automatically read from Tilt Angle
    indexer = ebsd_index.EBSDIndexer(filename="scan.h5oina", sampleTilt=70.0)

These defaults apply to ``index_pats`` and ``index_pats_distributed``
when they construct a new indexer. Reusing an indexer, including calling
``update_file``, retains its resolved geometry; construct a new indexer or
explicitly update ``indexer.camElev`` and ``indexer.sampleTilt`` when changing
acquisition geometry.

This is the elevation component of the existing scalar geometry model,
not support for the full three-angle detector orientation. Detector azimuth
and roll are not applied. Projection centers must still be supplied separately.
The scalar sample tilt assumes the existing tilt-axis convention; this change
does not apply ``Tilt Axis``, specimen orientation, or scanning rotation.

Projection center selection
---------------------------

The Oxford reader exposes the stored ``Pattern Center X``, ``Pattern Center Y``
and ``Detector Distance`` values through two explicit ``PC`` modes::

    # One arithmetic mean PC for every pattern (lower memory and simpler GPU work)
    fast = ebsd_index.EBSDIndexer(filename="scan.h5oina", PC="file_mean")

    # The stored PCx, PCy and detector distance for every pattern
    calibrated = ebsd_index.EBSDIndexer(
        filename="scan.h5oina", PC="file_per_pattern"
    )

Both modes require an Oxford ``.h5oina`` file with one finite PC triplet per
pattern and use the Oxford PC convention. Numeric values such as
``PC=[0.5, 0.5, 0.6]`` retain their existing meaning. The default remains the
historical fixed PC, so existing scripts do not silently change behavior.
For a subset starting at ``patstart``, ``file_per_pattern`` selects the matching
rows from the file. The distributed indexing entry point supports the same modes.
Invalid or incomplete file PC metadata raises ``ValueError`` instead of silently
falling back to a different geometry. The arithmetic mean is a convenience mode;
its speed advantage depends on detector size and the rest of the indexing setup.

Mapsweeper calibration arrays can be selected with ``PC="mapsweeper"`` (one
mean PC) or ``PC="mapsweeper_per_pattern"`` (the PC for each pattern). These
modes read ``/<acquisition>/Data Processing/Data/Pattern Center X`` and the
corresponding ``Pattern Center Y`` and ``Detector Distance`` datasets.

The metadata convention is specified in the
`Oxford H5OINA specification <https://github.com/oinanoanalysis/h5oina/blob/master/H5OINAFile.md>`_.
