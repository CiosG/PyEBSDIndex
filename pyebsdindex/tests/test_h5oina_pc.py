"""File-derived projection centers, including nonzero batch offsets."""
import h5py
import numpy as np
import pytest

from pyebsdindex import ebsd_index, ebsd_pattern, _ray_installed


PC_NAMES = ('Pattern Center X', 'Pattern Center Y', 'Detector Distance')


@pytest.fixture
def pc_file(tmp_path):
    path = tmp_path / 'pc.h5oina'
    pc = np.array([[.40, .35, .50], [.42, .36, .52], [.44, .37, .54],
                   [.46, .38, .56], [.48, .39, .58], [.50, .40, .60]])
    with h5py.File(path, 'w') as f:
        f['Format Version'] = [b'6.0']
        d = f.require_group('1/EBSD/Data')
        d['Processed Patterns'] = np.broadcast_to(np.arange(6, dtype='u1')[:, None, None], (6, 24, 32))
        for i, key in enumerate(PC_NAMES):
            d[key] = pc[:, i]
        for key in ('Beam Position X', 'Beam Position Y'):
            d[key] = np.arange(6)
        h = f.require_group('1/EBSD/Header')
        for key, value in [('X Cells', 3), ('Y Cells', 2), ('X Step', 1), ('Y Step', 1)]:
            h[key] = [value]
    return path, pc


def capture_pc(monkeypatch):
    captured = []
    def detect(self, pats, PC, **kwargs):
        captured.append((pats.copy(), PC.copy()))
        return None, None
    monkeypatch.setattr(ebsd_index.EBSDIndexer, '_detectbands', detect)
    monkeypatch.setattr(ebsd_index.EBSDIndexer, '_indexbandsphase',
                        lambda self, *args, **kwargs: (None, None))
    return captured


@pytest.mark.parametrize('mode', ['file_mean', 'file_per_pattern'])
def test_public_wrapper_nonzero_start(pc_file, monkeypatch, mode):
    path, pc = pc_file
    captured = capture_pc(monkeypatch)
    ebsd_index.index_pats(filename=path, PC=mode, patstart=2, npats=2, useCPU=True)
    pats, actual = captured[-1]
    np.testing.assert_array_equal(pats[:, 0, 0], [2, 3])
    expected = np.tile(pc.mean(axis=0), (2, 1)) if mode == 'file_mean' else pc[2:4]
    np.testing.assert_allclose(actual, expected)


def test_repeated_batches_and_end_of_file(pc_file, monkeypatch):
    path, pc = pc_file
    captured = capture_pc(monkeypatch)
    idx = ebsd_index.EBSDIndexer(filename=path, PC='file_per_pattern', useCPU=True)
    for start, count in [(0, 2), (2, 2), (4, 20), (3, -1)]:
        result = idx.index_pats(patstart=start, npats=count)
        expected = pc[start:] if count == -1 else pc[start:start+count]
        np.testing.assert_allclose(captured[-1][1], expected)
        assert result[-1] == len(expected)


def test_numeric_override_and_legacy_default(pc_file, monkeypatch):
    path, pc = pc_file
    captured = capture_pc(monkeypatch)
    legacy = ebsd_index.EBSDIndexer(filename=path, useCPU=True)
    np.testing.assert_allclose(legacy.PC, [.471659, .675044, .630139])
    idx = ebsd_index.EBSDIndexer(filename=path, PC='file_per_pattern', useCPU=True)
    idx.index_pats(patstart=2, npats=2, PC=[.5, .5, .6])
    np.testing.assert_allclose(captured[-1][1], [[.5, .5, .6]] * 2)
    idx.index_pats(patstart=2, npats=2, PC=pc[:2])
    np.testing.assert_allclose(captured[-1][1], pc[:2])
    legacy.index_pats(patstart=2, npats=2, PC='file_per_pattern')
    np.testing.assert_allclose(captured[-1][1], pc[2:4])


def test_selected_acquisition_and_update(pc_file):
    path, pc = pc_file
    with h5py.File(path, 'r+') as f:
        f.copy('1', '2')
        f['2/EBSD/Data/Pattern Center X'][...] = pc[:, 0] + .1
    reader = ebsd_pattern.get_pattern_file_obj([str(path), '2/EBSD/Data/Processed Patterns'])
    expected = pc.copy()
    expected[:, 0] += .1
    np.testing.assert_allclose(reader.read_pc(), expected)
    idx = ebsd_index.EBSDIndexer(filename=path, PC='file_per_pattern', useCPU=True)
    idx.update_file(reader)
    np.testing.assert_allclose(idx._fillPCarray(None, 2, patstart=2), expected[2:4])


def test_mapsweeper_pc_source(pc_file, monkeypatch):
    path, pc = pc_file
    mapsweeper = pc + [.07, -.03, .11]
    with h5py.File(path, 'r+') as f:
        d = f.require_group('1/Data Processing/Data')
        for i, key in enumerate(PC_NAMES):
            d[key] = mapsweeper[:, i]
    reader = ebsd_pattern.get_pattern_file_obj(path)
    np.testing.assert_allclose(reader.read_pc(source='mapsweeper'), mapsweeper)
    captured = capture_pc(monkeypatch)
    idx = ebsd_index.EBSDIndexer(filename=path, PC='mapsweeper_per_pattern', useCPU=True)
    idx.index_pats(patstart=2, npats=2)
    np.testing.assert_allclose(captured[-1][1], mapsweeper[2:4])


@pytest.mark.parametrize('kind', ['missing', 'nan', 'infinite', 'zero_dd', 'negative_dd', 'shape'])
def test_bad_metadata_is_not_silently_replaced(pc_file, kind):
    path, pc = pc_file
    with h5py.File(path, 'r+') as f:
        d = f['1/EBSD/Data']
        if kind in ('missing', 'shape'):
            del d[PC_NAMES[0]]
            if kind == 'shape':
                d[PC_NAMES[0]] = pc[:2, 0]
        else:
            values = {'nan': np.nan, 'infinite': np.inf, 'zero_dd': 0, 'negative_dd': -1}
            d[PC_NAMES[2]][0] = values[kind]
    with pytest.raises(ValueError, match='Oxford PC metadata'):
        ebsd_index.EBSDIndexer(filename=path, PC='file_per_pattern', useCPU=True)
    # Invalid file metadata must not affect explicit numeric PCs.
    idx = ebsd_index.EBSDIndexer(filename=path, PC=pc[0], useCPU=True)
    np.testing.assert_array_equal(idx.PC, pc[0])


def test_column_vectors(pc_file):
    path, pc = pc_file
    with h5py.File(path, 'r+') as f:
        d = f['1/EBSD/Data']
        for i, key in enumerate(PC_NAMES):
            del d[key]
            d[key] = pc[:, i:i+1]
    np.testing.assert_array_equal(ebsd_pattern.get_pattern_file_obj(path).read_pc(), pc)


def test_modes_require_file_and_oxford_convention(pc_file):
    path, pc = pc_file
    with pytest.raises(ValueError, match='require an Oxford'):
        ebsd_index.EBSDIndexer(PC='file_mean', useCPU=True)
    with pytest.raises(ValueError, match='OXFORD vendor'):
        ebsd_index.EBSDIndexer(filename=path, PC='file_mean', vendor='EDAX', useCPU=True)
    with pytest.raises(ValueError, match='PC mode'):
        ebsd_index.EBSDIndexer(filename=path, PC='typo', useCPU=True)
    idx = ebsd_index.EBSDIndexer(filename=path, PC='file_per_pattern', useCPU=True)
    with pytest.raises(ValueError, match='range'):
        idx._fillPCarray(None, 3, patstart=5)


def test_supplied_batch_global_offset(pc_file, monkeypatch):
    path, pc = pc_file
    captured = capture_pc(monkeypatch)
    idx = ebsd_index.EBSDIndexer(filename=path, PC='file_per_pattern', useCPU=True)
    idx.index_pats(patsin=np.zeros((2, 24, 32)), patstart=3)
    np.testing.assert_array_equal(captured[-1][1], pc[3:5])


@pytest.mark.parametrize('mode', ['file_mean', 'file_per_pattern'])
def test_saved_pc_mode(pc_file, tmp_path, mode):
    path, pc = pc_file
    idx = ebsd_index.EBSDIndexer(filename=path, PC=mode, useCPU=True)
    saved = tmp_path / 'indexer.pyindx'
    idx.saveindexer(saved)
    restored = ebsd_index.restoreindexer(saved)
    assert restored.PC_file_mode == mode
    np.testing.assert_allclose(restored._fillPCarray(None, 2, patstart=2),
                               idx._fillPCarray(None, 2, patstart=2))


def test_actual_indexing_matches_explicit_pc(pc_file, pattern_al_sim_20kv):
    path, pc = pc_file
    with h5py.File(path, 'r+') as f:
        d = f['1/EBSD/Data']
        del d['Processed Patterns']
        d['Processed Patterns'] = np.repeat(pattern_al_sim_20kv[None], 6, axis=0)
    auto = ebsd_index.EBSDIndexer(filename=path, PC='file_per_pattern', useCPU=True)
    manual = ebsd_index.EBSDIndexer(filename=path, useCPU=True)
    result = auto.index_pats(patstart=2, npats=2)[0]
    expected = manual.index_pats(patstart=2, npats=2, PC=pc[2:4])[0]
    np.testing.assert_allclose(result['quat'], expected['quat'])
    np.testing.assert_array_equal(result['nmatch'], expected['nmatch'])


@pytest.mark.skipif(not _ray_installed, reason='ray is not installed')
def test_distributed_nonzero_start(pc_file, pattern_al_sim_20kv):
    path, pc = pc_file
    with h5py.File(path, 'r+') as f:
        d = f['1/EBSD/Data']
        del d['Processed Patterns']
        d['Processed Patterns'] = np.repeat(pattern_al_sim_20kv[None], 6, axis=0)
    idx = ebsd_index.EBSDIndexer(filename=path, PC='file_per_pattern', useCPU=True)
    expected = idx.index_pats(patstart=2, npats=3)[0]
    result = ebsd_index.index_pats_distributed(filename=path, ebsd_indexer_obj=idx,
                                              patstart=2, npats=3, chunksize=2, ncpu=1)[0]
    np.testing.assert_allclose(result['quat'], expected['quat'])
