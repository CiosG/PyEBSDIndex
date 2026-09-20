"""Synthetic Oxford metadata tests; no experimental acquisition is required."""
import h5py
import numpy as np
import pytest

from pyebsdindex import ebsd_pattern, ebsd_index


@pytest.fixture
def h5oina(tmp_path):
    def create(elevation=8.0, orientation=None, name='scan.h5oina'):
        path = tmp_path / name
        with h5py.File(path, 'w') as f:
            f['Format Version'] = [b'6.0']
            f.create_dataset('1/EBSD/Data/Processed Patterns', shape=(1, 24, 32), dtype='u1')
            h = f.require_group('1/EBSD/Header')
            for key in ('X Cells', 'Y Cells', 'X Step', 'Y Step'):
                h[key] = [1]
            if orientation is not None:
                h['Detector Orientation Euler'] = orientation
            elif elevation is not None:
                h['Detector Orientation Euler'] = np.radians([[0, 90 + elevation, 0]])
        return path
    return create


@pytest.mark.parametrize('elevation', [-5.0, 0.0, 8.0, 12.0])
def test_read_elevation_and_reference_transform(h5oina, elevation):
    path = h5oina(elevation)
    reader = ebsd_pattern.get_pattern_file_obj(path)
    assert reader.camElev == pytest.approx(elevation)
    indexer = ebsd_index.EBSDIndexer(filename=path, useCPU=True)
    assert indexer.camElev == pytest.approx(elevation)
    # The existing Oxford reference transform rotates only about x.
    half_angle = np.radians(-(90 - 70 + elevation)) / 2
    np.testing.assert_allclose(indexer._detector2refframe(),
                               [np.cos(half_angle), np.sin(half_angle), 0, 0])


@pytest.mark.parametrize('override', [0.0, 5.3, -3.0])
def test_explicit_override(h5oina, override):
    indexer = ebsd_index.EBSDIndexer(filename=h5oina(), camElev=override, useCPU=True)
    assert indexer.camElev == override


def test_missing_metadata_and_array_input_keep_legacy_default(h5oina):
    path = h5oina(elevation=None)
    assert ebsd_pattern.get_pattern_file_obj(path).camElev is None
    assert ebsd_index.EBSDIndexer(filename=path, useCPU=True).camElev == 5.3
    assert ebsd_index.EBSDIndexer(patDim=(24, 32), useCPU=True).camElev == 5.3


@pytest.mark.parametrize('orientation', [[[0, np.nan, 0]], [[0, np.inf, 0]], [0, 1], [[0, 1, 0], [0, 1, 0]]])
def test_invalid_optional_metadata_warns_and_falls_back(h5oina, orientation):
    with pytest.warns(UserWarning, match='Invalid Oxford Detector Orientation Euler'):
        indexer = ebsd_index.EBSDIndexer(filename=h5oina(orientation=orientation), useCPU=True)
    assert indexer.camElev == 5.3


def test_unread_file_object(h5oina):
    reader = ebsd_pattern.OXFORDOINA(h5oina())
    assert reader.patternH is None
    indexer = ebsd_index.EBSDIndexer(filename=reader, useCPU=True)
    assert indexer.camElev == pytest.approx(8)


def test_selected_acquisition_and_reset(h5oina):
    path = h5oina()
    with h5py.File(path, 'r+') as f:
        f.copy('1', '2')
        del f['2/EBSD/Header/Detector Orientation Euler']
    reader = ebsd_pattern.get_pattern_file_obj(path)
    assert reader.camElev == pytest.approx(8)
    reader.set_data_path('2/EBSD/Data/Processed Patterns')
    reader.read_header()
    assert reader.camElev is None
    assert ebsd_index.EBSDIndexer(filename=reader, useCPU=True).camElev == 5.3


@pytest.mark.parametrize('override,expected', [(None, 8), (0, 0), (5.3, 5.3)])
def test_index_pats_wrapper(h5oina, monkeypatch, override, expected):
    # Exercise public argument plumbing without indexing a blank pattern.
    monkeypatch.setattr(ebsd_index.EBSDIndexer, 'index_pats',
                        lambda self, **kwargs: (None, None, 0, 1))
    _, _, indexer = ebsd_index.index_pats(filename=h5oina(), camElev=override,
                                        useCPU=True, return_indexer_obj=True)
    assert indexer.camElev == pytest.approx(expected)


def test_existing_indexer_geometry_is_preserved(h5oina, monkeypatch):
    indexer = ebsd_index.EBSDIndexer(filename=h5oina(), useCPU=True)
    monkeypatch.setattr(ebsd_index.EBSDIndexer, 'index_pats',
                        lambda self, **kwargs: (None, None, 0, 1))
    _, _, reused = ebsd_index.index_pats(filename=h5oina(12, name='other.h5oina'),
        ebsd_indexer_obj=indexer, return_indexer_obj=True)
    assert reused.camElev == pytest.approx(8)


def test_flat_orientation(h5oina):
    reader = ebsd_pattern.get_pattern_file_obj(h5oina(orientation=np.radians([0, 98, 0])))
    assert reader.camElev == pytest.approx(8)
