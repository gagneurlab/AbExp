import hashlib

import pytest

from abexp.pangolin import MODEL_SHA256, download_model
from abexp.pangolin import weights

FILE = 'final.1.0.3.v2'


def sha256(data):
    return hashlib.sha256(data).hexdigest()


@pytest.fixture
def served(tmp_path, monkeypatch):
    """A folder that MODEL_URL points to, as a file URL, so that the tests need no network."""
    served = tmp_path / 'served'
    served.mkdir()
    monkeypatch.setattr(weights, 'MODEL_URL', served.as_uri() + '/{}')
    return served


def test_download_model(served, tmp_path, monkeypatch):
    (served / FILE).write_bytes(b'weights')
    monkeypatch.setitem(MODEL_SHA256, FILE, sha256(b'weights'))
    download_model(FILE, tmp_path / 'models')
    assert [path.name for path in (tmp_path / 'models').iterdir()] == [FILE]
    assert (tmp_path / 'models' / FILE).read_bytes() == b'weights'


def test_download_model_wrong_sha256(served, tmp_path):
    # the served file is not the published model; the downloaded file is deleted
    (served / FILE).write_bytes(b'other weights')
    with pytest.raises(ValueError) as e:
        download_model(FILE, tmp_path / 'models')
    part = tmp_path / 'models' / f'{FILE}.part'
    assert str(e.value) == f'The SHA-256 sum of {part} is {sha256(b"other weights")}, not {MODEL_SHA256[FILE]}'
    assert list((tmp_path / 'models').iterdir()) == []


def test_download_model_keeps_existing_file(served, tmp_path, monkeypatch):
    # nothing is served, so a download would fail
    (tmp_path / FILE).write_bytes(b'weights')
    monkeypatch.setitem(MODEL_SHA256, FILE, sha256(b'weights'))
    download_model(FILE, tmp_path)
    assert (tmp_path / FILE).read_bytes() == b'weights'


def test_download_model_existing_file_wrong_sha256(served, tmp_path):
    # the existing file stays
    (tmp_path / FILE).write_bytes(b'other weights')
    with pytest.raises(ValueError) as e:
        download_model(FILE, tmp_path)
    path = tmp_path / FILE
    assert str(e.value) == f'The SHA-256 sum of {path} is {sha256(b"other weights")}, not {MODEL_SHA256[FILE]}'
    assert (tmp_path / FILE).read_bytes() == b'other weights'
