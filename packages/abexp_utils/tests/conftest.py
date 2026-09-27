import pytest
from pathlib import Path


@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(Path(request.fspath.dirname).parent.parent.parent)

@pytest.fixture
def output_dir():
    output_dir = Path('output/test/')
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir