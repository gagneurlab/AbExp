import os
import sys
import warnings

import pytest

from abexp.utils.common import set_polars_memory_env

# 4000 MB, as Snakemake passes it in resources.mem_mb
LIMIT = 4000 * 10**6


@pytest.fixture
def no_polars(monkeypatch):
    """Clears both variables and hides an import of polars by other tests."""
    monkeypatch.delenv("POLARS_OVERRIDE_TOTAL_MEMORY", raising=False)
    monkeypatch.delenv("POLARS_OVERRIDE_TOTAL_MEMORY_MB", raising=False)
    monkeypatch.delitem(sys.modules, "polars", raising=False)


def test_set_polars_memory_env_sets_both_variables(no_polars):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        set_polars_memory_env(LIMIT)

    assert os.environ["POLARS_OVERRIDE_TOTAL_MEMORY"] == "4000000000"
    assert os.environ["POLARS_OVERRIDE_TOTAL_MEMORY_MB"] == "4000"


def test_set_polars_memory_env_keeps_set_variables(no_polars, monkeypatch):
    monkeypatch.setenv("POLARS_OVERRIDE_TOTAL_MEMORY", "8000000000")
    monkeypatch.setenv("POLARS_OVERRIDE_TOTAL_MEMORY_MB", "8000")

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        set_polars_memory_env(LIMIT)

    assert os.environ["POLARS_OVERRIDE_TOTAL_MEMORY"] == "8000000000"
    assert os.environ["POLARS_OVERRIDE_TOTAL_MEMORY_MB"] == "8000"


def test_set_polars_memory_env_warns_after_polars_import(no_polars, monkeypatch):
    monkeypatch.setitem(sys.modules, "polars", object())

    with pytest.warns(UserWarning, match="runs after the import of polars"):
        set_polars_memory_env(LIMIT)

    assert os.environ["POLARS_OVERRIDE_TOTAL_MEMORY"] == "4000000000"
    assert os.environ["POLARS_OVERRIDE_TOTAL_MEMORY_MB"] == "4000"
