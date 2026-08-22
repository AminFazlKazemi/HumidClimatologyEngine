from pathlib import Path
import importlib.util
import numpy as np
import pytest

CORE = Path(__file__).resolve().parents[1] / "src" / "moisture_climatology_v6.py"

try:
    spec = importlib.util.spec_from_file_location("humid_engine", CORE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
except ImportError as exc:
    pytest.skip(f"runtime dependency unavailable: {exc}", allow_module_level=True)


def test_leap_day_mapping():
    assert mod.get_clim_doy(59, 1984) == 60
    assert mod.get_clim_doy(60, 1984) == 60
    assert mod.get_clim_doy(61, 1984) == 61
    assert mod.get_clim_doy(59, 1985) == 60
    assert mod.get_clim_doy(60, 1985) == 61


def test_moisture_physics():
    T = np.array([25.0], dtype=np.float32)
    Td = np.array([18.0], dtype=np.float32)
    P = np.array([1005.0], dtype=np.float32)
    d = mod.derive_moisture(T, Td, P)
    assert 0 <= d["rh"][0] <= 100
    assert d["e"][0] > 0
    assert d["r"][0] > 0
    assert 0 < d["q"][0] < d["r"][0]


def test_log_progress_single_definition():
    source = CORE.read_text(encoding="utf-8")
    assert source.count("def log_progress(") == 1
    assert 'log_progress("STAGE", stage=' not in source
