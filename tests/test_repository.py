from pathlib import Path
import importlib.util

import pytest

CORE = Path(__file__).resolve().parents[1] / "humid_climatology_engine_v11.5.py"


def load_core():
    if not CORE.exists():
        pytest.skip(f"v11.5 core not present at {CORE}")
    spec = importlib.util.spec_from_file_location("humid_engine_v11_5", CORE)
    if spec is None or spec.loader is None:
        pytest.fail("Unable to load v11.5 core module")
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except ImportError as exc:
        pytest.skip(f"runtime dependency unavailable: {exc}")
    return mod


def test_public_release_identity():
    mod = load_core()
    assert str(mod.ENGINE_VERSION).startswith("11.5")


def test_temporal_state_has_33_bins():
    mod = load_core()
    assert mod.LEVEL_BINS["L1"] == (0, 1)
    assert mod.LEVEL_BINS["L2"] == (1, 9)
    assert mod.LEVEL_BINS["L3"] == (9, 33)


def test_primary_input_families_and_products():
    mod = load_core()
    assert tuple(mod.VARIABLES) == ("rh", "e", "r", "q")
    assert ("rh", "q") in tuple(mod.PAIRS)


def test_v11_5_cli_commands_are_present():
    source = CORE.read_text(encoding="utf-8")
    for command in ("selftest", "validate-input", "pilot", "run", "audit", "merge-audit", "report", "benchmark"):
        assert f'"{command}"' in source
