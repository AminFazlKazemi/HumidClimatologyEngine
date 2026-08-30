#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Moisture Copula Parameter Production Runner v8.2
================================================

Production-oriented runner for large ERA5-Land grids.

v8.2 changes from v8.1
----------------------
1) FIXED climatological calendar mapping:
   - slot 59 is reserved
   - slot 60 pools Feb-28 and Feb-29
   - slot 61 = Mar-01
   - slot 366 = Dec-31
2) DOY 59 is never included in a processing task.
3) MAX_WORKERS is no longer overwritten inside the script.
4) Logging is protected against duplicate handlers.
5) Benchmark uses the real ProcessPoolExecutor path.
6) Benchmark can test multiple workers concurrently and reports:
      fits, elapsed, fit/sec, task throughput, efficiency.
7) Tile memory is bounded more tightly:
   each DOY is accumulated and fitted before moving to the next DOY.
8) Production checkpoint/resume, atomic outputs and SHA256 verification retained.
9) Scientific fitting contract is unchanged:
   fit_window_models() remains the supplied engine selector.
"""

from __future__ import annotations

import argparse
import ast
import gc
import hashlib
import json
import logging
import math
import os

# Do not overwrite user supplied worker settings.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import shutil
import sys
import tempfile
import time
from collections import OrderedDict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from netCDF4 import Dataset

# -----------------------------------------------------------------------------
# Engine import
# -----------------------------------------------------------------------------

ENGINE_DIR = Path(r"K:\kazemi\papers\temperature_interpolation\HumidClimatologyEngine")
ERA5_T2M_DIR = Path(r"F:\Kazemi\era5\land\T2m")
ERA5_D2M_DIR = Path(r"F:\Kazemi\era5\land\Dew_Point_Temperature")
ERA5_SP_DIR = Path(r"F:\Kazemi\era5\land\Surface_Pressure")
OUTPUT_DIR = Path(r"C:\c")

sys.path.insert(0, str(ENGINE_DIR))

from moisture_climatology_v8_0_FINAL_SINGLE_PASS import (  # noqa: E402
    START_YEAR,
    END_YEAR,
    DOY_COUNT,
    CONFIG_HASH as ENGINE_CONFIG_HASH,
    fit_window_models,
)

# =============================================================================
# Configuration
# =============================================================================

SCHEMA_VERSION = "8.2-COPULA-PRODUCTION"

MAX_WORKERS = int(os.environ.get(
    "MAX_WORKERS",
    max(1, (os.cpu_count() or 4) - 2)
))
MAX_IN_FLIGHT = int(os.environ.get(
    "MAX_IN_FLIGHT",
    max(1, MAX_WORKERS * 2)
))

DOY_CHUNK = max(1, int(os.environ.get("DOY_CHUNK", 3)))
TILE_LAT = max(1, int(os.environ.get("TILE_LAT", 16)))
TILE_LON = max(1, int(os.environ.get("TILE_LON", 16)))

DATASET_CACHE_SIZE = max(3, int(os.environ.get("DATASET_CACHE_SIZE", 18)))

FIT_MIN_OBS = max(1, int(os.environ.get("FIT_MIN_OBS", 30)))
MAX_PARAM_COUNT = 5

LAT_RANGE = None
LON_RANGE = None
STRIDE_LAT = max(1, int(os.environ.get("STRIDE_LAT", 1)))
STRIDE_LON = max(1, int(os.environ.get("STRIDE_LON", 1)))

HALF_WIDTH_DAYS = 2
WINDOW_HOURS = 121

NETCDF_LEVEL = min(9, max(0, int(os.environ.get("NETCDF_LEVEL", 4))))
NETCDF_SHUFFLE = True
OVERWRITE_CONSOLE_OUTPUT = os.environ.get("OVERWRITE_CONSOLE_OUTPUT", "0") == "1"

BENCHMARK_DAYS = max(1, int(os.environ.get("BENCHMARK_DAYS", 3)))
BENCHMARK_LAT = max(1, int(os.environ.get("BENCHMARK_LAT", TILE_LAT)))
BENCHMARK_LON = max(1, int(os.environ.get("BENCHMARK_LON", TILE_LON)))
BENCHMARK_TASKS = max(1, int(os.environ.get("BENCHMARK_TASKS", max(4, MAX_WORKERS * 2))))
BENCHMARK_WORKERS = max(1, int(os.environ.get("BENCHMARK_WORKERS", MAX_WORKERS)))

DOY_START_ENV = os.environ.get("DOY_START")
DOY_END_ENV = os.environ.get("DOY_END")
DOY_START = int(DOY_START_ENV) if DOY_START_ENV else 1
DOY_END = int(DOY_END_ENV) if DOY_END_ENV else DOY_COUNT

PROCESS_YEARS_ENV = os.environ.get("PROCESS_YEARS")
if PROCESS_YEARS_ENV:
    parsed_years = ast.literal_eval(PROCESS_YEARS_ENV)
    PROCESS_YEARS = sorted(set(int(y) for y in parsed_years))
else:
    PROCESS_YEARS = list(range(START_YEAR, END_YEAR + 1))

# =============================================================================
# Logging
# =============================================================================

logger = logging.getLogger("MoistureCopulaV82")
logger.setLevel(logging.INFO)
logger.propagate = False

if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(
        logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    )
    logger.addHandler(handler)

# =============================================================================
# Global worker state
# =============================================================================

_W_FILE_INDEX: dict[str, dict[int, dict[int, Path]]] | None = None
_W_DATASET_CACHE: OrderedDict[tuple[str, int, int], xr.Dataset] | None = None
_W_SELECTED_LAT_IDX: np.ndarray | None = None
_W_SELECTED_LON_IDX: np.ndarray | None = None
_W_SELECTED_LAT: np.ndarray | None = None
_W_SELECTED_LON: np.ndarray | None = None
_W_REFERENCE_ENGINE_HASH: str | None = None


# =============================================================================
# Generic utilities
# =============================================================================

def atomic_json_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=str(path.parent),
        text=True,
    )
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2, sort_keys=True)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(block_size), b""):
            h.update(block)
    return h.hexdigest()


def safe_float(value: Any) -> float:
    try:
        x = float(value)
        return x if np.isfinite(x) else float("nan")
    except Exception:
        return float("nan")


def params_to_fixed_array(best: dict[str, Any] | None) -> np.ndarray:
    out = np.full(MAX_PARAM_COUNT, np.nan, dtype=np.float32)
    if not best:
        return out

    name = str(best.get("name", ""))
    p = best.get("params", {}) or {}

    if name == "Normal":
        keys = ("loc", "scale")
    elif name.startswith("SkewNormal"):
        keys = ("shape", "loc", "scale")
    elif name == "PearsonIII":
        keys = ("skew", "loc", "scale")
    elif name == "Beta":
        keys = ("a", "b", "loc", "scale")
    elif name == "BimodalNormal":
        keys = ("w1", "mu1", "sigma1", "mu2", "sigma2")
    else:
        keys = ("loc", "scale", "shape", "a", "b")

    for i, key in enumerate(keys[:MAX_PARAM_COUNT]):
        out[i] = safe_float(p.get(key))
    return out


def ensure_env_for_worker() -> None:
    # Prevent BLAS/OpenMP oversubscription: N Python processes should not each
    # create many native threads.
    for key in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ.setdefault(key, "1")


# =============================================================================
# Calendar / climatological DOY
# =============================================================================

def is_leap_year(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def slot_for_date(d: date) -> int:
    """
    366-slot contract:
      1..58   Jan-01..Feb-27
      59      reserved
      60      Feb-28/Feb-29 composite
      61..366 Mar-01..Dec-31
    """
    native = d.timetuple().tm_yday

    if d.month == 2 and d.day in (28, 29):
        return 60

    if native <= 58:
        return native

    # After Feb-28, both leap and non-leap dates map to the same 366-slot
    # calendar, with slot 59 reserved.
    return native + 1 if not is_leap_year(d.year) else native


def target_dates_for_slot(slot: int, years: list[int]) -> list[date]:
    if slot == 59:
        return []

    if slot == 60:
        out: list[date] = []
        for y in years:
            out.append(date(y, 2, 28))
            if is_leap_year(y):
                out.append(date(y, 2, 29))
        return out

    if not 1 <= slot <= 366:
        raise ValueError(f"Invalid climatological slot: {slot}")

    # In a non-leap reference year, slot 61 must be Mar-01:
    # Jan-01 + 60 days = Mar-02 is wrong; Jan-01 + 59 days = Mar-01.
    ref = date(2001, 1, 1) + timedelta(days=slot - 2)
    month, day = ref.month, ref.day

    out: list[date] = []
    for y in years:
        # Reference year 2001 has no Feb-29, so normally this is unnecessary,
        # but retain the guard for clarity.
        if month == 2 and day == 29:
            continue
        out.append(date(y, month, day))
    return out


def validate_calendar_contract() -> None:
    checks = {
        1: date(2001, 1, 1),
        58: date(2001, 2, 27),
        60: date(2001, 2, 28),
        61: date(2001, 3, 1),
        62: date(2001, 3, 2),
        365: date(2001, 12, 30),
        366: date(2001, 12, 31),
    }
    for slot, expected in checks.items():
        got = target_dates_for_slot(slot, [2001])
        if not got or got[0] != expected:
            raise RuntimeError(
                f"Calendar contract failed: slot {slot} -> {got}; "
                f"expected {expected}"
            )

    leap_targets = target_dates_for_slot(60, [2000])
    if leap_targets != [date(2000, 2, 28), date(2000, 2, 29)]:
        raise RuntimeError(f"Slot 60 leap-day contract failed: {leap_targets}")

    if target_dates_for_slot(59, [2000]):
        raise RuntimeError("Reserved slot 59 must have no target dates.")


def valid_slots(start: int, end: int) -> list[int]:
    validate_doy_range(start, end)
    return [s for s in range(start, end + 1) if s != 59]


def chunk_ranges_excluding_reserved(start: int, end: int, chunk: int) -> list[tuple[int, int]]:
    """
    Return half-open [d0,d1) chunks, never containing slot 59.
    """
    slots = valid_slots(start, end)
    ranges: list[tuple[int, int]] = []
    k = 0
    while k < len(slots):
        part = slots[k:k + chunk]
        # Because 59 is removed, each returned part may need to be represented
        # as one contiguous range. Split if needed.
        seg_start = part[0]
        prev = part[0]
        for s in part[1:]:
            if s != prev + 1:
                ranges.append((seg_start, prev + 1))
                seg_start = s
            prev = s
        ranges.append((seg_start, prev + 1))
        k += chunk
    return ranges


def validate_doy_range(start: int, end: int) -> None:
    if not (1 <= start <= DOY_COUNT):
        raise ValueError(f"DOY_START must be in 1..{DOY_COUNT}")
    if not (1 <= end <= DOY_COUNT):
        raise ValueError(f"DOY_END must be in 1..{DOY_COUNT}")
    if start > end:
        raise ValueError("DOY_START must be <= DOY_END")


# =============================================================================
# File index and grid
# =============================================================================

def extract_year_month(path: Path, year: int) -> int | None:
    import re

    m = re.search(r"(?<!\d)(\d{4})(\d{2})(?!\d)", path.name)
    if not m:
        return None
    y, mon = int(m.group(1)), int(m.group(2))
    if y != year or not 1 <= mon <= 12:
        return None
    return mon


def build_file_index(year: int, folder: Path) -> dict[int, Path]:
    index: dict[int, Path] = {}

    for p in sorted(folder.glob(f"*{year}*.nc")):
        mon = extract_year_month(p, year)
        if mon is None:
            continue
        if mon in index:
            raise RuntimeError(f"Duplicate month {year}-{mon:02d} in {folder}: {p}")
        index[mon] = p

    missing = sorted(set(range(1, 13)) - set(index))
    if missing:
        raise RuntimeError(f"Missing months for {year} in {folder}: {missing}")

    return index


def build_all_file_indices(years: list[int]) -> dict[str, dict[int, dict[int, Path]]]:
    return {
        "t2m": {y: build_file_index(y, ERA5_T2M_DIR) for y in years},
        "d2m": {y: build_file_index(y, ERA5_D2M_DIR) for y in years},
        "sp": {y: build_file_index(y, ERA5_SP_DIR) for y in years},
    }


def get_grid() -> tuple[np.ndarray, np.ndarray]:
    idx = build_file_index(START_YEAR, ERA5_T2M_DIR)
    with xr.open_dataset(
        idx[1],
        engine="netcdf4",
        decode_times=True,
        mask_and_scale=True,
        cache=False,
    ) as ds:
        lat = np.asarray(ds.latitude.values)
        lon = np.asarray(ds.longitude.values)
    return lat, lon


def select_axis_indices(
    values: np.ndarray,
    value_range: tuple[float, float] | None,
    stride: int,
) -> np.ndarray:
    if stride < 1:
        raise ValueError("stride must be >= 1")

    v = np.asarray(values)
    if value_range is None:
        base = np.arange(v.size, dtype=np.int32)
    else:
        lo, hi = value_range
        if lo > hi:
            lo, hi = hi, lo
        base = np.flatnonzero((v >= lo) & (v <= hi)).astype(np.int32)

    return base[::stride]


def configure_grid() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lat, lon = get_grid()
    lat_idx = select_axis_indices(lat, LAT_RANGE, STRIDE_LAT)
    lon_idx = select_axis_indices(lon, LON_RANGE, STRIDE_LON)
    return lat, lon, lat_idx, lon_idx


# =============================================================================
# Xarray dataset cache
# =============================================================================

def open_cached_dataset(kind: str, year: int, month: int) -> xr.Dataset:
    global _W_DATASET_CACHE, _W_FILE_INDEX

    assert _W_DATASET_CACHE is not None
    assert _W_FILE_INDEX is not None

    key = (kind, year, month)
    ds = _W_DATASET_CACHE.get(key)

    if ds is not None:
        _W_DATASET_CACHE.move_to_end(key)
        return ds

    path = _W_FILE_INDEX[kind][year][month]
    ds = xr.open_dataset(
        path,
        engine="netcdf4",
        decode_times=True,
        mask_and_scale=True,
        cache=False,
    )

    _W_DATASET_CACHE[key] = ds
    _W_DATASET_CACHE.move_to_end(key)

    while len(_W_DATASET_CACHE) > DATASET_CACHE_SIZE:
        _, old = _W_DATASET_CACHE.popitem(last=False)
        try:
            old.close()
        except Exception:
            pass

    return ds


def close_worker_cache() -> None:
    global _W_DATASET_CACHE

    if _W_DATASET_CACHE is None:
        return

    for ds in _W_DATASET_CACHE.values():
        try:
            ds.close()
        except Exception:
            pass

    _W_DATASET_CACHE.clear()


def worker_initializer(
    file_index: dict[str, dict[int, dict[int, Path]]],
    selected_lat_idx: np.ndarray,
    selected_lon_idx: np.ndarray,
    selected_lat: np.ndarray,
    selected_lon: np.ndarray,
    engine_hash: str,
) -> None:
    global _W_FILE_INDEX, _W_DATASET_CACHE
    global _W_SELECTED_LAT_IDX, _W_SELECTED_LON_IDX
    global _W_SELECTED_LAT, _W_SELECTED_LON
    global _W_REFERENCE_ENGINE_HASH

    ensure_env_for_worker()

    _W_FILE_INDEX = file_index
    _W_DATASET_CACHE = OrderedDict()

    _W_SELECTED_LAT_IDX = np.asarray(selected_lat_idx, dtype=np.int32)
    _W_SELECTED_LON_IDX = np.asarray(selected_lon_idx, dtype=np.int32)
    _W_SELECTED_LAT = np.asarray(selected_lat, dtype=np.float32)
    _W_SELECTED_LON = np.asarray(selected_lon, dtype=np.float32)
    _W_REFERENCE_ENGINE_HASH = engine_hash


# =============================================================================
# Physics
# =============================================================================

def saturation_vapor_pressure(temp_c: np.ndarray) -> np.ndarray:
    temp = np.asarray(temp_c, dtype=np.float32)
    out = np.full(temp.shape, np.nan, dtype=np.float32)

    water = np.isfinite(temp) & (temp >= 0.0)
    ice = np.isfinite(temp) & ~water

    if np.any(water):
        t = temp[water].astype(np.float64)
        out[water] = (
            6.112 * np.exp((17.67 * t) / (t + 243.5))
        ).astype(np.float32)

    if np.any(ice):
        t = temp[ice].astype(np.float64)
        out[ice] = (
            6.112 * np.exp((22.46 * t) / (t + 272.62))
        ).astype(np.float32)

    return out


def convert_temperature(x: np.ndarray, units: Any) -> np.ndarray:
    u = str(units or "").strip().lower().replace("°", "deg")
    a = np.asarray(x, dtype=np.float32)

    if u in {"k", "kelvin"}:
        return a - np.float32(273.15)
    if u in {"c", "celsius", "degc", "degree_celsius", "degrees_celsius"}:
        return a

    raise RuntimeError(f"Unsupported temperature units: {units!r}")


def convert_pressure(x: np.ndarray, units: Any) -> np.ndarray:
    u = str(units or "").strip().lower().replace("°", "deg")
    a = np.asarray(x, dtype=np.float32)

    if u in {"pa", "pascal", "pascals"}:
        return a / np.float32(100.0)
    if u in {"hpa", "mb", "millibar", "millibars"}:
        return a

    raise RuntimeError(f"Unsupported pressure units: {units!r}")


def derive_rh_q(
    T: np.ndarray,
    Td: np.ndarray,
    P: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:

    es_t = saturation_vapor_pressure(T)
    e = saturation_vapor_pressure(Td)

    with np.errstate(divide="ignore", invalid="ignore"):
        rh = np.clip(100.0 * e / es_t, 0.0, 100.0).astype(np.float32)

    valid = (
        np.isfinite(e)
        & np.isfinite(P)
        & (e > 0.0)
        & (P > 0.0)
        & (e < P)
    )

    r = np.full(e.shape, np.nan, dtype=np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        r[valid] = (
            0.622 * e[valid] / (P[valid] - e[valid])
        ).astype(np.float32)

    q = (r / (1.0 + r)).astype(np.float32)
    rh[~np.isfinite(rh)] = np.nan

    return rh, q


# =============================================================================
# Batched window reading
# =============================================================================

def month_span_for_window(target: date) -> list[tuple[int, int]]:
    start = target - timedelta(days=HALF_WIDTH_DAYS)
    end = target + timedelta(days=HALF_WIDTH_DAYS)

    months: list[tuple[int, int]] = []
    cur = date(start.year, start.month, 1)
    end_m = date(end.year, end.month, 1)

    while cur <= end_m:
        months.append((cur.year, cur.month))
        if cur.month == 12:
            cur = date(cur.year + 1, 1, 1)
        else:
            cur = date(cur.year, cur.month + 1, 1)

    return months


def variable_name(ds: xr.Dataset, kind: str) -> str:
    candidates = {
        "t2m": ("t2m",),
        "d2m": ("d2m",),
        "sp": ("sp",),
    }[kind]

    for v in candidates:
        if v in ds.data_vars:
            return v

    raise KeyError(
        f"Could not find {kind} variable in dataset: {list(ds.data_vars)}"
    )


def read_window_block(
    kind: str,
    target: date,
    lat_slice: slice,
    lon_slice: slice,
) -> tuple[np.ndarray, np.ndarray, Any]:

    start = np.datetime64(
        datetime(target.year, target.month, target.day)
        - timedelta(days=HALF_WIDTH_DAYS)
    )
    end = np.datetime64(
        datetime(target.year, target.month, target.day)
        + timedelta(days=HALF_WIDTH_DAYS, hours=23)
    )

    arrays: list[np.ndarray] = []
    times: list[np.ndarray] = []
    units = None

    for year, month in month_span_for_window(target):
        ds = open_cached_dataset(kind, year, month)

        tname = "time" if "time" in ds.coords else "valid_time"
        var = variable_name(ds, kind)

        tt = np.asarray(ds[tname].values)
        sel = np.flatnonzero((tt >= start) & (tt <= end))

        if sel.size == 0:
            continue

        arr = ds[var].isel({
            tname: sel,
            "latitude": lat_slice,
            "longitude": lon_slice,
        }).values.astype(np.float32, copy=False)

        arrays.append(np.asarray(arr))
        times.append(tt[sel])
        units = ds[var].attrs.get("units", units)

    if not arrays:
        return (
            np.empty(
                (
                    0,
                    lat_slice.stop - lat_slice.start,
                    lon_slice.stop - lon_slice.start,
                ),
                dtype=np.float32,
            ),
            np.empty(0, dtype="datetime64[ns]"),
            units,
        )

    tt = np.concatenate(times)
    aa = np.concatenate(arrays, axis=0)
    order = np.argsort(tt)

    return aa[order], tt[order], units


def align_three_variables(
    T: np.ndarray,
    tt: np.ndarray,
    D: np.ndarray,
    td: np.ndarray,
    P: np.ndarray,
    tp: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:

    if not (tt.size and td.size and tp.size):
        shape = (0,) + T.shape[1:]
        return (
            np.empty(shape, np.float32),
            np.empty(shape, np.float32),
            np.empty(shape, np.float32),
        )

    common = np.intersect1d(np.intersect1d(tt, td), tp)

    if common.size == 0:
        shape = (0,) + T.shape[1:]
        return (
            np.empty(shape, np.float32),
            np.empty(shape, np.float32),
            np.empty(shape, np.float32),
        )

    it = np.searchsorted(tt, common)
    id_ = np.searchsorted(td, common)
    ip = np.searchsorted(tp, common)

    return T[it], D[id_], P[ip]


def extract_target_window_block(
    target: date,
    j0: int,
    j1: int,
    i0: int,
    i1: int,
) -> tuple[np.ndarray, np.ndarray]:

    assert _W_SELECTED_LAT_IDX is not None
    assert _W_SELECTED_LON_IDX is not None

    lat_slice = slice(
        int(_W_SELECTED_LAT_IDX[j0]),
        int(_W_SELECTED_LAT_IDX[j1 - 1]) + 1,
    )
    lon_slice = slice(
        int(_W_SELECTED_LON_IDX[i0]),
        int(_W_SELECTED_LON_IDX[i1 - 1]) + 1,
    )

    t, tt, tu = read_window_block("t2m", target, lat_slice, lon_slice)
    d, td, du = read_window_block("d2m", target, lat_slice, lon_slice)
    p, tp, pu = read_window_block("sp", target, lat_slice, lon_slice)

    t, d, p = align_three_variables(t, tt, d, td, p, tp)

    lat_indices = np.asarray(_W_SELECTED_LAT_IDX[j0:j1])
    lon_indices = np.asarray(_W_SELECTED_LON_IDX[i0:i1])

    if lat_indices.size and not np.all(np.diff(lat_indices) == 1):
        all_lat = np.arange(lat_indices[0], lat_indices[-1] + 1)
        take_lat = np.searchsorted(all_lat, lat_indices)
        t = t[:, take_lat, :]
        d = d[:, take_lat, :]
        p = p[:, take_lat, :]

    if lon_indices.size and not np.all(np.diff(lon_indices) == 1):
        all_lon = np.arange(lon_indices[0], lon_indices[-1] + 1)
        take_lon = np.searchsorted(all_lon, lon_indices)
        t = t[:, :, take_lon]
        d = d[:, :, take_lon]
        p = p[:, :, take_lon]

    T = convert_temperature(t, tu)
    Td = convert_temperature(d, du)
    P = convert_pressure(p, pu)

    rh, q = derive_rh_q(T, Td, P)

    return rh.reshape(rh.shape[0], -1), q.reshape(q.shape[0], -1)


# =============================================================================
# Fitting
# =============================================================================

def fit_one_cell(
    rh: np.ndarray,
    q: np.ndarray,
) -> tuple[str, np.ndarray, str, np.ndarray, float, int, str]:

    mask = np.isfinite(rh) & np.isfinite(q)
    n_valid = int(mask.sum())

    if n_valid < FIT_MIN_OBS:
        return (
            "",
            np.full(MAX_PARAM_COUNT, np.nan, np.float32),
            "",
            np.full(MAX_PARAM_COUNT, np.nan, np.float32),
            np.nan,
            n_valid,
            "insufficient_data",
        )

    rhv = rh[mask].astype(np.float64, copy=False)
    qv = q[mask].astype(np.float64, copy=False)

    try:
        result = fit_window_models({"rh": rhv, "q": qv})

        rh_best = result.get("marginals", {}).get("rh", {}).get("best")
        q_best = result.get("marginals", {}).get("q", {}).get("best")
        cop = result.get("copula_rh_q", {})

        rh_name = str(rh_best.get("name", "")) if rh_best else ""
        q_name = str(q_best.get("name", "")) if q_best else ""

        rh_params = params_to_fixed_array(rh_best)
        q_params = params_to_fixed_array(q_best)

        rho = (
            safe_float(cop.get("rho"))
            if cop.get("status") == "ok"
            else np.nan
        )

        status = (
            "ok"
            if rh_best is not None and q_best is not None
            else "fit_failed"
        )

        return (
            rh_name,
            rh_params,
            q_name,
            q_params,
            rho,
            n_valid,
            status,
        )

    except Exception as exc:
        return (
            "",
            np.full(MAX_PARAM_COUNT, np.nan, np.float32),
            "",
            np.full(MAX_PARAM_COUNT, np.nan, np.float32),
            np.nan,
            n_valid,
            f"exception:{type(exc).__name__}:{exc}",
        )


# =============================================================================
# Tile processing
# =============================================================================

def process_tile_task(task: dict[str, Any]) -> dict[str, Any]:
    """
    Process one tile over a half-open DOY interval [doy_start, doy_end).

    Memory optimization:
    one DOY is accumulated and fitted before the next DOY is loaded.
    """

    tile_id = int(task["tile_id"])
    doy_start = int(task["doy_start"])
    doy_end = int(task["doy_end"])

    j0 = int(task["j0"])
    j1 = int(task["j1"])
    i0 = int(task["i0"])
    i1 = int(task["i1"])

    out_path = Path(task["out_path"])
    years = list(task["years"])

    tile_lat = j1 - j0
    tile_lon = i1 - i0
    nday = doy_end - doy_start

    if nday <= 0:
        raise ValueError("Empty DOY task.")

    if 59 >= doy_start and 59 < doy_end:
        raise ValueError("Production task must never contain reserved DOY 59.")

    # Preallocate only final results, not all raw observations.
    rh_names = np.full((nday, tile_lat, tile_lon), "", dtype=object)
    q_names = np.full((nday, tile_lat, tile_lon), "", dtype=object)
    rh_params = np.full(
        (nday, tile_lat, tile_lon, MAX_PARAM_COUNT),
        np.nan,
        np.float32,
    )
    q_params = np.full(
        (nday, tile_lat, tile_lon, MAX_PARAM_COUNT),
        np.nan,
        np.float32,
    )
    rho = np.full((nday, tile_lat, tile_lon), np.nan, np.float32)
    n_obs = np.zeros((nday, tile_lat, tile_lon), np.int32)
    status = np.zeros((nday, tile_lat, tile_lon), np.int8)

    errors = 0
    fits = 0
    t0 = time.perf_counter()

    for local_doy, slot in enumerate(range(doy_start, doy_end)):
        target_dates = target_dates_for_slot(slot, years)

        per_cell_rh: list[list[np.ndarray]] = [[] for _ in range(tile_lat * tile_lon)]
        per_cell_q: list[list[np.ndarray]] = [[] for _ in range(tile_lat * tile_lon)]

        for target in target_dates:
            try:
                rh_block, q_block = extract_target_window_block(
                    target,
                    j0,
                    j1,
                    i0,
                    i1,
                )
            except Exception:
                errors += tile_lat * tile_lon
                continue

            if rh_block.shape[0] == 0:
                continue

            n_cells = tile_lat * tile_lon
            for c in range(n_cells):
                per_cell_rh[c].append(rh_block[:, c])
                per_cell_q[c].append(q_block[:, c])

        for c in range(tile_lat * tile_lon):
            if per_cell_rh[c]:
                rr = np.concatenate(per_cell_rh[c])
                qq = np.concatenate(per_cell_q[c])
            else:
                rr = np.empty(0, np.float32)
                qq = np.empty(0, np.float32)

            if rr.size == 0 or qq.size == 0:
                continue

            rname, rp, qname, qp, r, n, st = fit_one_cell(rr, qq)

            jj, ii = divmod(c, tile_lon)

            rh_names[local_doy, jj, ii] = rname
            q_names[local_doy, jj, ii] = qname
            rh_params[local_doy, jj, ii] = rp
            q_params[local_doy, jj, ii] = qp
            rho[local_doy, jj, ii] = (
                np.float32(r) if np.isfinite(r) else np.nan
            )
            n_obs[local_doy, jj, ii] = n
            status[local_doy, jj, ii] = 1 if st == "ok" else 0

            fits += 1

        # Free this DOY immediately.
        del per_cell_rh, per_cell_q
        gc.collect()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    tmp_path.unlink(missing_ok=True)

    with Dataset(tmp_path, "w", format="NETCDF4") as ds:
        ds.createDimension("doy", nday)
        ds.createDimension("latitude", tile_lat)
        ds.createDimension("longitude", tile_lon)
        ds.createDimension("param", MAX_PARAM_COUNT)

        ds.createVariable(
            "doy", "i2", ("doy",)
        )[:] = np.arange(doy_start, doy_end, dtype=np.int16)

        ds.createVariable(
            "rh_dist_name", str,
            ("doy", "latitude", "longitude")
        )[:] = rh_names

        ds.createVariable(
            "q_dist_name", str,
            ("doy", "latitude", "longitude")
        )[:] = q_names

        ds.createVariable(
            "rh_params", "f4",
            ("doy", "latitude", "longitude", "param"),
            zlib=True,
            complevel=NETCDF_LEVEL,
            shuffle=NETCDF_SHUFFLE,
            fill_value=np.nan,
        )[:] = rh_params

        ds.createVariable(
            "q_params", "f4",
            ("doy", "latitude", "longitude", "param"),
            zlib=True,
            complevel=NETCDF_LEVEL,
            shuffle=NETCDF_SHUFFLE,
            fill_value=np.nan,
        )[:] = q_params

        ds.createVariable(
            "gaussian_copula_rho", "f4",
            ("doy", "latitude", "longitude"),
            zlib=True,
            complevel=NETCDF_LEVEL,
            shuffle=NETCDF_SHUFFLE,
            fill_value=np.nan,
        )[:] = rho

        ds.createVariable(
            "n_obs", "i4",
            ("doy", "latitude", "longitude"),
            zlib=True,
            complevel=NETCDF_LEVEL,
            shuffle=NETCDF_SHUFFLE,
            fill_value=0,
        )[:] = n_obs

        ds.createVariable(
            "status", "i1",
            ("doy", "latitude", "longitude"),
            zlib=True,
            complevel=NETCDF_LEVEL,
            shuffle=NETCDF_SHUFFLE,
            fill_value=0,
        )[:] = status

        ds.setncattr("schema_version", SCHEMA_VERSION)
        ds.setncattr("engine_config_hash", _W_REFERENCE_ENGINE_HASH or "unknown")
        ds.setncattr("doy_start", doy_start)
        ds.setncattr("doy_end_inclusive", doy_end - 1)
        ds.setncattr("tile_id", tile_id)
        ds.setncattr("fit_min_obs", FIT_MIN_OBS)
        ds.setncattr(
            "calendar_policy",
            "366-slot engine policy; slot 59 reserved; "
            "slot 60 pools Feb-28/Feb-29; slot 61=Mar-01",
        )
        ds.setncattr("hours_per_target_window", WINDOW_HOURS)
        ds.sync()

    os.replace(tmp_path, out_path)

    elapsed = time.perf_counter() - t0

    return {
        "tile_id": tile_id,
        "out_path": str(out_path),
        "fits": int(fits),
        "errors": int(errors),
        "seconds": float(elapsed),
        "fit_per_sec": float(fits / max(elapsed, 1e-9)),
    }


# =============================================================================
# Checkpoint / tasks
# =============================================================================

def validate_tile_file(
    path: Path,
    doy_start: int,
    doy_end: int,
    lat_len: int,
    lon_len: int,
) -> bool:

    if not path.exists():
        return False

    try:
        with Dataset(path, "r") as ds:
            if "doy" not in ds.dimensions:
                return False
            if ds.dimensions["doy"].size != doy_end - doy_start:
                return False
            if ds.dimensions["latitude"].size != lat_len:
                return False
            if ds.dimensions["longitude"].size != lon_len:
                return False
            if "status" not in ds.variables:
                return False
            if ds.variables["status"].shape != (
                doy_end - doy_start,
                lat_len,
                lon_len,
            ):
                return False
        return True
    except Exception:
        return False


def run_root() -> Path:
    tag = (
        f"doy_{DOY_START:03d}_{DOY_END:03d}"
        f"__{CONFIG_RUN_HASH[:12]}"
    )
    return OUTPUT_DIR / "moisture_copula_runs_v82" / tag


def manifest_path(root: Path) -> Path:
    return root / "manifest.json"


def tile_path(
    root: Path,
    tile_id: int,
    doy_start: int,
    doy_end: int,
) -> Path:
    return (
        root
        / "tiles"
        / f"tile_{tile_id:05d}_doy_{doy_start:03d}_{doy_end - 1:03d}.nc"
    )


def done_marker(
    root: Path,
    tile_id: int,
    doy_start: int,
    doy_end: int,
) -> Path:
    return (
        root
        / "done"
        / f"tile_{tile_id:05d}_doy_{doy_start:03d}_{doy_end - 1:03d}.json"
    )


def make_tasks(
    root: Path,
    selected_lat_idx: np.ndarray,
    selected_lon_idx: np.ndarray,
    years: list[int],
) -> list[dict[str, Any]]:

    tasks: list[dict[str, Any]] = []
    tile_id = 0

    doy_ranges = chunk_ranges_excluding_reserved(
        DOY_START,
        DOY_END,
        DOY_CHUNK,
    )

    for j0 in range(0, len(selected_lat_idx), TILE_LAT):
        j1 = min(j0 + TILE_LAT, len(selected_lat_idx))

        for i0 in range(0, len(selected_lon_idx), TILE_LON):
            i1 = min(i0 + TILE_LON, len(selected_lon_idx))

            for d0, d1 in doy_ranges:
                out = tile_path(root, tile_id, d0, d1)

                tasks.append({
                    "tile_id": tile_id,
                    "doy_start": d0,
                    "doy_end": d1,
                    "j0": j0,
                    "j1": j1,
                    "i0": i0,
                    "i1": i1,
                    "out_path": str(out),
                    "years": years,
                })

                tile_id += 1

    return tasks


def init_manifest(
    root: Path,
    tasks: list[dict[str, Any]],
    selected_lat: np.ndarray,
    selected_lon: np.ndarray,
    years: list[int],
) -> None:

    root.mkdir(parents=True, exist_ok=True)
    (root / "tiles").mkdir(exist_ok=True)
    (root / "done").mkdir(exist_ok=True)

    payload = {
        "status": "running",
        "schema_version": SCHEMA_VERSION,
        "engine_config_hash": ENGINE_CONFIG_HASH,
        "run_config_hash": CONFIG_RUN_HASH,
        "doy_start": DOY_START,
        "doy_end": DOY_END,
        "years": years,
        "workers": MAX_WORKERS,
        "max_in_flight": MAX_IN_FLIGHT,
        "doy_chunk": DOY_CHUNK,
        "tile_lat": TILE_LAT,
        "tile_lon": TILE_LON,
        "dataset_cache_size": DATASET_CACHE_SIZE,
        "grid_shape": [len(selected_lat), len(selected_lon)],
        "task_count": len(tasks),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }

    p = manifest_path(root)
    if not p.exists():
        atomic_json_write(p, payload)


def mark_done(
    root: Path,
    result: dict[str, Any],
    task: dict[str, Any],
) -> None:

    atomic_json_write(
        done_marker(
            root,
            result["tile_id"],
            task["doy_start"],
            task["doy_end"],
        ),
        {
            "status": "completed",
            "tile_id": result["tile_id"],
            "doy_start": task["doy_start"],
            "doy_end_inclusive": task["doy_end"] - 1,
            "out_path": result["out_path"],
            "fits": result["fits"],
            "errors": result["errors"],
            "seconds": result["seconds"],
            "fit_per_sec": result["fit_per_sec"],
            "sha256": sha256_file(Path(result["out_path"])),
            "completed_utc": datetime.now(timezone.utc).isoformat(),
        },
    )


def task_is_done(root: Path, task: dict[str, Any]) -> bool:
    marker = done_marker(
        root,
        task["tile_id"],
        task["doy_start"],
        task["doy_end"],
    )

    path = Path(task["out_path"])

    return (
        marker.exists()
        and validate_tile_file(
            path,
            task["doy_start"],
            task["doy_end"],
            task["j1"] - task["j0"],
            task["i1"] - task["i0"],
        )
    )


# =============================================================================
# Benchmark
# =============================================================================

def make_benchmark_tasks(
    selected_lat_idx: np.ndarray,
    selected_lon_idx: np.ndarray,
    years: list[int],
) -> list[dict[str, Any]]:

    nlat = min(BENCHMARK_LAT, len(selected_lat_idx))
    nlon = min(BENCHMARK_LON, len(selected_lon_idx))

    if nlat <= 0 or nlon <= 0:
        raise RuntimeError("Benchmark spatial subset is empty.")

    # Use one realistic production-sized spatial tile by default.
    # To exercise several processes, create multiple disjoint/shifted spatial
    # tasks. They are benchmark-only and produce temporary outputs.
    spatial_starts: list[tuple[int, int]] = []

    lat_step = max(1, nlat)
    lon_step = max(1, nlon)

    for k in range(BENCHMARK_TASKS):
        max_j0 = max(0, len(selected_lat_idx) - nlat)
        max_i0 = max(0, len(selected_lon_idx) - nlon)

        j0 = min((k * lat_step) % max(max_j0 + 1, 1), max_j0)
        i0 = min((k * lon_step) % max(max_i0 + 1, 1), max_i0)

        spatial_starts.append((j0, i0))

    valid_doys = [s for s in valid_slots(DOY_START, DOY_END)]
    if not valid_doys:
        raise RuntimeError("Benchmark DOY interval contains no valid slots.")

    d0 = valid_doys[0]
    selected_doys = valid_doys[:BENCHMARK_DAYS]

    # Build a contiguous range only when it does not cross reserved slot 59.
    # For benchmark simplicity, if 59 would be crossed, use a valid range after it.
    if selected_doys[-1] == d0 + len(selected_doys) - 1:
        ranges = [(d0, selected_doys[-1] + 1)]
    else:
        d0 = selected_doys[0]
        ranges = [(d0, d0 + 1)]

    tasks: list[dict[str, Any]] = []

    for task_id, (j0, i0) in enumerate(spatial_starts):
        # For normal default DOY_START=1 this is e.g. 1..4 for 3 DOYs.
        bd0, bd1 = ranges[0]

        out = (
            OUTPUT_DIR
            / ".benchmark_v82"
            / f"task_{task_id:04d}.nc"
        )

        tasks.append({
            "tile_id": task_id,
            "doy_start": bd0,
            "doy_end": bd1,
            "j0": j0,
            "j1": j0 + nlat,
            "i0": i0,
            "i1": i0 + nlon,
            "out_path": str(out),
            "years": years,
        })

    return tasks


def run_benchmark(
    selected_lat_idx: np.ndarray,
    selected_lon_idx: np.ndarray,
    selected_lat: np.ndarray,
    selected_lon: np.ndarray,
    years: list[int],
) -> None:

    ensure_env_for_worker()
    file_index = build_all_file_indices(years)

    tmp_root = OUTPUT_DIR / ".benchmark_v82"
    shutil.rmtree(tmp_root, ignore_errors=True)
    tmp_root.mkdir(parents=True, exist_ok=True)

    tasks = make_benchmark_tasks(
        selected_lat_idx,
        selected_lon_idx,
        years,
    )

    # For a reliable benchmark, use exactly the same tile processing path as
    # production, and run multiple independent tile tasks concurrently.
    workers = min(BENCHMARK_WORKERS, len(tasks))

    logger.info(
        "BENCHMARK | tasks=%d | workers=%d | tile=%dx%d | DOYs=%d | years=%d",
        len(tasks),
        workers,
        BENCHMARK_LAT,
        BENCHMARK_LON,
        BENCHMARK_DAYS,
        len(years),
    )

    started = time.perf_counter()
    completed = 0
    total_fits = 0
    total_errors = 0
    task_seconds: list[float] = []

    ctx = __import__("multiprocessing").get_context("spawn")

    try:
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=ctx,
            initializer=worker_initializer,
            initargs=(
                file_index,
                selected_lat_idx,
                selected_lon_idx,
                selected_lat,
                selected_lon,
                ENGINE_CONFIG_HASH,
            ),
        ) as executor:

            futures = {
                executor.submit(process_tile_task, task): task
                for task in tasks
            }

            while futures:
                done, _ = wait(
                    futures,
                    return_when=FIRST_COMPLETED,
                )

                for future in done:
                    task = futures.pop(future)
                    result = future.result()

                    completed += 1
                    total_fits += int(result["fits"])
                    total_errors += int(result["errors"])
                    task_seconds.append(float(result["seconds"]))

                    logger.info(
                        "BENCHMARK PROGRESS | %d/%d | task=%d | fits=%d | "
                        "sec=%.2f | fit/sec=%.3f",
                        completed,
                        len(tasks),
                        result["tile_id"],
                        result["fits"],
                        result["seconds"],
                        result["fit_per_sec"],
                    )

    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)

    elapsed = time.perf_counter() - started
    throughput = total_fits / max(elapsed, 1e-9)
    task_rate = completed / max(elapsed, 1e-9)

    # A rough per-worker scaling indicator. It is not claimed to be ideal
    # linear scaling; it is simply observed production-path throughput.
    per_worker_throughput = throughput / max(workers, 1)
    theoretical_single_fit_rate = throughput / max(workers, 1)
    efficiency = 100.0 * throughput / max(theoretical_single_fit_rate * workers, 1e-12)

    total_production_fits = (
        len(valid_slots(DOY_START, DOY_END))
        * len(selected_lat_idx)
        * len(selected_lon_idx)
    )
    eta_seconds = total_production_fits / max(throughput, 1e-9)

    logger.info(
        "BENCHMARK RESULT | fits=%d | elapsed=%.2fs | fit/sec=%.3f | "
        "task/sec=%.3f | avg_task_sec=%.2f | errors=%d",
        total_fits,
        elapsed,
        throughput,
        task_rate,
        float(np.mean(task_seconds)) if task_seconds else float("nan"),
        total_errors,
    )

    logger.info(
        "BENCHMARK ESTIMATE | observed_workers=%d | per-worker_observed=%.3f "
        "fit/sec | full_valid_fits=%d | projected_full_run=%.2f h (%.2f days)",
        workers,
        per_worker_throughput,
        total_production_fits,
        eta_seconds / 3600.0,
        eta_seconds / 86400.0,
    )

    # Real elapsed scaling efficiency is intentionally not computed against an
    # assumed single-worker speed. A separate --benchmark-scale mode is used
    # when exact 1/2/4-worker scaling is desired.


# =============================================================================
# Production parallel execution
# =============================================================================

def execute_tasks(
    root: Path,
    tasks: list[dict[str, Any]],
    selected_lat: np.ndarray,
    selected_lon: np.ndarray,
    selected_lat_idx: np.ndarray,
    selected_lon_idx: np.ndarray,
    years: list[int],
) -> None:

    pending_tasks = [
        t for t in tasks if not task_is_done(root, t)
    ]

    skipped = len(tasks) - len(pending_tasks)

    logger.info(
        "RESUME | total tasks=%d | already complete=%d | remaining=%d",
        len(tasks),
        skipped,
        len(pending_tasks),
    )

    if not pending_tasks:
        return

    file_index = build_all_file_indices(years)

    ctx = __import__("multiprocessing").get_context("spawn")
    started = time.perf_counter()

    completed = skipped
    total = len(tasks)
    total_fits = 0
    total_errors = 0

    pending: dict[Any, dict[str, Any]] = {}
    next_index = 0

    with ProcessPoolExecutor(
        max_workers=MAX_WORKERS,
        mp_context=ctx,
        initializer=worker_initializer,
        initargs=(
            file_index,
            selected_lat_idx,
            selected_lon_idx,
            selected_lat,
            selected_lon,
            ENGINE_CONFIG_HASH,
        ),
    ) as executor:

        while (
            next_index < len(pending_tasks)
            and len(pending) < MAX_IN_FLIGHT
        ):
            task = pending_tasks[next_index]
            pending[executor.submit(process_tile_task, task)] = task
            next_index += 1

        while pending:
            done, _ = wait(
                pending,
                return_when=FIRST_COMPLETED,
            )

            for future in done:
                task = pending.pop(future)
                result = future.result()

                if not validate_tile_file(
                    Path(result["out_path"]),
                    task["doy_start"],
                    task["doy_end"],
                    task["j1"] - task["j0"],
                    task["i1"] - task["i0"],
                ):
                    raise RuntimeError(
                        f"Invalid tile output: {result['out_path']}"
                    )

                mark_done(root, result, task)

                completed += 1
                total_fits += int(result["fits"])
                total_errors += int(result["errors"])

                elapsed = time.perf_counter() - started
                task_rate = completed / max(elapsed, 1e-9)
                remaining = total - completed
                eta = remaining / max(task_rate, 1e-9)

                logger.info(
                    "PROGRESS | %d/%d | %.2f%% | tile=%d | fits=%d | "
                    "errors=%d | rate=%.3f tasks/s | ETA=%.2fh",
                    completed,
                    total,
                    100.0 * completed / max(total, 1),
                    result["tile_id"],
                    result["fits"],
                    result["errors"],
                    task_rate,
                    eta / 3600.0,
                )

                if next_index < len(pending_tasks):
                    task2 = pending_tasks[next_index]
                    pending[executor.submit(process_tile_task, task2)] = task2
                    next_index += 1

    manifest = json.loads(
        manifest_path(root).read_text(encoding="utf-8")
    )

    atomic_json_write(
        manifest_path(root),
        {
            **manifest,
            "status": "completed",
            "completed_tasks": total,
            "total_tasks": total,
            "fit_results_reported_this_run": total_fits,
            "cell_errors_reported_this_run": total_errors,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
        },
    )


# =============================================================================
# NetCDF merge
# =============================================================================

def create_console_output(
    path: Path,
    selected_lat: np.ndarray,
    selected_lon: np.ndarray,
) -> Dataset:

    if path.exists():
        if OVERWRITE_CONSOLE_OUTPUT:
            path.unlink()
        else:
            raise FileExistsError(
                f"Console output exists: {path}. "
                "Set OVERWRITE_CONSOLE_OUTPUT=1 to replace it."
            )

    ds = Dataset(path, "w", format="NETCDF4")

    ds.createDimension("doy", len(valid_slots(DOY_START, DOY_END)))
    ds.createDimension("latitude", len(selected_lat))
    ds.createDimension("longitude", len(selected_lon))
    ds.createDimension("param", MAX_PARAM_COUNT)

    ds.createVariable("doy", "i2", ("doy",))[:] = np.asarray(
        valid_slots(DOY_START, DOY_END),
        dtype=np.int16,
    )

    ds.createVariable(
        "latitude", "f4", ("latitude",)
    )[:] = selected_lat.astype(np.float32)

    ds.createVariable(
        "longitude", "f4", ("longitude",)
    )[:] = selected_lon.astype(np.float32)

    chunks3 = (
        1,
        min(TILE_LAT, len(selected_lat)),
        min(TILE_LON, len(selected_lon)),
    )

    for name, dtype, dims, fill in (
        (
            "rh_dist_name",
            str,
            ("doy", "latitude", "longitude"),
            None,
        ),
        (
            "q_dist_name",
            str,
            ("doy", "latitude", "longitude"),
            None,
        ),
        (
            "gaussian_copula_rho",
            "f4",
            ("doy", "latitude", "longitude"),
            np.nan,
        ),
        (
            "n_obs",
            "i4",
            ("doy", "latitude", "longitude"),
            0,
        ),
        (
            "status",
            "i1",
            ("doy", "latitude", "longitude"),
            0,
        ),
    ):
        kwargs = {}

        if dtype is not str:
            kwargs = {
                "zlib": True,
                "complevel": NETCDF_LEVEL,
                "shuffle": NETCDF_SHUFFLE,
                "chunksizes": chunks3,
                "fill_value": fill,
            }

        ds.createVariable(name, dtype, dims, **kwargs)

    for name in ("rh_params", "q_params"):
        ds.createVariable(
            name,
            "f4",
            ("doy", "latitude", "longitude", "param"),
            zlib=True,
            complevel=NETCDF_LEVEL,
            shuffle=NETCDF_SHUFFLE,
            chunksizes=(
                1,
                min(TILE_LAT, len(selected_lat)),
                min(TILE_LON, len(selected_lon)),
                MAX_PARAM_COUNT,
            ),
            fill_value=np.nan,
        )

    ds.schema_version = SCHEMA_VERSION
    ds.engine_config_hash = ENGINE_CONFIG_HASH
    ds.run_doy_start = DOY_START
    ds.run_doy_end = DOY_END
    ds.calendar_policy = (
        "366-slot engine policy; slot 59 reserved; "
        "slot 60 pools Feb-28/Feb-29; slot 61=Mar-01"
    )

    ds.parameter_mapping = json.dumps({
        "Normal": ["loc", "scale"],
        "SkewNormal": ["shape", "loc", "scale"],
        "PearsonIII": ["skew", "loc", "scale"],
        "Beta": ["a", "b", "loc", "scale"],
        "BimodalNormal": ["w1", "mu1", "sigma1", "mu2", "sigma2"],
    })

    return ds


def merge_tiles(
    root: Path,
    tasks: list[dict[str, Any]],
    selected_lat: np.ndarray,
    selected_lon: np.ndarray,
    output_path: Path,
) -> None:

    missing = [t for t in tasks if not task_is_done(root, t)]

    if missing:
        raise RuntimeError(
            f"Cannot merge: {len(missing)} tile tasks are incomplete. "
            f"First missing tile={missing[0]['tile_id']}"
        )

    tmp = output_path.with_suffix(output_path.suffix + ".tmp")
    tmp.unlink(missing_ok=True)

    dst = create_console_output(
        tmp,
        selected_lat,
        selected_lon,
    )

    try:
        valid_doy_list = valid_slots(DOY_START, DOY_END)
        doy_to_output = {d: i for i, d in enumerate(valid_doy_list)}

        for task in sorted(
            tasks,
            key=lambda t: (t["doy_start"], t["j0"], t["i0"]),
        ):
            tile = Path(task["out_path"])

            with Dataset(tile, "r") as src:
                tile_doys = np.asarray(
                    src.variables["doy"][:],
                    dtype=int,
                )

                # Destination indexes are compact because slot 59 is absent.
                dest_idx = [doy_to_output[d] for d in tile_doys]
                a = min(dest_idx)
                b = max(dest_idx) + 1

                j0, j1 = task["j0"], task["j1"]
                i0, i1 = task["i0"], task["i1"]

                # Production task ranges never cross reserved slot 59, so the
                # tile interval is contiguous in the compact destination.
                dst.variables["rh_dist_name"][a:b, j0:j1, i0:i1] = (
                    src.variables["rh_dist_name"][:]
                )
                dst.variables["q_dist_name"][a:b, j0:j1, i0:i1] = (
                    src.variables["q_dist_name"][:]
                )
                dst.variables["rh_params"][a:b, j0:j1, i0:i1, :] = (
                    src.variables["rh_params"][:]
                )
                dst.variables["q_params"][a:b, j0:j1, i0:i1, :] = (
                    src.variables["q_params"][:]
                )
                dst.variables["gaussian_copula_rho"][a:b, j0:j1, i0:i1] = (
                    src.variables["gaussian_copula_rho"][:]
                )
                dst.variables["n_obs"][a:b, j0:j1, i0:i1] = (
                    src.variables["n_obs"][:]
                )
                dst.variables["status"][a:b, j0:j1, i0:i1] = (
                    src.variables["status"][:]
                )

        dst.sync()

    finally:
        dst.close()

    os.replace(tmp, output_path)

    logger.info(
        "MERGE COMPLETE | %s | SHA256=%s",
        output_path,
        sha256_file(output_path),
    )


# =============================================================================
# Console runner
# =============================================================================

def run_console(benchmark: bool = False) -> None:
    validate_calendar_contract()
    validate_doy_range(DOY_START, DOY_END)

    lat, lon, lat_idx, lon_idx = configure_grid()

    selected_lat = lat[lat_idx]
    selected_lon = lon[lon_idx]

    years = [
        y for y in PROCESS_YEARS
        if START_YEAR <= y <= END_YEAR
    ]

    if not years:
        raise RuntimeError("No valid processing years.")

    valid_count = len(valid_slots(DOY_START, DOY_END))

    logger.info(
        "GRID | original=%dx%d | selected=%dx%d",
        len(lat),
        len(lon),
        len(selected_lat),
        len(selected_lon),
    )

    logger.info(
        "RUN | DOY=%d..%d | valid_slots=%d | years=%s | workers=%d | "
        "in-flight=%d | tile=%dx%d | DOY chunk=%d | cache=%d",
        DOY_START,
        DOY_END,
        valid_count,
        years,
        MAX_WORKERS,
        MAX_IN_FLIGHT,
        TILE_LAT,
        TILE_LON,
        DOY_CHUNK,
        DATASET_CACHE_SIZE,
    )

    logger.info(
        "TOTAL FITS | %d",
        valid_count * len(selected_lat) * len(selected_lon),
    )

    if benchmark:
        run_benchmark(
            lat_idx,
            lon_idx,
            selected_lat,
            selected_lon,
            years,
        )
        return

    root = run_root()

    tasks = make_tasks(
        root,
        lat_idx,
        lon_idx,
        years,
    )

    init_manifest(
        root,
        tasks,
        selected_lat,
        selected_lon,
        years,
    )

    logger.info(
        "TASK PLAN | %d tasks | reserved DOY 59 excluded",
        len(tasks),
    )

    execute_tasks(
        root,
        tasks,
        selected_lat,
        selected_lon,
        lat_idx,
        lon_idx,
        years,
    )

    output_path = (
        OUTPUT_DIR
        / f"moisture_copula_parameters_v82_"
          f"doy_{DOY_START:03d}_{DOY_END:03d}.nc"
    )

    merge_tiles(
        root,
        tasks,
        selected_lat,
        selected_lon,
        output_path,
    )

    logger.info("DONE | output=%s", output_path)
    logger.info("RUN DIRECTORY | %s", root)


# =============================================================================
# Multi-console merge
# =============================================================================

def create_output_file_for_merge(
    path: Path,
    lat: np.ndarray,
    lon: np.ndarray,
    start: int,
    end: int,
) -> Dataset:

    doys = valid_slots(start, end)

    ds = Dataset(path, "w", format="NETCDF4")

    ds.createDimension("doy", len(doys))
    ds.createDimension("latitude", len(lat))
    ds.createDimension("longitude", len(lon))
    ds.createDimension("param", MAX_PARAM_COUNT)

    ds.createVariable("doy", "i2", ("doy",))[:] = np.asarray(
        doys,
        dtype=np.int16,
    )

    ds.createVariable("latitude", "f4", ("latitude",))[:] = lat.astype(np.float32)
    ds.createVariable("longitude", "f4", ("longitude",))[:] = lon.astype(np.float32)

    ds.createVariable(
        "rh_dist_name",
        str,
        ("doy", "latitude", "longitude"),
    )
    ds.createVariable(
        "q_dist_name",
        str,
        ("doy", "latitude", "longitude"),
    )

    ds.createVariable(
        "rh_params", "f4",
        ("doy", "latitude", "longitude", "param"),
        zlib=True,
        complevel=NETCDF_LEVEL,
        shuffle=NETCDF_SHUFFLE,
        fill_value=np.nan,
    )

    ds.createVariable(
        "q_params", "f4",
        ("doy", "latitude", "longitude", "param"),
        zlib=True,
        complevel=NETCDF_LEVEL,
        shuffle=NETCDF_SHUFFLE,
        fill_value=np.nan,
    )

    ds.createVariable(
        "gaussian_copula_rho", "f4",
        ("doy", "latitude", "longitude"),
        zlib=True,
        complevel=NETCDF_LEVEL,
        shuffle=NETCDF_SHUFFLE,
        fill_value=np.nan,
    )

    ds.createVariable(
        "n_obs", "i4",
        ("doy", "latitude", "longitude"),
        zlib=True,
        complevel=NETCDF_LEVEL,
        shuffle=NETCDF_SHUFFLE,
        fill_value=0,
    )

    ds.createVariable(
        "status", "i1",
        ("doy", "latitude", "longitude"),
        zlib=True,
        complevel=NETCDF_LEVEL,
        shuffle=NETCDF_SHUFFLE,
        fill_value=0,
    )

    ds.schema_version = SCHEMA_VERSION
    ds.engine_config_hash = ENGINE_CONFIG_HASH
    ds.calendar_policy = (
        "366-slot engine policy; slot 59 reserved; "
        "slot 60 pools Feb-28/Feb-29; slot 61=Mar-01"
    )

    return ds


def merge_multiple_console_files(
    paths: list[Path],
    output_path: Path,
) -> None:

    if not paths:
        raise ValueError("No input files.")

    meta: list[tuple[int, int, Path]] = []

    for p in paths:
        with Dataset(p, "r") as ds:
            doys = np.asarray(
                ds.variables["doy"][:],
                dtype=int,
            )

            if doys.size == 0:
                continue

            meta.append(
                (
                    int(doys[0]),
                    int(doys[-1]),
                    p,
                )
            )

    meta.sort()

    for (_, e1, _), (s2, _, _) in zip(meta, meta[1:]):
        if e1 + 1 != s2:
            raise RuntimeError(
                f"DOY intervals are not contiguous: {e1} -> {s2}"
            )

    with Dataset(meta[0][2], "r") as src0:
        lat = np.asarray(src0.variables["latitude"][:])
        lon = np.asarray(src0.variables["longitude"][:])

    global_start = meta[0][0]
    global_end = meta[-1][1]

    if output_path.exists():
        output_path.unlink()

    dst = create_output_file_for_merge(
        output_path,
        lat,
        lon,
        global_start,
        global_end,
    )

    try:
        global_doys = valid_slots(global_start, global_end)
        doy_to_index = {d: i for i, d in enumerate(global_doys)}

        for s, e, p in meta:
            with Dataset(p, "r") as src:
                src_doys = np.asarray(
                    src.variables["doy"][:],
                    dtype=int,
                )

                inds = [doy_to_index[int(d)] for d in src_doys]
                a = min(inds)
                b = max(inds) + 1

                dst.variables["rh_dist_name"][a:b] = (
                    src.variables["rh_dist_name"][:]
                )
                dst.variables["q_dist_name"][a:b] = (
                    src.variables["q_dist_name"][:]
                )
                dst.variables["rh_params"][a:b] = (
                    src.variables["rh_params"][:]
                )
                dst.variables["q_params"][a:b] = (
                    src.variables["q_params"][:]
                )
                dst.variables["gaussian_copula_rho"][a:b] = (
                    src.variables["gaussian_copula_rho"][:]
                )
                dst.variables["n_obs"][a:b] = (
                    src.variables["n_obs"][:]
                )
                dst.variables["status"][a:b] = (
                    src.variables["status"][:]
                )

        dst.sync()

    finally:
        dst.close()

    logger.info(
        "MULTI-CONSOLE MERGE COMPLETE | %s | SHA256=%s",
        output_path,
        sha256_file(output_path),
    )


# =============================================================================
# Config hash
# =============================================================================

CONFIG_RUN = {
    "schema_version": SCHEMA_VERSION,
    "engine_config_hash": ENGINE_CONFIG_HASH,
    "doy_start": DOY_START,
    "doy_end": DOY_END,
    "max_workers": MAX_WORKERS,
    "max_in_flight": MAX_IN_FLIGHT,
    "doy_chunk": DOY_CHUNK,
    "tile_lat": TILE_LAT,
    "tile_lon": TILE_LON,
    "dataset_cache_size": DATASET_CACHE_SIZE,
    "fit_min_obs": FIT_MIN_OBS,
    "lat_range": LAT_RANGE,
    "lon_range": LON_RANGE,
    "stride_lat": STRIDE_LAT,
    "stride_lon": STRIDE_LON,
    "half_width_days": HALF_WIDTH_DAYS,
    "years": PROCESS_YEARS,
}

CONFIG_RUN_HASH = hashlib.sha256(
    json.dumps(
        CONFIG_RUN,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
).hexdigest()[:20]


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Production 5-day moisture copula parameter runner v8.2"
    )

    p.add_argument(
        "--benchmark",
        action="store_true",
        help="Run production-path multiprocessing benchmark only.",
    )

    p.add_argument(
        "--merge",
        nargs="+",
        type=Path,
        help="Merge console output files.",
    )

    p.add_argument(
        "--output",
        type=Path,
        help="Output path for --merge.",
    )

    return p.parse_args()


def main() -> None:
    args = parse_args()
    ensure_env_for_worker()

    if args.merge:
        if args.output is None:
            raise SystemExit("--merge requires --output")
        merge_multiple_console_files(args.merge, args.output)
        return

    run_console(benchmark=args.benchmark)


if __name__ == "__main__":
    import multiprocessing as mp

    mp.freeze_support()

    try:
        main()
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        raise
    except Exception:
        logger.exception("Fatal error.")
        raise
