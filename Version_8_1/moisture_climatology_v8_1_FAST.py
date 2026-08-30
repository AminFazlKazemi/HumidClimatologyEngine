#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ERA5-Land hourly empirical moisture climatology engine v8.0 fast optimization draft.

Scientific contract
-------------------
v8.0 FINAL SINGLE-PASS is a direct hourly empirical production engine. It does NOT use the v6 daily-statistic
Gaussian/Monte-Carlo architecture for its primary moisture products.

For each climatological DOY and grid cell, hourly T2m, D2m and surface pressure
are transformed to RH, vapor pressure e, mixing ratio r and specific humidity q.
The engine retains count, mean, central sums M2/M3/M4 and selected pairwise
covariance states. Year checkpoints are disk-backed NetCDF files written one DOY
at a time, so peak RAM does not scale with the number of climatological days.

The bivariate probability layer is empirical-first. A separate five-day centered window layer is provided for distribution fitting and station queries. The primary joint product is an
empirical 2-D probability-mass function (piecewise-constant PDF) built directly
from hourly paired observations for each climatological DOY and grid cell. No
Gaussian shape is imposed. A bivariate Gaussian evaluator is retained only as a
reference candidate; Beta/coplanula-style parametric models may be compared
later, but they are never silently treated as the truth.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import xarray as xr
from tqdm import tqdm

try:
    from netCDF4 import Dataset
except Exception:  # pragma: no cover - optional for pure scientific/unit tests
    Dataset = None  # type: ignore[assignment]


def require_netcdf4() -> None:
    if Dataset is None:
        raise ImportError(
            "netCDF4 is required for checkpointing and NetCDF finalization. "
            "Install the production dependency before running the full climatology."
        )

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

T2M_DIR = Path(r"F:\Kazemi\era5\land\T2m")
D2M_DIR = Path(r"F:\Kazemi\era5\land\Dew_Point_Temperature")
SP_DIR = Path(r"F:\Kazemi\era5\land\Surface_Pressure")
OUTPUT_DIR = Path(r"C:\c")

CHECKPOINT_DIR = OUTPUT_DIR / "checkpoints_moisture_v8_0"
YEAR_DIR = CHECKPOINT_DIR / "years"
RUN_MANIFEST_FILE = OUTPUT_DIR / "moisture_climatology_run_manifest_v8_0.json"

START_YEAR = 1981
END_YEAR = 2020
DOY_COUNT = 366

MAX_WORKERS = max(1, (os.cpu_count() or 4) - 2)
CHUNK_LAT = 64
CHUNK_LON = 128
PROGRESS_FLUSH_CHUNKS = 16
PROGRESS_LOG_EVERY_CHUNKS = 8
PROGRESS_REFRESH_SECONDS = 5.0

SCHEMA_VERSION = "8.0"
CHECKPOINT_VERSION = "8.0-FINAL-SINGLE-PASS"

# Primary pairwise probability-function parameterizations.
# Additional pairs can be added without changing the primary marginal statistics.
BIVARIATE_PAIRS = (("rh", "q"),)
BIVARIATE_NX = 8
BIVARIATE_NY = 8
BUILD_EMPIRICAL_BIVARIATE = True

# v8.0 optimization contract: cached metadata, larger spatial chunks, automatic workers.
V8_FAST_FEATURES = {"optimized_chunks": True, "automatic_workers": True, "single_pass_bivariate_target": True}
BIVARIATE_CHECKPOINT_FILE = CHECKPOINT_DIR / "bivariate_progress_v8_0.json"

OUTPUT_FILE = OUTPUT_DIR / "moisture_climatology_1981_2020_v8_0.nc"
DIAGNOSTIC_FILE = OUTPUT_DIR / "moisture_climatology_diagnostics_1981_2020_v8_0.nc"
BIVARIATE_FILE = OUTPUT_DIR / "moisture_climatology_bivariate_1981_2020_v8_0.nc"

for _d in (OUTPUT_DIR, CHECKPOINT_DIR, YEAR_DIR):
    _d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("Moisture_Climatology_v8_0")

# =============================================================================
# 2. CONFIG / HASHING / PROVENANCE
# =============================================================================

@dataclass(frozen=True)
class Config:
    start_year: int
    end_year: int
    doy_count: int
    chunk_lat: int
    chunk_lon: int
    workers: int
    progress_flush_chunks: int
    schema_version: str
    bivariate_pairs: tuple[tuple[str, str], ...]
    bivariate_nx: int
    bivariate_ny: int
    build_empirical_bivariate: bool


CONFIG = Config(
    start_year=START_YEAR,
    end_year=END_YEAR,
    doy_count=DOY_COUNT,
    chunk_lat=CHUNK_LAT,
    chunk_lon=CHUNK_LON,
    workers=MAX_WORKERS,
    progress_flush_chunks=PROGRESS_FLUSH_CHUNKS,
    schema_version=SCHEMA_VERSION,
    bivariate_pairs=tuple(BIVARIATE_PAIRS),
    bivariate_nx=BIVARIATE_NX,
    bivariate_ny=BIVARIATE_NY,
    build_empirical_bivariate=BUILD_EMPIRICAL_BIVARIATE,
)


def hash_dict(data: dict) -> str:
    payload = json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


CONFIG_HASH = hash_dict(asdict(CONFIG))


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(block_size), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent), text=True
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


def script_sha256() -> str:
    try:
        return sha256_file(Path(__file__))
    except Exception:
        return "unavailable"

# =============================================================================
# 3. CALENDAR CONTRACT
# =============================================================================

def is_leap_year(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def get_clim_doy(native_doy: int, year: int) -> int:
    """Map Gregorian DOY to the 366-slot project calendar.

    Slot 59 is reserved.
    Slot 60 pools Feb-28 and Feb-29.
    Slot 61 is Mar-01.
    """
    if not 1 <= native_doy <= 366:
        return -1
    if is_leap_year(year):
        if native_doy in (59, 60):
            return 60
        return native_doy
    if native_doy == 59:
        return 60
    if native_doy >= 60:
        return native_doy + 1
    return native_doy


def calendar_labels() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    import datetime as _dt

    month = np.full(DOY_COUNT, -1, dtype=np.int16)
    day = np.full(DOY_COUNT, -1, dtype=np.int16)
    label = np.full(DOY_COUNT, "RESERVED", dtype=object)
    reserved = np.zeros(DOY_COUNT, dtype=np.int8)

    d = _dt.date(2001, 1, 1)
    for idx in range(58):
        month[idx] = d.month
        day[idx] = d.day
        label[idx] = d.strftime("%b-%d")
        d += _dt.timedelta(days=1)

    reserved[58] = 1
    month[59] = 2
    day[59] = 29
    label[59] = "Feb-29-composite"

    d = _dt.date(2001, 3, 1)
    for idx in range(60, 366):
        month[idx] = d.month
        day[idx] = d.day
        label[idx] = d.strftime("%b-%d")
        d += _dt.timedelta(days=1)
    return month, day, np.asarray(label, dtype=str), reserved

# =============================================================================
# 4. INPUT FILE INDEX / GRID / UNITS
# =============================================================================

def extract_year_month(path: Path, year: int) -> Optional[int]:
    m = re.search(r"(?<!\d)(\d{4})(\d{2})(?!\d)", path.name)
    if not m:
        return None
    y, mon = int(m.group(1)), int(m.group(2))
    return mon if y == year and 1 <= mon <= 12 else None


def build_file_index(year: int, folder: Path) -> Dict[int, Path]:
    index: Dict[int, Path] = {}
    for p in sorted(folder.glob(f"*{year}*.nc")):
        mon = extract_year_month(p, year)
        if mon is None:
            continue
        if mon in index:
            raise RuntimeError(f"Duplicate month {year}-{mon:02d}: {p}")
        index[mon] = p
    if set(index) != set(range(1, 13)):
        missing = sorted(set(range(1, 13)) - set(index))
        raise RuntimeError(f"Missing months for {year} in {folder}: {missing}")
    return index


def open_dataset(path: Path) -> xr.Dataset:
    return xr.open_dataset(
        path,
        engine="netcdf4",
        decode_times=True,
        mask_and_scale=True,
        cache=False,
    )


def sort_dataset(ds: xr.Dataset) -> xr.Dataset:
    return ds.sortby(["latitude", "longitude"])


def _norm_units(units: object) -> str:
    return str(units or "").strip().lower().replace("°", "deg")


def convert_temperature(values: np.ndarray, units: object, var_name: str) -> np.ndarray:
    u = _norm_units(units)
    x = values.astype(np.float32, copy=False)
    if u in {"k", "kelvin"}:
        return x - np.float32(273.15)
    if u in {"degc", "degree_celsius", "degrees_celsius", "c", "celsius"}:
        return x
    raise RuntimeError(f"{var_name}: unsupported temperature units {units!r}")


def convert_pressure(values: np.ndarray, units: object, var_name: str) -> np.ndarray:
    u = _norm_units(units)
    x = values.astype(np.float32, copy=False)
    if u in {"pa", "pascal", "pascals"}:
        return x / np.float32(100.0)
    if u in {"hpa", "mb", "millibar", "millibars"}:
        return x
    raise RuntimeError(f"{var_name}: unsupported pressure units {units!r}")


def validate_time_axis(ds: xr.Dataset, year: int, month: int) -> None:
    times = ds.time.values
    if times.size == 0:
        raise RuntimeError(f"Empty time axis: {year}-{month:02d}")
    diffs = np.diff(times).astype("timedelta64[m]").astype(np.int64)
    if np.any(diffs <= 0):
        raise RuntimeError(f"Non-increasing or duplicate timestamps: {year}-{month:02d}")
    if np.any(diffs != 60):
        bad = np.unique(diffs[diffs != 60])[:8]
        raise RuntimeError(f"Expected hourly ERA5-Land time steps in {year}-{month:02d}; found minute deltas {bad.tolist()}")


def validate_grids_and_axes(
    ds_t: xr.Dataset, ds_d: xr.Dataset, ds_p: xr.Dataset, year: int, month: int
) -> None:
    required = ("time", "latitude", "longitude")
    for name, ds, var in (("T2m", ds_t, "t2m"), ("D2m", ds_d, "d2m"), ("SP", ds_p, "sp")):
        if not all(dim in ds.dims for dim in required):
            raise RuntimeError(f"{name}: required dimensions missing in {year}-{month:02d}")
        if var not in ds.data_vars:
            raise RuntimeError(f"{name}: expected variable {var!r} not found in {year}-{month:02d}")
        validate_time_axis(ds, year, month)

    pairs = [
        ("time T/D", ds_t.time.values, ds_d.time.values),
        ("time T/SP", ds_t.time.values, ds_p.time.values),
        ("latitude T/D", ds_t.latitude.values, ds_d.latitude.values),
        ("latitude T/SP", ds_t.latitude.values, ds_p.latitude.values),
        ("longitude T/D", ds_t.longitude.values, ds_d.longitude.values),
        ("longitude T/SP", ds_t.longitude.values, ds_p.longitude.values),
    ]
    for name, a, b in pairs:
        if a.shape != b.shape or not np.array_equal(a, b):
            raise RuntimeError(f"{name} coordinates differ for {year}-{month:02d}")

    for name in ("T2m", "D2m", "SP"):
        pass

    expected_year = year
    expected_month = month
    # xarray time decoding supports standard datetime64 for ERA5-Land.
    for ts in ds_t.time.values[[0, -1]]:
        text = np.datetime_as_string(ts, unit="m")
        if int(text[:4]) != expected_year or int(text[5:7]) != expected_month:
            raise RuntimeError(f"Time axis leaves requested month: {year}-{month:02d}, sample={text}")

# =============================================================================
# 5. PHYSICS
# =============================================================================

def es_water(temp_c: np.ndarray) -> np.ndarray:
    return 6.112 * np.exp((17.67 * temp_c) / (temp_c + 243.5))


def es_ice(temp_c: np.ndarray) -> np.ndarray:
    return 6.112 * np.exp((22.46 * temp_c) / (temp_c + 272.62))


def saturation_vapor_pressure(temp_c: np.ndarray) -> np.ndarray:
    out = np.full(temp_c.shape, np.nan, dtype=np.float32)
    water = np.isfinite(temp_c) & (temp_c >= 0.0)
    ice = np.isfinite(temp_c) & ~water
    if np.any(water):
        out[water] = es_water(temp_c[water]).astype(np.float32)
    if np.any(ice):
        out[ice] = es_ice(temp_c[ice]).astype(np.float32)
    return out


def derive_moisture(T: np.ndarray, Td: np.ndarray, P_hpa: np.ndarray) -> dict[str, np.ndarray]:
    es_t = saturation_vapor_pressure(T)
    e = saturation_vapor_pressure(Td)

    rh_raw = 100.0 * (e / es_t)
    supersat = np.isfinite(rh_raw) & (rh_raw > 100.0)
    rh = np.clip(rh_raw, 0.0, 100.0).astype(np.float32)

    valid_qr = np.isfinite(e) & np.isfinite(P_hpa) & (e > 0.0) & (P_hpa > 0.0) & (e < P_hpa)
    r = np.full_like(e, np.nan, dtype=np.float32)
    r[valid_qr] = (0.622 * e[valid_qr] / (P_hpa[valid_qr] - e[valid_qr])).astype(np.float32)
    q = (r / (1.0 + r)).astype(np.float32)

    valid_all = np.isfinite(rh) & np.isfinite(e) & np.isfinite(r) & np.isfinite(q)

    return {
        "rh": rh,
        "e": e.astype(np.float32),
        "r": r,
        "q": q,
        "supersat": supersat,
        "valid_all": valid_all,
        "invalid_e_over_p": np.isfinite(e) & np.isfinite(P_hpa) & (e >= P_hpa),
    }

# =============================================================================
# 6. ONLINE STATISTICS
# =============================================================================

def update_moments_4_order(
    n: np.ndarray,
    mean: np.ndarray,
    M2: np.ndarray,
    M3: np.ndarray,
    M4: np.ndarray,
    x: np.ndarray,
    mask: np.ndarray,
    *,
    increment_n: bool = True,
) -> None:
    """Update one variable's online moments using a shared paired-observation count.

    When ``increment_n`` is True, the shared count ``n`` is advanced once.
    Subsequent variables can use ``increment_n=False`` so the same observation
    contributes to their moments without multiplying the shared count.
    """
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return

    n_current = n[idx].astype(np.float64)
    n_old = n_current if increment_n else np.maximum(n_current - 1.0, 0.0)
    n_new = n_old + 1.0

    x_obs = x[idx].astype(np.float64)
    mean_old = mean[idx]
    M2_old = M2[idx]
    M3_old = M3[idx]
    M4_old = M4[idx]

    delta = x_obs - mean_old
    mean_new = mean_old + delta / n_new
    M2_new = M2_old + delta * (x_obs - mean_new)
    M3_new = (
        M3_old
        + delta**3 * n_old * (n_old - 1.0) / n_new**2
        - 3.0 * delta * M2_old / n_new
    )
    M4_new = (
        M4_old
        + delta**4 * n_old * (n_old**2 - n_old + 1.0) / n_new**3
        + 6.0 * delta**2 * M2_old / n_new**2
        - 4.0 * delta * M3_old / n_new
    )

    if increment_n:
        n[idx] = n_new.astype(np.int64)
    mean[idx] = mean_new
    M2[idx] = M2_new
    M3[idx] = M3_new
    M4[idx] = M4_new


def update_covariance(
    n_pair: np.ndarray,
    mean_x: np.ndarray,
    mean_y: np.ndarray,
    Cxy: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    mask: np.ndarray,
    *,
    increment_n: bool = True,
) -> None:
    """Update covariance using a shared observation count.

    With ``increment_n=False`` the covariance is aligned to the already-updated
    shared count and therefore does not inflate the observation count.
    """
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return

    n_current = n_pair[idx].astype(np.float64)
    n_old = n_current if increment_n else np.maximum(n_current - 1.0, 0.0)
    n_new = n_old + 1.0

    xo = x[idx].astype(np.float64)
    yo = y[idx].astype(np.float64)
    mx = mean_x[idx].astype(np.float64)
    my = mean_y[idx].astype(np.float64)

    dx = xo - mx
    mx_new = mx + dx / n_new
    my_new = my + (yo - my) / n_new
    C_new = Cxy[idx] + dx * (yo - my_new)

    if increment_n:
        n_pair[idx] = n_new.astype(np.int64)
    mean_x[idx] = mx_new
    mean_y[idx] = my_new
    Cxy[idx] = C_new


def combine_moments(
    n1: np.ndarray, m1: np.ndarray, M21: np.ndarray, M31: np.ndarray, M41: np.ndarray,
    n2: np.ndarray, m2: np.ndarray, M22: np.ndarray, M32: np.ndarray, M42: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n1f = n1.astype(np.float64)
    n2f = n2.astype(np.float64)
    nt = n1f + n2f
    out_n = nt.astype(np.int64)
    out_m = m1.copy()
    out_M2 = M21.copy()
    out_M3 = M31.copy()
    out_M4 = M41.copy()

    only2 = (n1 == 0) & (n2 > 0)
    both = (n1 > 0) & (n2 > 0)
    if np.any(only2):
        out_m[only2] = m2[only2]
        out_M2[only2] = M22[only2]
        out_M3[only2] = M32[only2]
        out_M4[only2] = M42[only2]
    if np.any(both):
        a = n1f[both]; b = n2f[both]; n_tot = nt[both]
        delta = m2[both] - m1[both]
        out_m[both] = m1[both] + delta * b / n_tot
        out_M2[both] = M21[both] + M22[both] + delta**2 * a * b / n_tot
        out_M3[both] = (
            M31[both] + M32[both]
            + delta**3 * a * b * (a - b) / n_tot**2
            + 3.0 * delta * (a * M22[both] - b * M21[both]) / n_tot
        )
        out_M4[both] = (
            M41[both] + M42[both]
            + delta**4 * a * b * (a*a - a*b + b*b) / n_tot**3
            + 6.0 * delta**2 * (a*a * M22[both] + b*b * M21[both]) / n_tot**2
            + 4.0 * delta * (a * M32[both] - b * M31[both]) / n_tot
        )
    return out_n, out_m, out_M2, out_M3, out_M4


def combine_covariance(
    n1: np.ndarray, mx1: np.ndarray, my1: np.ndarray, c1: np.ndarray,
    n2: np.ndarray, mx2: np.ndarray, my2: np.ndarray, c2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    a = n1.astype(np.float64)
    b = n2.astype(np.float64)
    nt = a + b
    out_mx = mx1.copy(); out_my = my1.copy(); out_c = c1.copy(); out_n = nt.astype(np.int64)
    only2 = (n1 == 0) & (n2 > 0)
    both = (n1 > 0) & (n2 > 0)
    if np.any(only2):
        out_mx[only2] = mx2[only2]
        out_my[only2] = my2[only2]
        out_c[only2] = c2[only2]
    if np.any(both):
        aa = a[both]; bb = b[both]; nn = nt[both]
        dx = mx2[both] - mx1[both]
        dy = my2[both] - my1[both]
        out_mx[both] = mx1[both] + dx * bb / nn
        out_my[both] = my1[both] + dy * bb / nn
        out_c[both] = c1[both] + c2[both] + dx * dy * aa * bb / nn
    return out_n, out_mx, out_my, out_c

# =============================================================================
# 7. BIVARIATE PROBABILITY FUNCTIONS
# =============================================================================

# Primary empirical product: RH in [0,100] %, q in [0,1].
# The output is a piecewise-constant PDF reconstructed from exact hourly counts.
BIVARIATE_RANGES = {
    "rh": (0.0, 100.0),
    "q": (0.0, 1.0),
}

def bivariate_histogram_edges(xname: str, yname: str) -> tuple[np.ndarray, np.ndarray]:
    if xname not in BIVARIATE_RANGES or yname not in BIVARIATE_RANGES:
        raise ValueError(f"No fixed physical histogram range for {xname},{yname}; empirical histogram pass supports only configured bounded variables.")
    xrng = BIVARIATE_RANGES[xname]
    yrng = BIVARIATE_RANGES[yname]
    return (
        np.linspace(xrng[0], xrng[1], BIVARIATE_NX + 1, dtype=np.float64),
        np.linspace(yrng[0], yrng[1], BIVARIATE_NY + 1, dtype=np.float64),
    )

def empirical_bivariate_histogram(values_x: np.ndarray, values_y: np.ndarray,
                                  x_edges: np.ndarray, y_edges: np.ndarray) -> tuple[np.ndarray, int]:
    """Return exact 2-D bin counts and paired-valid sample count."""
    x=np.asarray(values_x, dtype=np.float64).reshape(-1)
    y=np.asarray(values_y, dtype=np.float64).reshape(-1)
    mask=np.isfinite(x)&np.isfinite(y)
    mask &= (x >= x_edges[0]) & (x <= x_edges[-1]) & (y >= y_edges[0]) & (y <= y_edges[-1])
    if not np.any(mask):
        return np.zeros((len(x_edges)-1,len(y_edges)-1), dtype=np.uint16), 0
    xx=x[mask]; yy=y[mask]
    ix=np.searchsorted(x_edges, xx, side="right")-1
    iy=np.searchsorted(y_edges, yy, side="right")-1
    ix=np.minimum(ix, len(x_edges)-2); iy=np.minimum(iy, len(y_edges)-2)
    flat=ix*(len(y_edges)-1)+iy
    counts=np.bincount(flat, minlength=(len(x_edges)-1)*(len(y_edges)-1))
    if int(counts.max(initial=0)) > np.iinfo(np.uint16).max:
        raise OverflowError("Bivariate histogram count exceeds uint16 capacity; use a larger integer dtype.")
    return counts.reshape(len(x_edges)-1,len(y_edges)-1).astype(np.uint16), int(mask.sum())

def update_empirical_histogram_inplace(hist: np.ndarray, x: np.ndarray, y: np.ndarray) -> None:
    """v8 FINAL single-pass exact RH-q empirical histogram update."""
    x=np.asarray(x,dtype=np.float64).reshape(-1)
    y=np.asarray(y,dtype=np.float64).reshape(-1)
    valid=np.isfinite(x)&np.isfinite(y)&(x>=0)&(x<=100)&(y>=0)&(y<=1)
    if not np.any(valid):
        return
    ix=np.clip((x[valid] / 100.0 * BIVARIATE_NX).astype(np.int64),0,BIVARIATE_NX-1)
    iy=np.clip((y[valid] * BIVARIATE_NY).astype(np.int64),0,BIVARIATE_NY-1)
    hist += np.bincount(ix*BIVARIATE_NY+iy, minlength=BIVARIATE_NX*BIVARIATE_NY).reshape(BIVARIATE_NX,BIVARIATE_NY)


def empirical_bivariate_pdf(x: np.ndarray | float, y: np.ndarray | float,
                            counts: np.ndarray, n_valid: np.ndarray | int,
                            x_edges: np.ndarray, y_edges: np.ndarray) -> np.ndarray | float:
    """Evaluate the piecewise-constant empirical 2-D PDF represented by bin counts."""
    x,y,n=np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float),np.asarray(n_valid,float))
    out=np.full(x.shape,np.nan,float)
    den_area=np.diff(x_edges)[:,None]*np.diff(y_edges)[None,:]
    ix=np.searchsorted(x_edges,x,side="right")-1; iy=np.searchsorted(y_edges,y,side="right")-1
    ix=np.clip(ix,0,len(x_edges)-2); iy=np.clip(iy,0,len(y_edges)-2)
    valid=np.isfinite(x)&np.isfinite(y)&np.isfinite(n)&(n>0)&(x>=x_edges[0])&(x<=x_edges[-1])&(y>=y_edges[0])&(y<=y_edges[-1])
    if np.any(valid):
        c=np.asarray(counts)
        out[valid]=c[ix[valid],iy[valid]]/(n[valid]*den_area[ix[valid],iy[valid]])
    return out.item() if out.ndim==0 else out

def bivariate_gaussian_pdf(
    x: np.ndarray | float, y: np.ndarray | float, mean_x: np.ndarray | float,
    std_x: np.ndarray | float, mean_y: np.ndarray | float, std_y: np.ndarray | float,
    rho: np.ndarray | float,
) -> np.ndarray | float:
    """Reference parametric candidate only; v8.0 FINAL SINGLE-PASS does not assume Gaussianity."""
    x, y, mx, sx, my, sy, r = np.broadcast_arrays(
        np.asarray(x,dtype=float), np.asarray(y,dtype=float), np.asarray(mean_x,dtype=float),
        np.asarray(std_x,dtype=float), np.asarray(mean_y,dtype=float), np.asarray(std_y,dtype=float), np.asarray(rho,dtype=float))
    valid=np.isfinite(x)&np.isfinite(y)&np.isfinite(mx)&np.isfinite(my)&np.isfinite(sx)&np.isfinite(sy)&np.isfinite(r)&(sx>0)&(sy>0)&(np.abs(r)<1)
    out=np.full(x.shape,np.nan,float)
    if np.any(valid):
        rr=r[valid]; zx=(x[valid]-mx[valid])/sx[valid]; zy=(y[valid]-my[valid])/sy[valid]
        den=2*np.pi*sx[valid]*sy[valid]*np.sqrt(1-rr**2)
        out[valid]=np.exp(-(zx**2-2*rr*zx*zy+zy**2)/(2*(1-rr**2)))/den
    return out.item() if out.ndim==0 else out

# =============================================================================
# 8. CHECKPOINT CONTRACT
# =============================================================================

STATE_VARS = {
    "n": "i8",
    "mean": "f8",
    "M2": "f8",
    "M3": "f8",
    "M4": "f8",
}


def year_paths(year: int) -> tuple[Path, Path, Path]:
    base = YEAR_DIR / f"year_{year:04d}_{CONFIG_HASH}"
    return base.with_suffix(".nc"), base.with_suffix(".json"), base.with_suffix(".part.nc")


def is_year_complete(year: int) -> bool:
    nc_path, js_path, _ = year_paths(year)
    if not nc_path.exists() or not js_path.exists():
        return False
    try:
        meta = json.loads(js_path.read_text(encoding="utf-8"))
        return (
            meta.get("status") == "completed"
            and meta.get("year") == year
            and meta.get("schema_version") == CHECKPOINT_VERSION
            and meta.get("config_hash") == CONFIG_HASH
            and meta.get("sha256") == sha256_file(nc_path)
        )
    except Exception:
        return False


def create_year_checkpoint(path: Path, lat: np.ndarray, lon: np.ndarray) -> Dataset:
    require_netcdf4()
    if path.exists():
        return Dataset(path, "r+")
    path.parent.mkdir(parents=True, exist_ok=True)
    ds = Dataset(path, "w", format="NETCDF4")
    ds.createDimension("doy", DOY_COUNT)
    ds.createDimension("latitude", len(lat))
    ds.createDimension("longitude", len(lon))
    ds.createDimension("rh_bin", BIVARIATE_NX)
    ds.createDimension("q_bin", BIVARIATE_NY)
    ds.createDimension("y_chunk", (len(lat) + CHUNK_LAT - 1) // CHUNK_LAT)
    ds.createDimension("x_chunk", (len(lon) + CHUNK_LON - 1) // CHUNK_LON)
    ds.createVariable("doy", "i2", ("doy",))[:] = np.arange(1, DOY_COUNT + 1, dtype=np.int16)
    ds.createVariable("latitude", "f4", ("latitude",))[:] = lat.astype(np.float32)
    ds.createVariable("longitude", "f4", ("longitude",))[:] = lon.astype(np.float32)
    month, day, label, reserved = calendar_labels()
    ds.createVariable("month", "i2", ("doy",))[:] = month
    ds.createVariable("day", "i2", ("doy",))[:] = day
    ds.createVariable("reserved_day", "i1", ("doy",))[:] = reserved

    chunks = (1, min(CHUNK_LAT, len(lat)), min(CHUNK_LON, len(lon)))
    ds.createVariable("completed_chunk", "i1", ("doy", "y_chunk", "x_chunk"),
                      zlib=True, complevel=4, shuffle=True, fill_value=0)
    ds.createVariable("n_obs", "i8", ("doy", "latitude", "longitude"),
                      zlib=True, complevel=4, shuffle=True, chunksizes=chunks, fill_value=0)
    for key in ("rh", "e", "r", "q"):
        for field in ("mean", "M2", "M3", "M4"):
            ds.createVariable(f"{field}_{key}", "f8", ("doy", "latitude", "longitude"),
                              zlib=True, complevel=4, shuffle=True, chunksizes=chunks, fill_value=0.0)
    for name in ("total_supersat", "total_invalid_ep"):
        ds.createVariable(name, "i8", ("doy", "latitude", "longitude"),
                          zlib=True, complevel=4, shuffle=True, chunksizes=chunks, fill_value=0)

    for xname, yname in BIVARIATE_PAIRS:
        tag = f"{xname}__{yname}"
        for field in ("mean_x", "mean_y", "Cxy"):
            ds.createVariable(f"pair_{tag}_{field}", "f8", ("doy", "latitude", "longitude"),
                              zlib=True, complevel=4, shuffle=True, chunksizes=chunks, fill_value=0.0)
        # v8 FINAL: exact empirical histogram accumulated in the same hourly pass
        ds.createVariable(f"pair_{tag}_hist", "i8",
                          ("doy", "latitude", "longitude", "rh_bin", "q_bin"),
                          zlib=True, complevel=4, shuffle=True, fill_value=0)

    ds.schema_version = CHECKPOINT_VERSION
    ds.config_hash = CONFIG_HASH
    ds.script_sha256 = script_sha256()
    ds.purpose = "Annual disk-backed hourly empirical sufficient-statistic checkpoint; spatial blocks and completion bitmap provide power-failure restart"
    ds.chunk_lat = CHUNK_LAT
    ds.chunk_lon = CHUNK_LON
    ds.progress_flush_chunks = PROGRESS_FLUSH_CHUNKS
    ds.sync()
    return ds

# =============================================================================
# 9. YEAR PROCESSING - DAILY MEMORY, ANNUAL CHECKPOINT
# =============================================================================

def _empty_daily_states(ncells: int) -> dict:
    out = {"n": np.zeros(ncells, np.int64), "total_supersat": np.zeros(ncells, np.int64),
           "total_invalid_ep": np.zeros(ncells, np.int64)}
    for key in ("rh", "e", "r", "q"):
        for field in ("mean", "M2", "M3", "M4"):
            out[f"{field}_{key}"] = np.zeros(ncells, np.float64)
    for xname, yname in BIVARIATE_PAIRS:
        tag = f"{xname}__{yname}"
        out[f"pair_{tag}_mean_x"] = np.zeros(ncells, np.float64)
        out[f"pair_{tag}_mean_y"] = np.zeros(ncells, np.float64)
        out[f"pair_{tag}_Cxy"] = np.zeros(ncells, np.float64)
        out[f"pair_{tag}_hist"] = np.zeros((ncells, BIVARIATE_NX, BIVARIATE_NY), np.int64)
    return out


def _write_daily_state(ds: Dataset, doy0: int, s: dict, j0: int, j1: int, i0: int, i1: int) -> None:
    ds.variables["n_obs"][doy0, j0:j1, i0:i1] = s["n"].reshape(j1 - j0, i1 - i0)
    for key in ("rh", "e", "r", "q"):
        for field in ("mean", "M2", "M3", "M4"):
            ds.variables[f"{field}_{key}"][doy0, j0:j1, i0:i1] = s[f"{field}_{key}"].reshape(j1 - j0, i1 - i0)
    ds.variables["total_supersat"][doy0, j0:j1, i0:i1] = s["total_supersat"].reshape(j1 - j0, i1 - i0)
    ds.variables["total_invalid_ep"][doy0, j0:j1, i0:i1] = s["total_invalid_ep"].reshape(j1 - j0, i1 - i0)
    for xname, yname in BIVARIATE_PAIRS:
        tag = f"{xname}__{yname}"
        ds.variables[f"pair_{tag}_mean_x"][doy0, j0:j1, i0:i1] = s[f"pair_{tag}_mean_x"].reshape(j1 - j0, i1 - i0)
        ds.variables[f"pair_{tag}_mean_y"][doy0, j0:j1, i0:i1] = s[f"pair_{tag}_mean_y"].reshape(j1 - j0, i1 - i0)
        ds.variables[f"pair_{tag}_Cxy"][doy0, j0:j1, i0:i1] = s[f"pair_{tag}_Cxy"].reshape(j1 - j0, i1 - i0)
        ds.variables[f"pair_{tag}_hist"][doy0, j0:j1, i0:i1, :, :] = s[f"pair_{tag}_hist"].reshape(j1-j0, i1-i0, BIVARIATE_NX, BIVARIATE_NY)



def _safe_completed_flag(var, doy_index: int, y_chunk_index: int, x_chunk_index: int) -> int:
    """Read a completion flag safely even when netCDF4 exposes an unset fill value as masked."""
    raw = var[doy_index, y_chunk_index, x_chunk_index]
    return int(np.ma.filled(raw, 0))


def _read_progress(js_path: Path) -> dict:
    if not js_path.exists():
        return {}
    try:
        return json.loads(js_path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _write_progress(
    json_path: Path,
    *,
    year: int,
    ny: int,
    nx: int,
    completed_dates: set[str],
    completed_units: int,
    total_units: int,
    last_doy: int | None = None,
) -> None:
    pct = 100.0 * completed_units / max(total_units, 1)
    atomic_json_write(
        json_path,
        {
            "status": "running",
            "year": year,
            "schema_version": CHECKPOINT_VERSION,
            "config_hash": CONFIG_HASH,
            "shape": [ny, nx],
            "chunk_shape": [CHUNK_LAT, CHUNK_LON],
            "completed_native_dates": sorted(completed_dates),
            "completed_units": int(completed_units),
            "total_units": int(total_units),
            "remaining_units": int(max(total_units - completed_units, 0)),
            "progress_percent": float(pct),
            "last_completed_doy": last_doy,
            "updated_utc": datetime.now(timezone.utc).isoformat(),
        },
    )


def _log_progress(year: int, completed_units: int, total_units: int, elapsed: float, context: str) -> None:
    pct = 100.0 * completed_units / max(total_units, 1)
    remaining = max(total_units - completed_units, 0)
    rate = completed_units / max(elapsed, 1e-9)
    eta = remaining / rate if rate > 0 else float("inf")
    eta_txt = f"{eta/3600:.2f} h" if np.isfinite(eta) else "unknown"
    logger.info(
        "PROGRESS | Year %04d | %s | %.2f%% | %d/%d units | remaining %d | rate %.2f units/s | ETA %s",
        year, context, pct, completed_units, total_units, remaining, rate, eta_txt,
    )

def process_year_empirical(year: int) -> tuple[int, Optional[Path]]:
    """Process one year with chunk-first monthly reads and deferred NetCDF sync.

    Scientific results and checkpoint schema are unchanged. The optimization changes
    only I/O scheduling and histogram accumulation:
      * each monthly spatial chunk is read once instead of once per DOY/hour;
      * np.add.at is replaced by equivalent np.bincount accumulation;
      * NetCDF sync is performed every PROGRESS_FLUSH_CHUNKS completed chunks
        and at DOY/month/year boundaries, rather than after every chunk.
    """
    logger.info("Starting hourly empirical processing for year %d", year)
    t0 = time.time()
    final_path, json_path, part_path = year_paths(year)
    if is_year_complete(year):
        return year, final_path

    ds_ckpt = None
    try:
        t_idx = build_file_index(year, T2M_DIR)
        d_idx = build_file_index(year, D2M_DIR)
        p_idx = build_file_index(year, SP_DIR)

        with open_dataset(t_idx[1]) as ds0:
            ds0 = sort_dataset(ds0)
            ny, nx = ds0.sizes["latitude"], ds0.sizes["longitude"]
            lat, lon = ds0.latitude.values, ds0.longitude.values
        require_netcdf4()

        n_y_chunks = (ny + CHUNK_LAT - 1) // CHUNK_LAT
        n_x_chunks = (nx + CHUNK_LON - 1) // CHUNK_LON
        total_slots_per_year = 365
        total_units = total_slots_per_year * n_y_chunks * n_x_chunks

        ds_ckpt = Dataset(part_path, "r+") if part_path.exists() else create_year_checkpoint(part_path, lat, lon)
        completed_dates = set(_read_progress(json_path).get("completed_native_dates", []))
        completed_units = int(np.asarray(ds_ckpt.variables["completed_chunk"][:], dtype=np.int64).sum())
        _write_progress(
            json_path, year=year, ny=ny, nx=nx, completed_dates=completed_dates,
            completed_units=completed_units, total_units=total_units,
        )
        _log_progress(year, completed_units, total_units, time.time() - t0, "resume/start")

        for month in range(1, 13):
            logger.info("Year %d: processing month %02d", year, month)
            with open_dataset(t_idx[month]) as ds_t, open_dataset(d_idx[month]) as ds_d, open_dataset(p_idx[month]) as ds_p:
                ds_t = sort_dataset(ds_t)
                ds_d = sort_dataset(ds_d)
                ds_p = sort_dataset(ds_p)
                validate_grids_and_axes(ds_t, ds_d, ds_p, year, month)

                tu = ds_t["t2m"].attrs.get("units")
                du = ds_d["d2m"].attrs.get("units")
                pu = ds_p["sp"].attrs.get("units")
                native_doys = ds_t.time.dt.dayofyear.values.astype(np.int16)
                date_texts = np.asarray([np.datetime_as_string(t, unit="D") for t in ds_t.time.values])

                # Build the climatological-slot index once per month.
                slot_time_indices: dict[int, np.ndarray] = {}
                slot_dates_map: dict[int, list[str]] = {}
                for native_doy in np.unique(native_doys):
                    cdoy = get_clim_doy(int(native_doy), year)
                    if cdoy < 1 or cdoy > DOY_COUNT or cdoy == 59:
                        continue
                    idx = np.flatnonzero(np.asarray(
                        [get_clim_doy(int(d), year) == cdoy for d in native_doys],
                        dtype=bool
                    ))
                    if idx.size:
                        slot_time_indices[cdoy] = idx
                        slot_dates_map[cdoy] = sorted(set(date_texts[idx].tolist()))

                if not slot_time_indices:
                    continue

                slot_completed_counts = {
                    cdoy: sum(
                        _safe_completed_flag(
                            ds_ckpt.variables["completed_chunk"], cdoy - 1, yyci, xxci
                        )
                        for yyci in range(n_y_chunks)
                        for xxci in range(n_x_chunks)
                    )
                    for cdoy in slot_time_indices
                }

                for yci, j0 in enumerate(range(0, ny, CHUNK_LAT)):
                    j1 = min(j0 + CHUNK_LAT, ny)
                    for xci, i0 in enumerate(range(0, nx, CHUNK_LON)):
                        i1 = min(i0 + CHUNK_LON, nx)

                        pending_slots = [
                            cdoy for cdoy in sorted(slot_time_indices)
                            if _safe_completed_flag(
                                ds_ckpt.variables["completed_chunk"], cdoy - 1, yci, xci
                            ) == 0
                        ]
                        if not pending_slots:
                            continue

                        # Critical I/O optimization: load this monthly spatial chunk once.
                        T_month = ds_t["t2m"].isel(
                            latitude=slice(j0, j1), longitude=slice(i0, i1)
                        ).values
                        Td_month = ds_d["d2m"].isel(
                            latitude=slice(j0, j1), longitude=slice(i0, i1)
                        ).values
                        P_month = ds_p["sp"].isel(
                            latitude=slice(j0, j1), longitude=slice(i0, i1)
                        ).values

                        for cdoy in pending_slots:
                            time_idx = slot_time_indices[cdoy]
                            slot_dates = slot_dates_map[cdoy]
                            s = _empty_daily_states((j1 - j0) * (i1 - i0))

                            for ti in time_idx:
                                T = convert_temperature(
                                    T_month[int(ti)], tu, "t2m"
                                ).reshape(-1)
                                Td = convert_temperature(
                                    Td_month[int(ti)], du, "d2m"
                                ).reshape(-1)
                                P = convert_pressure(
                                    P_month[int(ti)], pu, "sp"
                                ).reshape(-1)
                                phys = derive_moisture(T, Td, P)
                                mask = phys["valid_all"]

                                update_moments_4_order(
                                    s["n"], s["mean_rh"], s["M2_rh"], s["M3_rh"], s["M4_rh"],
                                    phys["rh"], mask, increment_n=True
                                )
                                update_moments_4_order(
                                    s["n"], s["mean_e"], s["M2_e"], s["M3_e"], s["M4_e"],
                                    phys["e"], mask, increment_n=False
                                )
                                update_moments_4_order(
                                    s["n"], s["mean_r"], s["M2_r"], s["M3_r"], s["M4_r"],
                                    phys["r"], mask, increment_n=False
                                )
                                update_moments_4_order(
                                    s["n"], s["mean_q"], s["M2_q"], s["M3_q"], s["M4_q"],
                                    phys["q"], mask, increment_n=False
                                )

                                for xname, yname in BIVARIATE_PAIRS:
                                    tag = f"{xname}__{yname}"
                                    update_covariance(
                                        s["n"],
                                        s[f"pair_{tag}_mean_x"],
                                        s[f"pair_{tag}_mean_y"],
                                        s[f"pair_{tag}_Cxy"],
                                        phys[xname],
                                        phys[yname],
                                        mask,
                                        increment_n=False,
                                    )

                                    flat_hist = s[f"pair_{tag}_hist"]
                                    h = flat_hist.reshape(-1, BIVARIATE_NX, BIVARIATE_NY)
                                    valid_idx = np.flatnonzero(mask)
                                    if valid_idx.size:
                                        xx = phys[xname][valid_idx]
                                        yy = phys[yname][valid_idx]
                                        bins_x = np.clip(
                                            (xx / 100.0 * BIVARIATE_NX).astype(np.int64),
                                            0, BIVARIATE_NX - 1,
                                        )
                                        bins_y = np.clip(
                                            (yy * BIVARIATE_NY).astype(np.int64),
                                            0, BIVARIATE_NY - 1,
                                        )
                                        flat_idx = (
                                            valid_idx * (BIVARIATE_NX * BIVARIATE_NY)
                                            + bins_x * BIVARIATE_NY
                                            + bins_y
                                        )
                                        h += np.bincount(
                                            flat_idx,
                                            minlength=h.size,
                                        ).reshape(h.shape)

                                s["total_supersat"] += phys["supersat"].astype(np.int64)
                                s["total_invalid_ep"] += phys["invalid_e_over_p"].astype(np.int64)

                            _write_daily_state(
                                ds_ckpt, cdoy - 1, s, j0, j1, i0, i1
                            )
                            ds_ckpt.variables["completed_chunk"][cdoy - 1, yci, xci] = 1
                            completed_units += 1
                            slot_completed_counts[cdoy] += 1

                            if completed_units % PROGRESS_FLUSH_CHUNKS == 0:
                                ds_ckpt.sync()
                                _write_progress(
                                    json_path, year=year, ny=ny, nx=nx,
                                    completed_dates=completed_dates,
                                    completed_units=completed_units,
                                    total_units=total_units, last_doy=cdoy,
                                )

                            if (
                                completed_units % PROGRESS_LOG_EVERY_CHUNKS == 0
                                or completed_units == total_units
                            ):
                                _log_progress(
                                    year, completed_units, total_units,
                                    time.time() - t0,
                                    f"DOY {cdoy:03d} chunk {yci + 1}/{n_y_chunks},{xci + 1}/{n_x_chunks}",
                                )

                            # Mark a DOY complete only after every spatial chunk is committed.
                            if slot_completed_counts[cdoy] == n_y_chunks * n_x_chunks:
                                completed_dates.update(slot_dates)
                                ds_ckpt.sync()
                                _write_progress(
                                    json_path, year=year, ny=ny, nx=nx,
                                    completed_dates=completed_dates,
                                    completed_units=completed_units,
                                    total_units=total_units, last_doy=cdoy,
                                )
                                logger.info(
                                    "CHECKPOINT COMMIT | Year %04d | DOY %03d | %.2f%% complete | remaining %d units",
                                    year, cdoy, 100.0 * completed_units / total_units,
                                    total_units - completed_units,
                                )

                        # Monthly chunk arrays go out of scope before the next chunk.
                        del T_month, Td_month, P_month

                # Make month boundary durable.
                ds_ckpt.sync()

        ds_ckpt.sync()
        import calendar as _calendar
        expected_dates = {
            f"{year}-{month:02d}-{day:02d}"
            for month in range(1, 13)
            for day in range(1, _calendar.monthrange(year, month)[1] + 1)
        }
        if completed_dates != expected_dates:
            missing = sorted(expected_dates - completed_dates)
            raise RuntimeError(
                f"Year {year} incomplete; missing {len(missing)} native dates (first={missing[:3]})"
            )
        if completed_units != total_units:
            raise RuntimeError(
                f"Year {year} incomplete checkpoint: {completed_units}/{total_units} spatial-day units committed"
            )

        ds_ckpt.close()
        ds_ckpt = None
        os.replace(part_path, final_path)
        meta = {
            "status": "completed", "year": year, "schema_version": CHECKPOINT_VERSION,
            "config_hash": CONFIG_HASH, "sha256": sha256_file(final_path), "shape": [ny, nx],
            "chunk_shape": [CHUNK_LAT, CHUNK_LON], "completed_native_dates": sorted(completed_dates),
            "completed_units": int(completed_units), "total_units": int(total_units),
            "progress_percent": 100.0, "remaining_units": 0,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        atomic_json_write(json_path, meta)
        logger.info(
            "POWER-FAILURE SAFE COMMIT | Year %04d | 100.00%% | %d/%d units | remaining 0 | elapsed %.1fs",
            year, completed_units, total_units, time.time() - t0
        )
        return year, final_path
    except Exception:
        logger.exception("Year %d failed; disk checkpoint is retained for restart", year)
        try:
            if ds_ckpt is not None:
                ds_ckpt.sync()
                ds_ckpt.close()
        except Exception:
            pass
        return year, None
# =============================================================================
# 10. NETCDF OUTPUT HELPERS
# =============================================================================

def create_output_file(path: Path, lat: np.ndarray, lon: np.ndarray, title: str) -> Dataset:
    require_netcdf4()
    if path.exists():
        path.unlink()
    ds = Dataset(path, "w", format="NETCDF4")
    ds.createDimension("doy", DOY_COUNT)
    ds.createDimension("latitude", len(lat))
    ds.createDimension("longitude", len(lon))
    ds.createDimension("rh_bin", BIVARIATE_NX)
    ds.createDimension("q_bin", BIVARIATE_NY)
    ds.createVariable("doy", "i2", ("doy",))[:] = np.arange(1, DOY_COUNT + 1, dtype=np.int16)
    ds.createVariable("latitude", "f4", ("latitude",))[:] = lat.astype(np.float32)
    ds.createVariable("longitude", "f4", ("longitude",))[:] = lon.astype(np.float32)
    month, day, label, reserved = calendar_labels()
    ds.createVariable("month", "i2", ("doy",))[:] = month
    ds.createVariable("day", "i2", ("doy",))[:] = day
    ds.createVariable("reserved_day", "i1", ("doy",))[:] = reserved
    ds.title = title
    ds.period = f"{START_YEAR}-{END_YEAR}"
    ds.calendar = "366-slot; Feb-28 and Feb-29 pooled into slot 60; slot 59 reserved"
    ds.schema_version = SCHEMA_VERSION
    ds.config_hash = CONFIG_HASH
    ds.script_sha256 = script_sha256()
    return ds


def stats_from_state(n: np.ndarray, mean: np.ndarray, M2: np.ndarray, M3: np.ndarray, M4: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n64 = n.astype(np.float64)
    std = np.full(n.shape, np.nan, dtype=np.float64)
    skew = np.full(n.shape, np.nan, dtype=np.float64)
    kurt = np.full(n.shape, np.nan, dtype=np.float64)
    ok2 = n >= 2
    std[ok2] = np.sqrt(np.maximum(M2[ok2] / (n64[ok2] - 1.0), 0.0))
    ok3 = (n >= 3) & (M2 > 0)
    if np.any(ok3):
        m2 = M2[ok3] / n64[ok3]
        m3 = M3[ok3] / n64[ok3]
        nn = n64[ok3]
        skew[ok3] = np.sqrt(nn * (nn - 1.0)) / (nn - 2.0) * m3 / np.power(m2, 1.5)
    ok4 = (n >= 4) & (M2 > 0)
    if np.any(ok4):
        nn = n64[ok4]
        b2 = nn * M4[ok4] / np.square(M2[ok4])
        kurt[ok4] = ((nn - 1.0) / ((nn - 2.0) * (nn - 3.0))) * ((nn + 1.0) * b2 - 3.0 * (nn - 1.0))
    return mean.astype(np.float32), std.astype(np.float32), skew.astype(np.float32), kurt.astype(np.float32)

# =============================================================================
# 11. FINAL MERGE - DOY x CELL CHUNK
# =============================================================================

def finalize_all(years: Iterable[int], lat: np.ndarray, lon: np.ndarray) -> None:
    main = create_output_file(OUTPUT_FILE, lat, lon, "ERA5-Land Moisture Climatology 1981-2020 (v8.0 FINAL SINGLE-PASS hourly empirical)")
    diag = create_output_file(DIAGNOSTIC_FILE, lat, lon, "ERA5-Land Moisture Diagnostics 1981-2020 (v8.0 FINAL SINGLE-PASS)")
    biv = create_output_file(BIVARIATE_FILE, lat, lon, "ERA5-Land Bivariate Probability Parameters 1981-2020 (v8.0 FINAL SINGLE-PASS)")

    out_chunks = (1, min(CHUNK_LAT, len(lat)), min(CHUNK_LON, len(lon)))
    for key in ("rh", "e", "r", "q"):
        for field in ("mean", "std", "skew", "kurt"):
            main.createVariable(
                f"{field}_{key}", "f4", ("doy", "latitude", "longitude"), zlib=True, complevel=4,
                shuffle=True, chunksizes=out_chunks, fill_value=-9999.0
            )
            v = main.variables[f"{field}_{key}"]
            units = {"rh": "%", "e": "hPa", "r": "kg kg-1", "q": "kg kg-1"}[key]
            v.units = units
            v.long_name = f"{field} of {key} from hourly empirical accumulation"

    diag.createVariable("valid_observation_count", "i8", ("doy", "latitude", "longitude"),
                       zlib=True, complevel=4, shuffle=True, chunksizes=out_chunks, fill_value=-9999)
    diag.variables["valid_observation_count"].long_name = "Number of paired-valid hourly observations"
    diag.variables["valid_observation_count"].units = "1"
    diag.createVariable("supersaturation_fraction", "f4", ("doy", "latitude", "longitude"),
                       zlib=True, complevel=4, shuffle=True, chunksizes=out_chunks, fill_value=-9999.0)
    diag.variables["supersaturation_fraction"].long_name = "Fraction of valid-pressure hourly states with RH above 100 percent"
    diag.variables["supersaturation_fraction"].units = "1"
    diag.createVariable("invalid_e_over_p_fraction", "f4", ("doy", "latitude", "longitude"),
                       zlib=True, complevel=4, shuffle=True, chunksizes=out_chunks, fill_value=-9999.0)
    diag.variables["invalid_e_over_p_fraction"].long_name = "Fraction of hourly states with vapor pressure not below surface pressure"
    diag.variables["invalid_e_over_p_fraction"].units = "1"

    for xname, yname in BIVARIATE_PAIRS:
        tag = f"{xname}__{yname}"
        for field in ("n", "mean_x", "mean_y", "std_x", "std_y", "cov", "corr"):
            dtype = "i8" if field == "n" else "f4"
            biv.createVariable(
                f"{tag}_{field}", dtype, ("doy", "latitude", "longitude"),
                zlib=True, complevel=4, shuffle=True, chunksizes=out_chunks,
                fill_value=-9999 if field == "n" else -9999.0
            )
            vv = biv.variables[f"{tag}_{field}"]
            if field == "n":
                vv.long_name = f"Paired-valid hourly count for {xname} and {yname}"
                vv.units = "1"
            elif field in ("mean_x", "std_x"):
                vv.long_name = f"{field} for {xname}"
            elif field in ("mean_y", "std_y"):
                vv.long_name = f"{field} for {yname}"
            elif field == "cov":
                vv.long_name = f"Sample covariance of {xname} and {yname}"
            elif field == "corr":
                vv.long_name = f"Pearson correlation of {xname} and {yname}"
                vv.units = "1"

    main.method = "Direct hourly empirical accumulation; no Monte Carlo in primary products"
    main.bivariate_probability = "Empirical 2-D PDF is the primary joint product; Gaussian parameters are reference candidates only"
    main.bivariate_pairs = json.dumps([list(p) for p in BIVARIATE_PAIRS])
    diag.description = "Counts and physical-validity diagnostics for the same hourly empirical sample"
    biv.description = "Evaluable bivariate Gaussian reference parameters per DOY and grid cell; reference candidate only"
    biv.bivariate_pairs = json.dumps([list(p) for p in BIVARIATE_PAIRS])

    year_datasets = {year: Dataset(year_paths(year)[0], "r") for year in years}
    try:
        ny, nx = len(lat), len(lon)
        for doy0 in tqdm(range(DOY_COUNT), desc="Final DOY merge", unit="doy"):
            if doy0 == 58:
                continue
            for j0 in range(0, ny, 32):
                j1 = min(j0 + 32, ny)
                for i0 in range(0, nx, 64):
                    i1 = min(i0 + 64, nx)
                    shp = (j1 - j0) * (i1 - i0)
                    n = np.zeros(shp, np.int64)
                    total_sup = np.zeros(shp, np.int64)
                    total_bad = np.zeros(shp, np.int64)
                    states = {k: {f: np.zeros(shp, np.float64) for f in ("mean", "M2", "M3", "M4")} for k in ("rh", "e", "r", "q")}
                    pairs = {f"{x}__{y}": {"mx": np.zeros(shp, np.float64), "my": np.zeros(shp, np.float64), "c": np.zeros(shp, np.float64)} for x, y in BIVARIATE_PAIRS}

                    for year in years:
                        ds = year_datasets[year]
                        n2 = np.asarray(ds.variables["n_obs"][doy0, j0:j1, i0:i1], dtype=np.int64).reshape(-1)
                        n_before = n.copy()
                        for key in ("rh", "e", "r", "q"):
                            n, states[key]["mean"], states[key]["M2"], states[key]["M3"], states[key]["M4"] = combine_moments(
                                n_before, states[key]["mean"], states[key]["M2"], states[key]["M3"], states[key]["M4"],
                                n2,
                                np.asarray(ds.variables[f"mean_{key}"][doy0, j0:j1, i0:i1], dtype=np.float64).reshape(-1),
                                np.asarray(ds.variables[f"M2_{key}"][doy0, j0:j1, i0:i1], dtype=np.float64).reshape(-1),
                                np.asarray(ds.variables[f"M3_{key}"][doy0, j0:j1, i0:i1], dtype=np.float64).reshape(-1),
                                np.asarray(ds.variables[f"M4_{key}"][doy0, j0:j1, i0:i1], dtype=np.float64).reshape(-1),
                            )
                        total_sup += np.asarray(ds.variables["total_supersat"][doy0, j0:j1, i0:i1], dtype=np.int64).reshape(-1)
                        total_bad += np.asarray(ds.variables["total_invalid_ep"][doy0, j0:j1, i0:i1], dtype=np.int64).reshape(-1)
                        for xname, yname in BIVARIATE_PAIRS:
                            tag = f"{xname}__{yname}"
                            _, pairs[tag]["mx"], pairs[tag]["my"], pairs[tag]["c"] = combine_covariance(
                                n_before, pairs[tag]["mx"], pairs[tag]["my"], pairs[tag]["c"],
                                n2,
                                np.asarray(ds.variables[f"pair_{tag}_mean_x"][doy0, j0:j1, i0:i1], dtype=np.float64).reshape(-1),
                                np.asarray(ds.variables[f"pair_{tag}_mean_y"][doy0, j0:j1, i0:i1], dtype=np.float64).reshape(-1),
                                np.asarray(ds.variables[f"pair_{tag}_Cxy"][doy0, j0:j1, i0:i1], dtype=np.float64).reshape(-1),
                            )

                    for key in ("rh", "e", "r", "q"):
                        mean, std, skew, kurt = stats_from_state(n, states[key]["mean"], states[key]["M2"], states[key]["M3"], states[key]["M4"])
                        block = np.empty((j1 - j0, i1 - i0), dtype=np.float32)
                        for field, arr in (("mean", mean), ("std", std), ("skew", skew), ("kurt", kurt)):
                            block[:] = np.where(np.isfinite(arr), arr, -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)
                            main.variables[f"{field}_{key}"][doy0, j0:j1, i0:i1] = block

                    diag.variables["valid_observation_count"][doy0, j0:j1, i0:i1] = n.reshape(j1 - j0, i1 - i0)
                    denom = np.maximum(n, 1).astype(np.float64)
                    fs = total_sup / denom
                    fb = total_bad / denom
                    diag.variables["supersaturation_fraction"][doy0, j0:j1, i0:i1] = np.where(n > 0, fs, -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)
                    diag.variables["invalid_e_over_p_fraction"][doy0, j0:j1, i0:i1] = np.where(n > 0, fb, -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)

                    for xname, yname in BIVARIATE_PAIRS:
                        tag = f"{xname}__{yname}"
                        p = pairs[tag]
                        cov = np.full(shp, np.nan)
                        ok = n >= 2
                        cov[ok] = p["c"][ok] / (n[ok].astype(np.float64) - 1.0)
                        stdx = np.full(shp, np.nan); stdy = np.full(shp, np.nan)
                        stdx[ok] = np.sqrt(np.maximum(states[xname]["M2"][ok] / (n[ok] - 1), 0.0))
                        stdy[ok] = np.sqrt(np.maximum(states[yname]["M2"][ok] / (n[ok] - 1), 0.0))
                        corr = np.full(shp, np.nan)
                        okc = ok & (stdx > 0) & (stdy > 0)
                        corr[okc] = np.clip(cov[okc] / (stdx[okc] * stdy[okc]), -0.999999, 0.999999)
                        sl = np.s_[doy0, j0:j1, i0:i1]
                        biv.variables[f"{tag}_n"][sl] = n.reshape(j1 - j0, i1 - i0)
                        biv.variables[f"{tag}_mean_x"][sl] = np.where(np.isfinite(p["mx"]), p["mx"], -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)
                        biv.variables[f"{tag}_mean_y"][sl] = np.where(np.isfinite(p["my"]), p["my"], -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)
                        biv.variables[f"{tag}_std_x"][sl] = np.where(np.isfinite(stdx), stdx, -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)
                        biv.variables[f"{tag}_std_y"][sl] = np.where(np.isfinite(stdy), stdy, -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)
                        biv.variables[f"{tag}_cov"][sl] = np.where(np.isfinite(cov), cov, -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)
                        biv.variables[f"{tag}_corr"][sl] = np.where(np.isfinite(corr), corr, -9999.0).reshape(j1 - j0, i1 - i0).astype(np.float32)

            main.sync(); diag.sync(); biv.sync()
    finally:
        for ds in year_datasets.values():
            ds.close()
        main.close(); diag.close(); biv.close()

    manifest = {
        "project": "HumidClimatologyEngine",
        "implementation": "moisture_climatology_v8_0.py",
        "schema_version": SCHEMA_VERSION,
        "config_hash": CONFIG_HASH,
        "script_sha256": script_sha256(),
        "period": f"{START_YEAR}-{END_YEAR}",
        "calendar_slots": DOY_COUNT,
        "leap_day_policy": "pool_feb28_feb29_into_slot_60",
        "input_directories": {"t2m": str(T2M_DIR), "d2m": str(D2M_DIR), "sp": str(SP_DIR)},
        "bivariate_pairs": [list(p) for p in BIVARIATE_PAIRS],
        "outputs": {
            "main": {"path": str(OUTPUT_FILE), "sha256": sha256_file(OUTPUT_FILE)},
            "diagnostics": {"path": str(DIAGNOSTIC_FILE), "sha256": sha256_file(DIAGNOSTIC_FILE)},
            "bivariate_reference_parameters": {"path": str(BIVARIATE_FILE), "sha256": sha256_file(BIVARIATE_FILE)},
        },
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json_write(RUN_MANIFEST_FILE, manifest)

# =============================================================================
# 12. OPTIONAL SECOND-PASS EMPIRICAL BIVARIATE PDF
# =============================================================================

def _bivariate_output_path(pair: tuple[str,str]) -> Path:
    return OUTPUT_DIR / f"moisture_bivariate_empirical_{pair[0]}__{pair[1]}_1981_2020_v8_0.nc"

def _bivariate_pair_progress_path(pair: tuple[str,str]) -> Path:
    return CHECKPOINT_DIR / f"bivariate_{pair[0]}__{pair[1]}_progress_v8_0.json"

def _bivariate_progress_read(path: Path) -> dict:
    if not path.exists(): return {}
    try: return json.loads(path.read_text(encoding="utf-8"))
    except Exception: return {}

def build_empirical_bivariate_pair(pair: tuple[str, str], years: list[int], lat: np.ndarray, lon: np.ndarray) -> Path:
    """ساخت یا ادامه‌ی ساخت PDF دوبعدی تجربی برای یک جفت متغیر.

    از چک‌پوینت داخلی (next_year) برای ادامه از نقطه‌ی قطع استفاده می‌کند.
    """
    require_netcdf4()
    xname, yname = pair
    x_edges, y_edges = bivariate_histogram_edges(xname, yname)
    out_path = _bivariate_output_path(pair)
    ny, nx = len(lat), len(lon)
    y_ranges = list(range(0, ny, CHUNK_LAT))
    x_ranges = list(range(0, nx, CHUNK_LON))
    n_y_chunks, n_x_chunks = len(y_ranges), len(x_ranges)
    total_units = DOY_COUNT * n_y_chunks * n_x_chunks

    # باز کردن یا ایجاد فایل خروجی
    if out_path.exists():
        ds = Dataset(out_path, "r+")
    else:
        ds = Dataset(out_path, "w", format="NETCDF4")
        ds.createDimension("doy", DOY_COUNT)
        ds.createDimension("latitude", ny)
        ds.createDimension("longitude", nx)
        ds.createDimension("y_chunk", n_y_chunks)
        ds.createDimension("x_chunk", n_x_chunks)
        ds.createDimension("x_bin", BIVARIATE_NX)
        ds.createDimension("y_bin", BIVARIATE_NY)

        ds.createVariable("doy", "i2", ("doy",))[:] = np.arange(1, DOY_COUNT + 1, dtype=np.int16)
        ds.createVariable("latitude", "f4", ("latitude",))[:] = lat.astype(np.float32)
        ds.createVariable("longitude", "f4", ("longitude",))[:] = lon.astype(np.float32)
        ds.createVariable("x_bin_left", "f8", ("x_bin",))[:] = x_edges[:-1]
        ds.createVariable("x_bin_right", "f8", ("x_bin",))[:] = x_edges[1:]
        ds.createVariable("y_bin_left", "f8", ("y_bin",))[:] = y_edges[:-1]
        ds.createVariable("y_bin_right", "f8", ("y_bin",))[:] = y_edges[1:]

        chunks = (1, min(CHUNK_LAT, ny), min(CHUNK_LON, nx), BIVARIATE_NX, BIVARIATE_NY)
        ds.createVariable("count", "u2", ("doy", "latitude", "longitude", "x_bin", "y_bin"),
                          zlib=True, complevel=4, shuffle=True, chunksizes=chunks, fill_value=0)
        ds.createVariable("n_valid", "u2", ("doy", "latitude", "longitude"),
                          zlib=True, complevel=4, shuffle=True,
                          chunksizes=(1, min(CHUNK_LAT, ny), min(CHUNK_LON, nx)), fill_value=0)
        ds.createVariable("next_year", "i2", ("doy", "y_chunk", "x_chunk"),
                          zlib=True, complevel=4, fill_value=0)

        # متادیتا
        ds.purpose = "Empirical 2-D piecewise-constant PDF from hourly observations; no Gaussian assumption"
        ds.pair = f"{xname}__{yname}"
        ds.schema_version = SCHEMA_VERSION
        ds.config_hash = CONFIG_HASH
        ds.x_range = x_edges[[0, -1]].tolist()
        ds.y_range = y_edges[[0, -1]].tolist()
        ds.x_bin_count = BIVARIATE_NX
        ds.y_bin_count = BIVARIATE_NY
        ds.pdf_definition = "f(x,y)=count(xbin,ybin)/(N_valid * bin_area) within each bin"
        ds.restart_definition = "next_year[yday,y_chunk,x_chunk] is committed transactionally with count/n_valid"
        ds.sync()

    # ایندکس‌های سال‌ها
    indices = {year: (build_file_index(year, T2M_DIR),
                      build_file_index(year, D2M_DIR),
                      build_file_index(year, SP_DIR))
               for year in years}
    year_pos = {year: pos for pos, year in enumerate(years)}

    try:
        # خواندن وضعیت پیشرفت – با دقت نسبت به ماسک
        next_year_data = ds.variables["next_year"][:]
        # تبدیل آرایه ماسک‌دار به آرایه معمولی با پر کردن ۰ به جای ماسک
        next_year_filled = np.ma.filled(next_year_data, 0).astype(np.int16)
        done_units = int(np.count_nonzero(next_year_filled >= len(years)) if DOY_COUNT else 0)
        logger.info("BIVARIATE START | %s | empirical 2-D PDF | committed units %d/%d",
                    pair, done_units, total_units)

        # حلقه بر روی سال‌ها
        for year in years:
            yi = year_pos[year]
            t_idx, d_idx, p_idx = indices[year]

            for month in range(1, 13):
                with open_dataset(t_idx[month]) as dt, open_dataset(d_idx[month]) as dd, open_dataset(p_idx[month]) as dp:
                    dt = sort_dataset(dt)
                    dd = sort_dataset(dd)
                    dp = sort_dataset(dp)
                    validate_grids_and_axes(dt, dd, dp, year, month)

                    tu = dt["t2m"].attrs.get("units")
                    du = dd["d2m"].attrs.get("units")
                    pu = dp["sp"].attrs.get("units")

                    native_doys = dt.time.dt.dayofyear.values.astype(np.int16)
                    # لیست روزهای اقلیمی معتبر (به جز روز ۵۹)
                    cdoys = sorted({
                        get_clim_doy(int(d), year)
                        for d in native_doys
                        if get_clim_doy(int(d), year) not in (-1, 59)
                    })

                    for cdoy in cdoys:
                        time_idx = np.flatnonzero([
                            get_clim_doy(int(d), year) == cdoy
                            for d in native_doys
                        ])
                        if time_idx.size == 0:
                            continue

                        for yci, j0 in enumerate(y_ranges):
                            j1 = min(j0 + CHUNK_LAT, ny)
                            for xci, i0 in enumerate(x_ranges):
                                i1 = min(i0 + CHUNK_LON, nx)

                                # ====== اصلاح اصلی: خواندن امن next_year ======
                                raw_val = ds.variables["next_year"][cdoy - 1, yci, xci]
                                # اگر ماسک شده باشد یعنی هنوز پردازش نشده → مقدار ۰ در نظر گرفته می‌شود
                                current_progress = int(np.ma.filled(raw_val, 0))
                                if current_progress >= yi + 1:
                                    continue
                                # ==============================================

                                cells = (j1 - j0) * (i1 - i0)
                                counts = np.zeros((cells, BIVARIATE_NX, BIVARIATE_NY), dtype=np.uint32)
                                nvalid = np.zeros(cells, dtype=np.uint32)

                                for ti in time_idx:
                                    T = convert_temperature(
                                        dt["t2m"].isel(time=int(ti), latitude=slice(j0, j1), longitude=slice(i0, i1)).values,
                                        tu, "t2m"
                                    ).reshape(-1)
                                    Td = convert_temperature(
                                        dd["d2m"].isel(time=int(ti), latitude=slice(j0, j1), longitude=slice(i0, i1)).values,
                                        du, "d2m"
                                    ).reshape(-1)
                                    P = convert_pressure(
                                        dp["sp"].isel(time=int(ti), latitude=slice(j0, j1), longitude=slice(i0, i1)).values,
                                        pu, "sp"
                                    ).reshape(-1)

                                    ph = derive_moisture(T, Td, P)
                                    xv = np.asarray(ph[xname], dtype=np.float64).reshape(-1)
                                    yv = np.asarray(ph[yname], dtype=np.float64).reshape(-1)

                                    valid = np.isfinite(xv) & np.isfinite(yv)
                                    valid &= (xv >= x_edges[0]) & (xv <= x_edges[-1])
                                    valid &= (yv >= y_edges[0]) & (yv <= y_edges[-1])

                                    idx_valid = np.flatnonzero(valid)
                                    if idx_valid.size == 0:
                                        continue

                                    xx = xv[idx_valid]
                                    yy = yv[idx_valid]
                                    ix = np.searchsorted(x_edges, xx, side="right") - 1
                                    iy = np.searchsorted(y_edges, yy, side="right") - 1
                                    ix = np.minimum(ix, BIVARIATE_NX - 1)
                                    iy = np.minimum(iy, BIVARIATE_NY - 1)

                                    nvalid += np.bincount(
                                        idx_valid, minlength=nvalid.size
                                    ).astype(nvalid.dtype, copy=False)
                                    flat_idx = (
                                        idx_valid * (BIVARIATE_NX * BIVARIATE_NY)
                                        + ix * BIVARIATE_NY
                                        + iy
                                    )
                                    counts += np.bincount(
                                        flat_idx, minlength=counts.size
                                    ).reshape(counts.shape).astype(counts.dtype, copy=False)

                                # خواندن مقادیر قدیمی با دقت نسبت به ماسک
                                old_counts = np.asarray(
                                    ds.variables["count"][cdoy - 1, j0:j1, i0:i1, :, :],
                                    dtype=np.uint32
                                ).reshape(cells, BIVARIATE_NX, BIVARIATE_NY)
                                old_n = np.asarray(
                                    ds.variables["n_valid"][cdoy - 1, j0:j1, i0:i1],
                                    dtype=np.uint32
                                ).reshape(-1)

                                merged_counts = old_counts + counts
                                merged_n = old_n + nvalid

                                if merged_counts.max(initial=0) > np.iinfo(np.uint16).max or \
                                   merged_n.max(initial=0) > np.iinfo(np.uint16).max:
                                    raise OverflowError("Empirical bivariate histogram exceeds uint16 capacity.")

                                # ذخیره‌سازی تراکنشی (داده + وضعیت)
                                ds.variables["count"][cdoy - 1, j0:j1, i0:i1, :, :] = \
                                    merged_counts.astype(np.uint16).reshape(j1 - j0, i1 - i0, BIVARIATE_NX, BIVARIATE_NY)
                                ds.variables["n_valid"][cdoy - 1, j0:j1, i0:i1] = \
                                    merged_n.astype(np.uint16).reshape(j1 - j0, i1 - i0)
                                ds.variables["next_year"][cdoy - 1, yci, xci] = yi + 1
                                ds.sync()

                                # گزارش پیشرفت هر از گاهی
                                if (yi == len(years) - 1) and \
                                   ((cdoy % 10 == 0) or
                                    (cdoy == DOY_COUNT - 1 and xci == len(x_ranges) - 1 and yci == len(y_ranges) - 1)):
                                    # خواندن دوباره next_year با پر کردن ماسک
                                    next_year_filled = np.ma.filled(ds.variables["next_year"][:], 0).astype(np.int16)
                                    committed = int(np.count_nonzero(next_year_filled >= len(years)))
                                    pct = 100.0 * committed / max(total_units, 1)
                                    logger.info("BIVARIATE PROGRESS | %s | %.2f%% | %d/%d spatial-day units | remaining %d | through year %d",
                                                pair, pct, committed, total_units, total_units - committed, year)

        # بررسی تکمیل نهایی
        next_year_filled = np.ma.filled(ds.variables["next_year"][:], 0).astype(np.int16)
        committed = int(np.count_nonzero(next_year_filled >= len(years)))
        if committed != total_units:
            raise RuntimeError(f"Empirical bivariate checkpoint incomplete: {committed}/{total_units} units committed")
        logger.info("BIVARIATE COMPLETE | %s | 100.00%% | %d/%d units", pair, committed, total_units)

        return out_path

    finally:
        ds.close()

# =============================================================================
# 12. TESTS
# =============================================================================

def test_calendar() -> None:
    assert get_clim_doy(59, 1984) == 60
    assert get_clim_doy(60, 1984) == 60
    assert get_clim_doy(61, 1984) == 61
    assert get_clim_doy(59, 1985) == 60
    assert get_clim_doy(60, 1985) == 61
    assert get_clim_doy(365, 1985) == 366


def test_moments_against_numpy() -> None:
    xs = np.asarray([10.0, 20.0, 30.0, 40.0, 25.0])
    n = np.zeros(1, np.int64); m = np.zeros(1); M2 = np.zeros(1); M3 = np.zeros(1); M4 = np.zeros(1)
    for x in xs:
        update_moments_4_order(n, m, M2, M3, M4, np.asarray([x]), np.asarray([True]))
    assert n[0] == xs.size
    assert np.isclose(m[0], xs.mean())
    assert np.isclose(M2[0], np.sum((xs - xs.mean())**2))
    assert np.isclose(M3[0], np.sum((xs - xs.mean())**3))
    assert np.isclose(M4[0], np.sum((xs - xs.mean())**4))


def test_covariance_against_numpy() -> None:
    xs = np.asarray([1., 2., 3., 4.]); ys = np.asarray([2., 4., 5., 8.])
    n = np.zeros(1, np.int64); mx = np.zeros(1); my = np.zeros(1); c = np.zeros(1)
    for x, y in zip(xs, ys):
        update_covariance(n, mx, my, c, np.asarray([x]), np.asarray([y]), np.asarray([True]))
    assert n[0] == 4
    assert np.isclose(mx[0], xs.mean())
    assert np.isclose(my[0], ys.mean())
    assert np.isclose(c[0] / 3.0, np.cov(xs, ys, ddof=1)[0, 1])


def test_pebay_merge_equivalence() -> None:
    a = np.asarray([1., 2., 3.]); b = np.asarray([4., 7., 9.])
    def acc(xs):
        n = np.zeros(1, np.int64); m = np.zeros(1); M2 = np.zeros(1); M3 = np.zeros(1); M4 = np.zeros(1)
        for x in xs:
            update_moments_4_order(n, m, M2, M3, M4, np.asarray([x]), np.asarray([True]))
        return n, m, M2, M3, M4
    aa = acc(a); bb = acc(b)
    merged = combine_moments(*aa, *bb)
    full = acc(np.concatenate([a, b]))
    for u, v in zip(merged, full):
        assert np.allclose(u, v)


def test_physics() -> None:
    T = np.asarray([20.0, 0.0, -10.0], np.float32)
    Td = np.asarray([15.0, -1.0, -12.0], np.float32)
    P = np.asarray([1013.25, 900.0, 800.0], np.float32)
    d = derive_moisture(T, Td, P)
    assert np.all((d["rh"] >= 0) & (d["rh"] <= 100))
    assert np.all(d["e"] > 0)
    assert np.all(d["r"] >= 0)
    assert np.all((d["q"] >= 0) & (d["q"] < 1))


def test_bivariate_pdf() -> None:
    f0 = bivariate_gaussian_pdf(0., 0., 0., 1., 0., 1., 0.)
    assert np.isclose(f0, 1.0 / (2.0 * np.pi), rtol=1e-12, atol=1e-12)
    f = bivariate_gaussian_pdf(1., 1., 0., 1., 0., 1., 0.5)
    assert np.isfinite(f) and f > 0
    xe, ye = bivariate_histogram_edges("rh", "q")
    h,n = empirical_bivariate_histogram(np.array([50.,50.,100.]), np.array([0.2,0.2,0.9]), xe, ye)
    assert n == 3 and int(h.sum()) == 3
    pdf0 = empirical_bivariate_pdf(50.,0.2,h,n,xe,ye)
    assert np.isfinite(pdf0) and pdf0 > 0


def run_tests() -> None:
    logger.info("Running v8.0 FINAL SINGLE-PASS scientific/unit tests...")
    test_calendar(); test_moments_against_numpy(); test_covariance_against_numpy()
    test_pebay_merge_equivalence(); test_physics(); test_bivariate_pdf()
    for xname, yname in BIVARIATE_PAIRS:
        assert xname in {"rh", "e", "r", "q"} and yname in {"rh", "e", "r", "q"} and xname != yname
    logger.info("All v8.0 FINAL SINGLE-PASS tests passed.")

# =============================================================================
# 13. MAIN
# =============================================================================

def _global_progress(years: list[int]) -> tuple[int, int, float, int]:
    done_units = 0
    total_units = 0
    # Always include the full work plan, even before any progress JSON exists.
    # This prevents an initial misleading 0/0 global progress report.
    for year in years:
        final_path, json_path, part_path = year_paths(year)
        meta = _read_progress(json_path)
        if meta:
            done_units += int(meta.get("completed_units", 0))
            total_units += int(meta.get("total_units", 0))
        elif final_path.exists() and is_year_complete(year):
            # Completed years will be represented in metadata, but keep a safe fallback.
            with open_dataset(final_path) as ds:
                ny = ds.sizes["latitude"]; nx = ds.sizes["longitude"]
            total_units += 365 * ((ny + CHUNK_LAT - 1) // CHUNK_LAT) * ((nx + CHUNK_LON - 1) // CHUNK_LON)
            done_units += total_units
    remaining = max(total_units - done_units, 0)
    pct = 100.0 * done_units / max(total_units, 1)
    return done_units, total_units, pct, remaining


def main() -> None:
    logger.info("=" * 90)
    logger.info("MOISTURE CLIMATOLOGY v8.0 - PRODUCTION OPTIMIZED HOURLY EMPIRICAL ENGINE")
    logger.info("Period: %d-%d | Config hash: %s", START_YEAR, END_YEAR, CONFIG_HASH)
    logger.info("Bivariate empirical pairs: %s | bins=%dx%d | enabled=%s", BIVARIATE_PAIRS, BIVARIATE_NX, BIVARIATE_NY, BUILD_EMPIRICAL_BIVARIATE)
    logger.info("Spatial chunk: %dx%d cells | Workers: %d", CHUNK_LAT, CHUNK_LON, MAX_WORKERS)
    logger.info("Checkpoint: per-DOY/per-spatial-chunk completion bitmap; power-failure restart enabled")
    logger.info("Progress: explicit units, percent complete, remaining units, rate and ETA")
    logger.info("=" * 90)
    run_tests()

    sample_file = build_file_index(START_YEAR, T2M_DIR)[1]
    with open_dataset(sample_file) as ds:
        ds = sort_dataset(ds)
        if "t2m" not in ds.data_vars:
            raise RuntimeError("Expected variable 't2m' not found.")
        lat = ds.latitude.values
        lon = ds.longitude.values
    years = list(range(START_YEAR, END_YEAR + 1))
    ny, nx = len(lat), len(lon)
    per_year_units = 365 * ((ny + CHUNK_LAT - 1) // CHUNK_LAT) * ((nx + CHUNK_LON - 1) // CHUNK_LON)
    grand_total = per_year_units * len(years)
    logger.info("GLOBAL WORK PLAN | %d years | %d units/year | %d total units", len(years), per_year_units, grand_total)

    remaining_years = [y for y in years if not is_year_complete(y)]
    initial_done = grand_total - len(remaining_years) * per_year_units
    initial_remaining = grand_total - initial_done
    logger.info("GLOBAL PROGRESS | %.2f%% | %d/%d units | remaining %d | active years %d", 100.0 * initial_done / max(grand_total, 1), initial_done, grand_total, initial_remaining, len(remaining_years))
    if remaining_years:
        logger.info("Annual checkpoints to process: %d/%d", len(remaining_years), len(years))
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futures = {ex.submit(process_year_empirical, y): y for y in remaining_years}
            pending = set(futures)
            last_report = 0.0
            while pending:
                done = {f for f in pending if f.done()}
                if done:
                    for fut in done:
                        y = futures[fut]
                        result_y, path = fut.result()
                        if path is None or not is_year_complete(result_y):
                            raise RuntimeError(f"Year {result_y} processing failed; inspect retained checkpoint.")
                        pending.remove(fut)
                        du, tu, pct, rem = _global_progress(years)
                        logger.info("GLOBAL PROGRESS | %.2f%% | %d/%d units | remaining %d | completed years %d/%d", pct, du, tu, rem, len(years)-len(pending), len(years))
                now = time.time()
                if pending and now - last_report >= PROGRESS_REFRESH_SECONDS:
                    du, tu, pct, rem = _global_progress(years)
                    logger.info("GLOBAL PROGRESS | %.2f%% | %d/%d units | remaining %d | active years %d", pct, du, tu, rem, len(pending))
                    last_report = now
                if pending:
                    time.sleep(0.5)
    else:
        logger.info("All annual checkpoints already valid; no accumulation work remains.")

    bad = [y for y in years if not is_year_complete(y)]
    if bad:
        raise RuntimeError(f"Invalid/missing annual checkpoints: {bad}")

    du, tu, pct, rem = _global_progress(years)
    logger.info("ACCUMULATION COMPLETE | %.2f%% | %d/%d units | remaining %d", pct, du, tu, rem)
    logger.info("Finalizing outputs with DOY and spatial chunk streaming...")
    finalize_all(years, lat, lon)
    if BUILD_EMPIRICAL_BIVARIATE:
        for pair in BIVARIATE_PAIRS:
            if pair[0] in BIVARIATE_RANGES and pair[1] in BIVARIATE_RANGES:
                empirical_path = build_empirical_bivariate_pair(pair, years, lat, lon)
                logger.info("Empirical bivariate PDF: %s", empirical_path)
                manifest_now = {}
                if RUN_MANIFEST_FILE.exists():
                    try:
                        manifest_now = json.loads(RUN_MANIFEST_FILE.read_text(encoding="utf-8"))
                    except Exception:
                        manifest_now = {}
                manifest_now.setdefault("empirical_bivariate_pdfs", {})[f"{pair[0]}__{pair[1]}"] = {
                    "path": str(empirical_path),
                    "sha256": sha256_file(empirical_path),
                    "bins": [BIVARIATE_NX, BIVARIATE_NY],
                    "definition": "piecewise-constant empirical PDF from exact hourly bin counts",
                }
                manifest_now["updated_utc"] = datetime.now(timezone.utc).isoformat()
                atomic_json_write(RUN_MANIFEST_FILE, manifest_now)
            else:
                logger.warning("Skipping empirical 2-D PDF for %s: no fixed physical range configured", pair)
    logger.info("SUCCESS: v8.0 FINAL SINGLE-PASS completed.")
    logger.info("Main: %s", OUTPUT_FILE)
    logger.info("Diagnostics: %s", DIAGNOSTIC_FILE)
    logger.info("Bivariate: %s", BIVARIATE_FILE)
    logger.info("Run manifest: %s", RUN_MANIFEST_FILE)



# =============================================================================
# 13. FIVE-DAY CENTERED WINDOW + DISTRIBUTION/COPULA MODEL-SELECTION LAYER
# =============================================================================

WINDOW_HALF_WIDTH_DAYS = 2
WINDOW_SIZE_DAYS = 5
EDGE_PADDING_FILES = {
    # Optional combined T2m/D2m edge file supplied by the user for the first
    # two days needed by the 1981-01-01 centred window.
    "1980-12-30_to_1980-12-31": Path(r"K:\kazemi\papers\temperature_interpolation\19801230-19801231.nc"),
}

FIT_MIN_OBS = 30
FIT_AICC_DELTA_ACCEPT = 2.0
BIMODAL_N_INIT = 10
BIMODAL_MAX_ITER = 1000
BIMODAL_TOL = 1e-4
BIMODAL_REG_COVAR = 1e-6


def canonical_time_name(ds: xr.Dataset) -> str:
    """Return the actual time coordinate name; never infer it from filename."""
    for name in ("valid_time", "time"):
        if name in ds.coords or name in ds.dims:
            return name
    for name in ds.coords:
        try:
            if np.issubdtype(ds[name].dtype, np.datetime64):
                return name
        except TypeError:
            pass
    raise RuntimeError("No datetime-like time coordinate found in dataset")


def canonical_grid_coords(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    if "latitude" not in ds.coords or "longitude" not in ds.coords:
        raise RuntimeError("Dataset must contain latitude/longitude coordinates")
    return np.asarray(ds.latitude.values), np.asarray(ds.longitude.values)


def _month_file(folder: Path, year: int, month: int) -> Path:
    idx = build_file_index(year, folder)
    return idx[month]


def _special_edge_file() -> Optional[Path]:
    p = EDGE_PADDING_FILES["1980-12-30_to_1980-12-31"]
    return p if p.exists() else None


def _window_file_specs(start: datetime, end: datetime, folder: Path) -> list[tuple[Path, Optional[str]]]:
    """Return unique monthly files intersecting [start,end], plus the special edge file."""
    specs: list[tuple[Path, Optional[str]]] = []
    cur = datetime(start.year, start.month, 1)
    last = datetime(end.year, end.month, 1)
    while cur <= last:
        # For Dec-1980 the special combined edge file is authoritative for T/D.
        if cur.year == 1980 and cur.month == 12:
            edge = _special_edge_file()
            if edge is not None:
                specs.append((edge, "edge_combined"))
            cur = datetime(cur.year + (cur.month == 12), 1 if cur.month == 12 else cur.month + 1, 1)
            continue
        try:
            specs.append((_month_file(folder, cur.year, cur.month), None))
        except Exception:
            pass
        if cur.month == 12:
            cur = datetime(cur.year + 1, 1, 1)
        else:
            cur = datetime(cur.year, cur.month + 1, 1)
    # De-duplicate paths while preserving order.
    seen=set(); out=[]
    for item in specs:
        if item[0] not in seen:
            out.append(item); seen.add(item[0])
    return out


def _read_variable_time_slice(path: Path, variable_candidates: tuple[str, ...], start: np.datetime64, end: np.datetime64,
                              lat_slice: slice, lon_slice: slice) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Read one variable from one file using actual time coordinate; supports valid_time and combined edge files."""
    with xr.open_dataset(path, engine="netcdf4", decode_times=True, mask_and_scale=True, cache=False) as ds:
        tname=canonical_time_name(ds)
        var=next((v for v in variable_candidates if v in ds.data_vars), None)
        if var is None:
            raise KeyError(f"None of {variable_candidates} found in {path.name}")
        times=np.asarray(ds[tname].values)
        mask=(times>=start)&(times<=end)
        idx=np.flatnonzero(mask)
        if idx.size==0:
            return np.empty((0, lat_slice.stop-lat_slice.start, lon_slice.stop-lon_slice.start), np.float32), np.empty(0, dtype='datetime64[ns]'), np.asarray(ds.latitude.values), np.asarray(ds.longitude.values), ds[var].attrs.get("units", "")
        arr=ds[var].isel({tname: idx, "latitude":lat_slice, "longitude":lon_slice}).values
        return np.asarray(arr), times[idx], np.asarray(ds.latitude.values), np.asarray(ds.longitude.values), ds[var].attrs.get("units", "")


def extract_centered_five_day_window(target_date: str | datetime, variable: str, lat_index: int, lon_index: int,
                                    folder: Path, half_width_days: int = WINDOW_HALF_WIDTH_DAYS) -> tuple[np.ndarray, np.ndarray, dict]:
    """Extract one station/grid-cell hourly window using the project's 5-day centred rule.

    The target date is the centre; returned timestamps are the actual timestamps from files.
    The 1980-12-30/31 combined edge file is supported explicitly for the start of 1981.
    """
    target = datetime.fromisoformat(str(target_date)[:10])
    start = np.datetime64(target - timedelta(days=half_width_days))
    end = np.datetime64(target + timedelta(days=half_width_days, hours=23))
    arrays=[]; times_all=[]; units=None; lat_vals=None; lon_vals=None
    for path, mode in _window_file_specs(target - timedelta(days=half_width_days), target + timedelta(days=half_width_days), folder):
        with xr.open_dataset(path, engine="netcdf4", decode_times=True, mask_and_scale=True, cache=False) as ds:
            tname=canonical_time_name(ds)
            var = variable if variable in ds.data_vars else {"t2m":"t2m", "d2m":"d2m", "sp":"sp"}.get(variable, variable)
            if var not in ds.data_vars:
                # In the combined edge file, the requested variable may be absent (e.g. SP).
                continue
            times=np.asarray(ds[tname].values)
            mask=(times>=start)&(times<=end)
            idx=np.flatnonzero(mask)
            if idx.size:
                if lat_vals is None:
                    lat_vals=np.asarray(ds.latitude.values); lon_vals=np.asarray(ds.longitude.values)
                arr=ds[var].isel({tname:idx, "latitude":lat_index, "longitude":lon_index}).values.astype(np.float32)
                arrays.append(arr.reshape(-1)); times_all.append(times[idx]); units=ds[var].attrs.get("units", units)
    if arrays:
        order=np.argsort(np.concatenate(times_all))
        data=np.concatenate(arrays)[order]; times=np.concatenate(times_all)[order]
    else:
        data=np.empty(0,np.float32); times=np.empty(0,dtype='datetime64[ns]')
    expected=((end-start).astype('timedelta64[h]').astype(int)+1)
    meta={"target_date":target.strftime("%Y-%m-%d"), "window_start":str(start), "window_end":str(end),
          "expected_hour_count":int(expected), "available_hour_count":int(len(times)),
          "completeness_fraction":float(len(times)/max(expected,1)), "lat":None if lat_vals is None else float(lat_vals[lat_index]),
          "lon":None if lon_vals is None else float(lon_vals[lon_index]), "units":units}
    return data, times, meta


def extract_centered_five_day_moisture_window(target_date: str | datetime, lat_index: int, lon_index: int) -> tuple[dict[str, np.ndarray], dict]:
    """Extract T2m/D2m/SP for one grid cell over the 5-day centred window.

    Variables are aligned by actual timestamp. T2m/D2m can use the supplied combined
    1980-12-30/31 edge file. Surface pressure is read independently from SP_DIR; if
    those two edge days are unavailable there, RH/e can still be fitted from T/D while
    r/q/coplanarity carry their reduced pairwise coverage in diagnostics.
    """
    def gather(varname: str, folder: Path) -> tuple[np.ndarray, np.ndarray, dict]:
        aliases={"t2m":("t2m",),"d2m":("d2m",),"sp":("sp",)}
        vals=[]; ts=[]; units=None; latv=None; lonv=None
        target=datetime.fromisoformat(str(target_date)[:10])
        start=np.datetime64(target-timedelta(days=WINDOW_HALF_WIDTH_DAYS))
        end=np.datetime64(target+timedelta(days=WINDOW_HALF_WIDTH_DAYS, hours=23))
        for path, _ in _window_file_specs(target-timedelta(days=WINDOW_HALF_WIDTH_DAYS), target+timedelta(days=WINDOW_HALF_WIDTH_DAYS), folder):
            with xr.open_dataset(path, engine="netcdf4", decode_times=True, mask_and_scale=True, cache=False) as ds:
                tname=canonical_time_name(ds)
                var=next((v for v in aliases[varname] if v in ds.data_vars), None)
                if var is None:
                    continue
                times=np.asarray(ds[tname].values)
                m=(times>=start)&(times<=end); idx=np.flatnonzero(m)
                if not idx.size: continue
                if latv is None:
                    latv=np.asarray(ds.latitude.values); lonv=np.asarray(ds.longitude.values)
                a=ds[var].isel({tname:idx,"latitude":lat_index,"longitude":lon_index}).values.astype(np.float32).reshape(-1)
                vals.append(a); ts.append(times[idx]); units=ds[var].attrs.get("units",units)
        if not vals:
            return np.empty(0,np.float32), np.empty(0,dtype="datetime64[ns]"), {"units":units,"available":0,"expected":121,"completeness":0.0}
        tt=np.concatenate(ts); aa=np.concatenate(vals); order=np.argsort(tt); tt=tt[order]; aa=aa[order]
        expected=121
        return aa,tt,{"units":units,"available":int(tt.size),"expected":expected,"completeness":float(tt.size/expected)}
    out={}; meta={"target_date":str(target_date)[:10],"window_days":WINDOW_SIZE_DAYS,"variables":{}}
    for var,folder in (("t2m",T2M_DIR),("d2m",D2M_DIR),("sp",SP_DIR)):
        a,t,m=gather(var,folder); out[var]=a; meta["variables"][var]={**m,"first_time":None if t.size==0 else str(t[0]),"last_time":None if t.size==0 else str(t[-1])}
    # Exact timestamp intersection for paired states.
    if out["t2m"].size and out["d2m"].size:
        target=datetime.fromisoformat(str(target_date)[:10]); start=np.datetime64(target-timedelta(days=WINDOW_HALF_WIDTH_DAYS)); end=np.datetime64(target+timedelta(days=WINDOW_HALF_WIDTH_DAYS, hours=23))
        # Re-read timestamps using the generic single-variable extractor at the target grid point.
        _,tt,mt=extract_centered_five_day_window(target_date,"t2m",lat_index,lon_index,T2M_DIR)
        _,td,md=extract_centered_five_day_window(target_date,"d2m",lat_index,lon_index,D2M_DIR)
        common=np.intersect1d(tt,td); meta["paired_T_D_count"]=int(common.size)
    else:
        meta["paired_T_D_count"]=0
    return out,meta


def _aic_bic(loglik: float, n: int, k: int) -> tuple[float, float, float]:
    aic=2*k-2*loglik
    bic=k*np.log(max(n,1))-2*loglik
    if n <= k+1:
        aicc=np.inf
    else:
        aicc=aic+(2*k*(k+1))/(n-k-1)
    return float(aic), float(aicc), float(bic)


def _safe_logpdf_sum(dist, data, params) -> float:
    try:
        lp=np.asarray(dist.logpdf(data,*params),dtype=float)
        if not np.all(np.isfinite(lp)):
            return float("-inf")
        return float(lp.sum())
    except Exception:
        return float("-inf")


def _fit_standard_candidates(data: np.ndarray, bounded: bool = False) -> list[dict]:
    """Fit Normal, Skew-Normal, Pearson III and Beta where scientifically valid."""
    from scipy import stats
    x=np.asarray(data,dtype=float); x=x[np.isfinite(x)]
    if x.size < FIT_MIN_OBS:
        return []
    out=[]
    # Normal
    loc,scale=stats.norm.fit(x)
    if scale>0:
        ll=_safe_logpdf_sum(stats.norm,x,(loc,scale)); a,aicc,b=_aic_bic(ll,x.size,2)
        out.append({"name":"Normal","params":{"loc":float(loc),"scale":float(scale)},"loglik":ll,"aic":a,"aicc":aicc,"bic":b,"n_params":2})
    # Skew-Normal: retain both the project's moment-based fit and a likelihood-refined fit.
    try:
        proj=fit_skewnormal_project_style(x)
        if proj.get("status")=="ok":
            out.append({"name":"SkewNormal_ProjectStyle", **{k:v for k,v in proj.items() if k!="status"}})
    except Exception:
        pass
    try:
        alpha,loc,scale=stats.skewnorm.fit(x)
        if scale>0:
            ll=_safe_logpdf_sum(stats.skewnorm,x,(alpha,loc,scale)); a,aicc,b=_aic_bic(ll,x.size,3)
            out.append({"name":"SkewNormal_MLE","params":{"shape":float(alpha),"loc":float(loc),"scale":float(scale)},"loglik":ll,"aic":a,"aicc":aicc,"bic":b,"n_params":3})
    except Exception:
        pass
    # Pearson III
    try:
        skew,loc,scale=stats.pearson3.fit(x)
        if scale>0:
            ll=_safe_logpdf_sum(stats.pearson3,x,(skew,loc,scale)); a,aicc,b=_aic_bic(ll,x.size,3)
            out.append({"name":"PearsonIII","params":{"skew":float(skew),"loc":float(loc),"scale":float(scale)},"loglik":ll,"aic":a,"aicc":aicc,"bic":b,"n_params":3})
    except Exception:
        pass
    if bounded and np.nanmin(x) >= 0.0 and np.nanmax(x) <= 1.0 and np.nanmax(x) > np.nanmin(x):
        eps=max(1e-6, 0.5/max(x.size,1)); z=np.clip(x,eps,1-eps)
        try:
            aa,bb,loc,scale=stats.beta.fit(z,floc=0.0,fscale=1.0)
            ll=_safe_logpdf_sum(stats.beta,z,(aa,bb,loc,scale)); a,aicc,b=_aic_bic(ll,z.size,2)
            out.append({"name":"Beta","params":{"a":float(aa),"b":float(bb),"loc":0.0,"scale":1.0,"endpoint_epsilon":eps},"loglik":ll,"aic":a,"aicc":aicc,"bic":b,"n_params":2})
        except Exception:
            pass
    return out


def _fit_bimodal_normal(data: np.ndarray) -> Optional[dict]:
    from sklearn.mixture import GaussianMixture
    from scipy import stats
    x=np.asarray(data,dtype=float); x=x[np.isfinite(x)]
    if x.size < FIT_MIN_OBS:
        return None
    gm=GaussianMixture(n_components=2,n_init=BIMODAL_N_INIT,max_iter=BIMODAL_MAX_ITER,tol=BIMODAL_TOL,reg_covar=BIMODAL_REG_COVAR,random_state=20260821)
    gm.fit(x.reshape(-1,1))
    means=gm.means_.ravel(); scales=np.sqrt(gm.covariances_.ravel()); weights=gm.weights_.ravel()
    order=np.argsort(means); means=means[order]; scales=scales[order]; weights=weights[order]
    ll=float(gm.score(x.reshape(-1,1))*x.size); a,aicc,b=_aic_bic(ll,x.size,5)
    separation=abs(means[1]-means[0])/max(np.sqrt(0.5*(scales[0]**2+scales[1]**2)),1e-12)
    overlap=float(np.trapz(np.minimum(weights[0]*stats.norm.pdf(np.linspace(min(means)-4*scales.max(),max(means)+4*scales.max(),4000),means[0],scales[0]),weights[1]*stats.norm.pdf(np.linspace(min(means)-4*scales.max(),max(means)+4*scales.max(),4000),means[1],scales[1])),np.linspace(min(means)-4*scales.max(),max(means)+4*scales.max(),4000)))
    return {"name":"BimodalNormal","params":{"w1":float(weights[0]),"mu1":float(means[0]),"sigma1":float(scales[0]),"mu2":float(means[1]),"sigma2":float(scales[1])},"loglik":ll,"aic":a,"aicc":aicc,"bic":b,"n_params":5,"ashman_d":float(separation),"overlap_coefficient":overlap,"em_n_init":BIMODAL_N_INIT,"em_max_iter":BIMODAL_MAX_ITER,"em_tol":BIMODAL_TOL,"em_reg_covar":BIMODAL_REG_COVAR}


def fit_skewnormal_project_style(data: np.ndarray) -> dict:
    """Replicate the ClimateProcessingEngine plugin's moment-based parameterization.

    The public plugin uses sample mean/std/skew to derive (xi, omega, alpha), then
    evaluates the log-likelihood. We retain this as an auditable candidate alongside
    scipy's MLE fit instead of silently replacing the project's fitting convention.
    """
    from scipy import stats
    x=np.asarray(data,dtype=float); x=x[np.isfinite(x)]
    n=x.size
    if n < FIT_MIN_OBS:
        return {"status":"insufficient_data","n":int(n)}
    mean=float(np.mean(x)); std=float(np.std(x))
    if not np.isfinite(std) or std <= 0:
        return {"status":"degenerate","n":int(n)}
    g=float(np.mean(((x-mean)/std)**3))
    delta=g/np.sqrt(1.0+(2.0/np.pi-1.0)*g*g) if abs(g)<0.99 else 0.0
    omega=std/np.sqrt(max(1.0-2.0*delta*delta/np.pi,1e-12))
    xi=mean-omega*delta*np.sqrt(2.0/np.pi)
    alpha=delta/np.sqrt(max(1.0-delta*delta,1e-12)) if abs(delta)<0.999 else 0.0
    ll=float(np.sum(stats.skewnorm.logpdf(x,alpha,loc=xi,scale=omega)))
    a,aicc,b=_aic_bic(ll,n,3)
    return {"status":"ok","method":"ClimateProcessingEngine_plugin_moment_fit","params":{"shape":alpha,"loc":xi,"scale":omega},"loglik":ll,"aic":a,"aicc":aicc,"bic":b,"n_params":3,"sample_skew":g}


def fit_univariate_distribution_set(data: np.ndarray, variable: str) -> dict:
    bounded=variable in {"rh_fraction","q"}
    x=np.asarray(data,dtype=float); x=x[np.isfinite(x)]
    candidates=_fit_standard_candidates(x,bounded=bounded)
    bi=_fit_bimodal_normal(x)
    if bi is not None: candidates.append(bi)
    if not candidates:
        return {"status":"insufficient_data","n":int(x.size),"candidates":[],"best":None}
    candidates=sorted(candidates,key=lambda d:(d["aicc"],d["bic"]))
    best=candidates[0]
    return {"status":"ok","n":int(x.size),"candidates":candidates,"best":best,"selection_rule":"minimum AICc, then BIC; diagnostics retain multimodality metrics"}


def pseudo_observations(x: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata
    a=np.asarray(x,dtype=float); n=a.size
    r=rankdata(a,method="average")
    return (r-0.5)/n

def fit_gaussian_copula(x: np.ndarray, y: np.ndarray) -> dict:
    from scipy.special import ndtri
    mask=np.isfinite(x)&np.isfinite(y); x=np.asarray(x)[mask]; y=np.asarray(y)[mask]
    if x.size < FIT_MIN_OBS:
        return {"status":"insufficient_data","n":int(x.size)}
    u=pseudo_observations(x); v=pseudo_observations(y)
    zx=ndtri(np.clip(u,1e-10,1-1e-10)); zy=ndtri(np.clip(v,1e-10,1-1e-10))
    rho=float(np.corrcoef(zx,zy)[0,1]); rho=float(np.clip(rho,-0.9999,0.9999))
    den=np.sqrt(1-rho*rho); q=np.log1p(-rho*rho)
    ll=float(np.sum(-0.5*q-(rho*rho*(zx*zx+zy*zy)-2*rho*zx*zy)/(2*(1-rho*rho))))
    a,aicc,b=_aic_bic(ll,x.size,1)
    return {"status":"ok","n":int(x.size),"copula":"Gaussian","rho":rho,"loglik":ll,"aic":a,"aicc":aicc,"bic":b}


def fit_window_models(window_values: dict[str,np.ndarray]) -> dict:
    """Fit marginals and Gaussian copula to a single 5-day centred window."""
    result={"window":{k:int(np.isfinite(v).sum()) for k,v in window_values.items()}}
    result["marginals"]={}
    for var,vals in window_values.items():
        result["marginals"][var]=fit_univariate_distribution_set(vals,var)
    if "rh" in window_values and "q" in window_values:
        result["copula_rh_q"]=fit_gaussian_copula(window_values["rh"],window_values["q"])
    elif "rh_fraction" in window_values and "q" in window_values:
        result["copula_rh_q"]=fit_gaussian_copula(window_values["rh_fraction"],window_values["q"])
    return result


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        raise
    except Exception:
        logger.exception("Fatal error.")
        raise
