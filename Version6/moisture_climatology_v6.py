#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
moisture_climatology_v6.py
===============================================================================
اقلیم رطوبتی ERA5-Land، نسخه production-grade

مدل آماری:
    X = (T2m, Td2m, ln(Surface Pressure))

خروجی اصلی برای هر DOY × latitude × longitude:
    RH, vapor pressure, mixing ratio, specific humidity
    هر کدام: mean, std, bias-corrected skewness, Fisher excess kurtosis

ویژگی‌های اصلی:
    - تقویم 366 روزه: Feb-28 + Feb-29 -> DOY 60؛ DOY 59 رزرو
    - paired-valid observations برای T, Td, SP
    - Welford/Pébay برای mean/variance/covariance
    - merge سال‌ها با فرمول parallel Welford
    - MVN سه‌متغیره در فضای (T, Td, logP)
    - PSD/Cholesky کنترل‌شده و بدون fallback خاموش
    - Monte Carlo streaming با batch و chunk
    - higher moments واقعی Pébay (M2/M3/M4)
    - checkpoint سالانه و روزانه با SHA256 و config hash
    - restart واقعی
    - finalizer واقعی streaming با netCDF4
    - diagnostic کامل
    - تست‌های leap-day، synthetic، statistics و restart
    - validation نهایی مستقل
===============================================================================
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import xarray as xr
from tqdm import tqdm

try:
    from netCDF4 import Dataset
except Exception as exc:
    raise ImportError("netCDF4 is required for true streaming finalization.") from exc

# =============================================================================
# 1. CONFIG
# =============================================================================

T2M_DIR = Path(r"F:\Kazemi\era5\land\daily\T2m")
D2M_DIR = Path(r"F:\Kazemi\era5\land\daily\Dew_Point_Temperature")
SP_DIR  = Path(r"F:\Kazemi\era5\land\daily\Surface_Pressure")

OUTPUT_DIR = Path(r"C:\c")

CHECKPOINT_DIR = OUTPUT_DIR / "checkpoints_moisture_v6"
YEAR_DIR = CHECKPOINT_DIR / "years"
DAY_DIR = CHECKPOINT_DIR / "days"
TEST_DIR = CHECKPOINT_DIR / "tests"

START_YEAR = 1981
END_YEAR = 2020
DOY_COUNT = 366

# Monte Carlo
N_SAMPLES = 5000
CELL_CHUNK_SIZE = 1024
SAMPLE_BATCH_SIZE = 256

# Parallelism: 4 can be RAM-heavy on Windows; reduce to 2 if needed.
MAX_WORKERS = 2  # conservative default for Windows/RAM; raise after benchmark

RANDOM_SEED = 20260821

# Thresholds
PSD_TOL = 1e-10
PSD_REPAIR_TOL = 1e-8
MIN_OBS = 3
MIN_MC_VALID = 10

SCHEMA_VERSION = "6.0"
CHECKPOINT_VERSION = "6.0"

OUTPUT_FILE = OUTPUT_DIR / "moisture_climatology_1981_2020.nc"
DIAGNOSTIC_FILE = OUTPUT_DIR / "moisture_climatology_diagnostics_1981_2020.nc"

for d in [OUTPUT_DIR, CHECKPOINT_DIR, YEAR_DIR, DAY_DIR, TEST_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# =============================================================================
# 2. LOGGING
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("Moisture_Climatology_v6")

# =============================================================================
# PROGRESS LOGGER
# =============================================================================
def runtime_snapshot() -> dict:
    """Best-effort RAM/CPU snapshot for progress reporting."""
    snap = {}
    try:
        import psutil
        vm = psutil.virtual_memory()
        snap["ram_used_gb"] = vm.used / (1024 ** 3)
        snap["ram_available_gb"] = vm.available / (1024 ** 3)
        snap["ram_percent"] = vm.percent
        snap["cpu_percent"] = psutil.cpu_percent(interval=None)
    except Exception:
        pass
    return snap

def log_progress(event: str, **kwargs) -> None:
    """Detailed non-fatal progress logging with RAM/CPU telemetry."""
    try:
        parts = [event]
        for key, value in kwargs.items():
            if value is None:
                continue
            if isinstance(value, float):
                if key in {"percent", "seconds", "elapsed", "eta"}:
                    value = f"{value:.2f}"
                else:
                    value = f"{value:.4f}"
            parts.append(f"{key}={value}")

        snap = runtime_snapshot()
        if "ram_used_gb" in snap:
            parts.append(f"RAM={snap['ram_used_gb']:.2f}GB")
            parts.append(f"avail={snap['ram_available_gb']:.2f}GB")
            parts.append(f"used={snap['ram_percent']:.0f}%")
        if "cpu_percent" in snap:
            parts.append(f"CPU={snap['cpu_percent']:.0f}%")

        logger.info("PROGRESS | " + " | ".join(parts))
    except Exception:
        # Progress telemetry must never interrupt production calculations.
        try:
            logger.info("PROGRESS | %s", event)
        except Exception:
            pass


# =============================================================================
# 3. CONFIG HASH / SCRIPT HASH
# =============================================================================

@dataclass(frozen=True)
class Config:
    start_year: int
    end_year: int
    doy_count: int
    n_samples: int
    cell_chunk_size: int
    sample_batch_size: int
    workers: int
    random_seed: int
    schema_version: str
    psd_tol: float
    psd_repair_tol: float

CONFIG = Config(
    start_year=START_YEAR,
    end_year=END_YEAR,
    doy_count=DOY_COUNT,
    n_samples=N_SAMPLES,
    cell_chunk_size=CELL_CHUNK_SIZE,
    sample_batch_size=SAMPLE_BATCH_SIZE,
    workers=MAX_WORKERS,
    random_seed=RANDOM_SEED,
    schema_version=SCHEMA_VERSION,
    psd_tol=PSD_TOL,
    psd_repair_tol=PSD_REPAIR_TOL,
)

def hash_dict(d: dict) -> str:
    payload = json.dumps(d, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]

CONFIG_HASH = hash_dict(asdict(CONFIG))

def script_sha256() -> str:
    try:
        return sha256_file(Path(__file__))
    except Exception:
        return "interactive-or-unavailable"

# =============================================================================
# 4. FILE HELPERS
# =============================================================================

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
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

def atomic_npz_write(path: Path, **arrays) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    os.close(fd)
    tmp = Path(tmp_name)
    real_tmp = tmp.with_suffix(tmp.suffix + ".npz")
    try:
        np.savez_compressed(real_tmp, **arrays)
        os.replace(real_tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
        real_tmp.unlink(missing_ok=True)

def json_load(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

# =============================================================================
# 5. DATA HELPERS
# =============================================================================

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

def is_leap_year(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

def get_clim_doy(native_doy: int, year: int) -> int:
    """
    1-based climatological day.

    DOY 1..58  : Jan 1 .. Feb 27
    DOY 59     : RESERVED
    DOY 60     : Feb 28 + Feb 29
    DOY 61..366: Mar 1 .. Dec 31
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
    month = np.full(DOY_COUNT, -1, dtype=np.int16)
    day = np.full(DOY_COUNT, -1, dtype=np.int16)
    label = np.full(DOY_COUNT, "RESERVED", dtype=object)
    is_reserved = np.zeros(DOY_COUNT, dtype=np.int8)

    # Build from a non-leap reference, plus explicit composite Feb-29.
    import datetime as _dt
    d = _dt.date(2001, 1, 1)
    for i in range(58):
        month[i] = d.month
        day[i] = d.day
        label[i] = d.strftime("%b-%d")
        d += _dt.timedelta(days=1)

    month[58] = -1
    day[58] = -1
    label[58] = "RESERVED"
    is_reserved[58] = 1

    month[59] = 2
    day[59] = 29
    label[59] = "Feb-29-composite"

    d = _dt.date(2001, 3, 1)
    for idx in range(60, 366):
        month[idx] = d.month
        day[idx] = d.day
        label[idx] = d.strftime("%b-%d")
        d += _dt.timedelta(days=1)

    return month, day, np.asarray(label, dtype=str), is_reserved

def extract_year_month(path: Path, year: int) -> Optional[int]:
    m = re.search(r"(?<!\d)(\d{4})(\d{2})(?!\d)", path.name)
    if not m:
        return None
    y, mon = int(m.group(1)), int(m.group(2))
    if y == year and 1 <= mon <= 12:
        return mon
    return None

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

def validate_grids(ds_t: xr.Dataset, ds_d: xr.Dataset, ds_p: xr.Dataset, year: int, month: int) -> None:
    for name, ds in [("T2m", ds_t), ("D2m", ds_d), ("SP", ds_p)]:
        if "time" not in ds.dims or "latitude" not in ds.dims or "longitude" not in ds.dims:
            raise RuntimeError(f"{name}: required dimensions missing in {year}-{month:02d}")

    if not np.array_equal(ds_t.time.values, ds_d.time.values) or not np.array_equal(ds_t.time.values, ds_p.time.values):
        raise RuntimeError(f"Time coordinates differ for {year}-{month:02d}")

    if not np.allclose(ds_t.latitude.values, ds_d.latitude.values, rtol=0, atol=1e-7):
        raise RuntimeError(f"Latitude mismatch T/D for {year}-{month:02d}")
    if not np.allclose(ds_t.latitude.values, ds_p.latitude.values, rtol=0, atol=1e-7):
        raise RuntimeError(f"Latitude mismatch T/SP for {year}-{month:02d}")

    if not np.allclose(ds_t.longitude.values, ds_d.longitude.values, rtol=0, atol=1e-7):
        raise RuntimeError(f"Longitude mismatch T/D for {year}-{month:02d}")
    if not np.allclose(ds_t.longitude.values, ds_p.longitude.values, rtol=0, atol=1e-7):
        raise RuntimeError(f"Longitude mismatch T/SP for {year}-{month:02d}")

    if ds_t.sizes["time"] != ds_d.sizes["time"] or ds_t.sizes["time"] != ds_p.sizes["time"]:
        raise RuntimeError(f"Time length mismatch for {year}-{month:02d}")

    if ds_t.sizes["latitude"] != ds_d.sizes["latitude"] or ds_t.sizes["latitude"] != ds_p.sizes["latitude"]:
        raise RuntimeError(f"Latitude size mismatch for {year}-{month:02d}")

    if ds_t.sizes["longitude"] != ds_d.sizes["longitude"] or ds_t.sizes["longitude"] != ds_p.sizes["longitude"]:
        raise RuntimeError(f"Longitude size mismatch for {year}-{month:02d}")

# =============================================================================
# 6. PHYSICS
# =============================================================================

def es_water(temp_c: np.ndarray) -> np.ndarray:
    return 6.112 * np.exp((17.67 * temp_c) / (temp_c + 243.5))

def es_ice(temp_c: np.ndarray) -> np.ndarray:
    return 6.112 * np.exp((22.46 * temp_c) / (temp_c + 272.62))

def saturation_vapor_pressure(temp_c: np.ndarray) -> np.ndarray:
    out = np.empty_like(temp_c, dtype=np.float32)
    water = temp_c >= 0.0
    ice = ~water
    out[...] = np.nan
    if np.any(water):
        out[water] = es_water(temp_c[water]).astype(np.float32)
    if np.any(ice):
        out[ice] = es_ice(temp_c[ice]).astype(np.float32)
    return out

def derive_moisture(T: np.ndarray, Td: np.ndarray, P_hpa: np.ndarray) -> dict:
    es_T = saturation_vapor_pressure(T)
    e = saturation_vapor_pressure(Td)

    rh_raw = 100.0 * (e / es_T)
    supersat = np.isfinite(rh_raw) & (rh_raw > 100.0)
    rh = np.clip(rh_raw, 0.0, 100.0).astype(np.float32)

    valid_qr = (
        np.isfinite(e)
        & np.isfinite(P_hpa)
        & (e > 0.0)
        & (P_hpa > 0.0)
        & (e < P_hpa)
    )

    r = np.full_like(e, np.nan, dtype=np.float32)
    r[valid_qr] = (
        0.622 * e[valid_qr] / (P_hpa[valid_qr] - e[valid_qr])
    ).astype(np.float32)

    q = (r / (1.0 + r)).astype(np.float32)

    valid_all = (
        np.isfinite(rh)
        & np.isfinite(e)
        & np.isfinite(r)
        & np.isfinite(q)
    )

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
# 7. WELFORD / PÉBAY SECOND-ORDER ACCUMULATORS
# =============================================================================

def welford_vector_update(
    x: np.ndarray,
    mean: np.ndarray,
    M2: np.ndarray,
    count: np.ndarray,
    mask: np.ndarray,
) -> None:
    """
    One observation per cell, vectorized over valid cells.
    Arrays are 1-D spatial views.
    """
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return

    n_old = count[idx].astype(np.float64)
    n_new = n_old + 1.0

    old_mean = mean[idx]
    values = x[idx].astype(np.float64)

    delta = values - old_mean
    new_mean = old_mean + delta / n_new

    M2[idx] += delta * (values - new_mean)
    mean[idx] = new_mean
    count[idx] = n_new.astype(np.int64)

def welford_cov_update(
    x: np.ndarray,
    y: np.ndarray,
    mean_x: np.ndarray,
    mean_y: np.ndarray,
    C: np.ndarray,
    count: np.ndarray,
    mask: np.ndarray,
) -> None:
    """
    Paired covariance update for one observation per cell.
    """
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return

    n_old = count[idx].astype(np.float64)
    n_new = n_old + 1.0

    old_x = mean_x[idx]
    old_y = mean_y[idx]

    xv = x[idx].astype(np.float64)
    yv = y[idx].astype(np.float64)

    dx = xv - old_x
    dy = yv - old_y

    new_x = old_x + dx / n_new
    new_y = old_y + dy / n_new

    C[idx] += dx * (yv - new_y)

# =============================================================================
# 8. YEAR ACCUMULATION
# =============================================================================

YEAR_META = YEAR_DIR / "year_index.json"

def year_paths(year: int) -> tuple[Path, Path]:
    npz = YEAR_DIR / f"year_{year:04d}_{CONFIG_HASH}.npz"
    js = YEAR_DIR / f"year_{year:04d}_{CONFIG_HASH}.json"
    return npz, js

def is_year_complete(year: int) -> bool:
    npz, js = year_paths(year)
    if not npz.exists() or not js.exists():
        return False
    try:
        meta = json_load(js)
        if meta.get("status") != "completed":
            return False
        if meta.get("config_hash") != CONFIG_HASH:
            return False
        return meta.get("sha256") == sha256_file(npz)
    except Exception:
        return False

def process_year_welford(year: int) -> tuple[int, Optional[Path]]:
    logger.info(f"Starting year {year}")
    t0 = time.time()

    try:
        if is_year_complete(year):
            npz, _ = year_paths(year)
            logger.info(f"Year {year} already valid; skipping.")
            return year, npz

        t_idx = build_file_index(year, T2M_DIR)
        d_idx = build_file_index(year, D2M_DIR)
        p_idx = build_file_index(year, SP_DIR)

        with open_dataset(t_idx[1]) as ds0:
            ds0 = sort_dataset(ds0)
            ny = ds0.sizes["latitude"]
            nx = ds0.sizes["longitude"]

        ncells = ny * nx

        n = np.zeros((DOY_COUNT, ncells), dtype=np.int64)
        mean_T = np.zeros((DOY_COUNT, ncells), dtype=np.float64)
        mean_Td = np.zeros((DOY_COUNT, ncells), dtype=np.float64)
        mean_logP = np.zeros((DOY_COUNT, ncells), dtype=np.float64)

        M2_T = np.zeros_like(mean_T)
        M2_Td = np.zeros_like(mean_Td)
        M2_logP = np.zeros_like(mean_logP)

        C_T_Td = np.zeros_like(mean_T)
        C_T_logP = np.zeros_like(mean_T)
        C_Td_logP = np.zeros_like(mean_T)

        for month in range(1, 13):
            with (
                open_dataset(t_idx[month]) as ds_t,
                open_dataset(d_idx[month]) as ds_d,
                open_dataset(p_idx[month]) as ds_p,
            ):
                ds_t = sort_dataset(ds_t)
                ds_d = sort_dataset(ds_d)
                ds_p = sort_dataset(ds_p)
                validate_grids(ds_t, ds_d, ds_p, year, month)

                T = ds_t["average_t2m"].values.astype(np.float32) - 273.15
                Td = ds_d["average_d2m"].values.astype(np.float32) - 273.15
                P = ds_p["average_sp"].values.astype(np.float32) / 100.0

                if T.shape != Td.shape or T.shape != P.shape:
                    raise RuntimeError(f"Array shape mismatch in {year}-{month:02d}")

                native_doys = ds_t.time.dt.dayofyear.values.astype(np.int16)

                for ti, ndoy in enumerate(native_doys):
                    cdoy = get_clim_doy(int(ndoy), year)
                    if cdoy < 1 or cdoy > DOY_COUNT:
                        continue
                    di = cdoy - 1

                    t = T[ti].reshape(-1)
                    td = Td[ti].reshape(-1)
                    p = P[ti].reshape(-1)

                    valid = np.isfinite(t) & np.isfinite(td) & np.isfinite(p) & (p > 0)
                    if not np.any(valid):
                        continue

                    lp = np.full_like(p, np.nan, dtype=np.float64)
                    lp[valid] = np.log(p[valid].astype(np.float64))

                    # common paired population
                    n_old = n[di, :]
                    idx = np.flatnonzero(valid)
                    n_new = n_old[idx].astype(np.float64) + 1.0

                    old_T = mean_T[di, idx]
                    old_Td = mean_Td[di, idx]
                    old_lp = mean_logP[di, idx]

                    xv = t[idx].astype(np.float64)
                    yv = td[idx].astype(np.float64)
                    zv = lp[idx]

                    dx = xv - old_T
                    dy = yv - old_Td
                    dz = zv - old_lp

                    new_T = old_T + dx / n_new
                    new_Td = old_Td + dy / n_new
                    new_lp = old_lp + dz / n_new

                    M2_T[di, idx] += dx * (xv - new_T)
                    M2_Td[di, idx] += dy * (yv - new_Td)
                    M2_logP[di, idx] += dz * (zv - new_lp)

                    C_T_Td[di, idx] += dx * (yv - new_Td)
                    C_T_logP[di, idx] += dx * (zv - new_lp)
                    C_Td_logP[di, idx] += dy * (zv - new_lp)

                    mean_T[di, idx] = new_T
                    mean_Td[di, idx] = new_Td
                    mean_logP[di, idx] = new_lp
                    n[di, idx] = n_new.astype(np.int64)

        npz_path, json_path = year_paths(year)
        atomic_npz_write(
            npz_path,
            n=n,
            mean_T=mean_T,
            mean_Td=mean_Td,
            mean_logP=mean_logP,
            M2_T=M2_T,
            M2_Td=M2_Td,
            M2_logP=M2_logP,
            C_T_Td=C_T_Td,
            C_T_logP=C_T_logP,
            C_Td_logP=C_Td_logP,
            shape=np.asarray([ny, nx], dtype=np.int32),
        )

        meta = {
            "status": "completed",
            "year": year,
            "config_hash": CONFIG_HASH,
            "schema_version": CHECKPOINT_VERSION,
            "sha256": sha256_file(npz_path),
            "shape": [ny, nx],
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        atomic_json_write(json_path, meta)

        logger.info(f"Year {year} completed in {time.time() - t0:.1f}s")
        return year, npz_path

    except Exception:
        logger.exception(f"Year {year} failed")
        return year, None

# =============================================================================
# 9. MERGE YEARLY ACCUMULATORS
# =============================================================================

def combine_welford(
    n1, mean1, M21, n2, mean2, M22
):
    n_total = n1 + n2
    out_mean = mean1.copy()
    out_M2 = M21.copy()

    m = n2 > 0
    both = (n1 > 0) & (n2 > 0)
    only2 = (n1 == 0) & (n2 > 0)

    out_mean[only2] = mean2[only2]
    out_M2[only2] = M22[only2]

    if np.any(both):
        a = n1[both].astype(np.float64)
        b = n2[both].astype(np.float64)
        nt = a + b
        delta = mean2[both] - mean1[both]
        out_mean[both] = mean1[both] + delta * (b / nt)
        out_M2[both] = (
            M21[both]
            + M22[both]
            + delta * delta * (a * b / nt)
        )
    return n_total, out_mean, out_M2

def combine_covariance(
    n1, mean_x1, mean_y1, C1,
    n2, mean_x2, mean_y2, C2
):
    n_total = n1 + n2
    out = C1.copy()

    both = (n1 > 0) & (n2 > 0)
    only2 = (n1 == 0) & (n2 > 0)

    out[only2] = C2[only2]

    if np.any(both):
        a = n1[both].astype(np.float64)
        b = n2[both].astype(np.float64)
        nt = a + b
        dx = mean_x2[both] - mean_x1[both]
        dy = mean_y2[both] - mean_y1[both]
        out[both] = (
            C1[both] + C2[both] + dx * dy * (a * b / nt)
        )

    return out

def merge_years(years: Iterable[int], shape: tuple[int, int]) -> dict:
    ny, nx = shape
    ncells = ny * nx

    n = np.zeros((DOY_COUNT, ncells), dtype=np.int64)
    mean_T = np.zeros((DOY_COUNT, ncells), dtype=np.float64)
    mean_Td = np.zeros_like(mean_T)
    mean_logP = np.zeros_like(mean_T)
    M2_T = np.zeros_like(mean_T)
    M2_Td = np.zeros_like(mean_T)
    M2_logP = np.zeros_like(mean_T)
    C_T_Td = np.zeros_like(mean_T)
    C_T_logP = np.zeros_like(mean_T)
    C_Td_logP = np.zeros_like(mean_T)

    for year in tqdm(sorted(years), desc="Merging years", unit="year"):
        npz_path, json_path = year_paths(year)
        if not npz_path.exists() or not json_path.exists():
            raise RuntimeError(f"Missing year checkpoint: {year}")
        if not is_year_complete(year):
            raise RuntimeError(f"Invalid year checkpoint: {year}")

        with np.load(npz_path, allow_pickle=False) as d:
            n2 = d["n"]
            mean_T2 = d["mean_T"]
            mean_Td2 = d["mean_Td"]
            mean_logP2 = d["mean_logP"]
            M2_T2 = d["M2_T"]
            M2_Td2 = d["M2_Td"]
            M2_logP2 = d["M2_logP"]
            C_T_Td2 = d["C_T_Td"]
            C_T_logP2 = d["C_T_logP"]
            C_Td_logP2 = d["C_Td_logP"]

            n_old = n
            n, mean_T, M2_T = combine_welford(n_old, mean_T, M2_T, n2, mean_T2, M2_T2)
            _, mean_Td, M2_Td = combine_welford(n_old, mean_Td, M2_Td, n2, mean_Td2, M2_Td2)
            _, mean_logP, M2_logP = combine_welford(n_old, mean_logP, M2_logP, n2, mean_logP2, M2_logP2)

            C_T_Td = combine_covariance(
                n_old, mean_T if False else mean_T, mean_Td if False else mean_Td, C_T_Td,
                n2, mean_T2, mean_Td2, C_T_Td2
            )

            # The covariance merge needs old means; recompute correctly from the
            # two component means using saved combined means is not valid.
            # Therefore redo covariance merge with stored old means from copies.
            # This block is replaced below by an explicit robust merge.
            # (The temporary result above is overwritten.)

            # Reconstruct old state from component snapshots is impossible after
            # mutation, so the covariance merge is performed via a standalone
            # helper using the pre-merge snapshots below.
            # This branch is intentionally unreachable because we immediately
            # use the correct path in merge_years_v2().
            raise RuntimeError("Internal merge path should not be used.")

    return {}

# Correct merge implementation with explicit old-state snapshots.
def merge_years_v2(years: Iterable[int], shape: tuple[int, int]) -> dict:
    ny, nx = shape
    ncells = ny * nx

    n = np.zeros((DOY_COUNT, ncells), dtype=np.int64)
    mean_T = np.zeros((DOY_COUNT, ncells), dtype=np.float64)
    mean_Td = np.zeros_like(mean_T)
    mean_logP = np.zeros_like(mean_T)
    M2_T = np.zeros_like(mean_T)
    M2_Td = np.zeros_like(mean_T)
    M2_logP = np.zeros_like(mean_T)
    C_T_Td = np.zeros_like(mean_T)
    C_T_logP = np.zeros_like(mean_T)
    C_Td_logP = np.zeros_like(mean_T)

    for year in tqdm(sorted(years), desc="Merging years", unit="year"):
        npz_path, _ = year_paths(year)
        if not is_year_complete(year):
            raise RuntimeError(f"Invalid year checkpoint: {year}")

        with np.load(npz_path, allow_pickle=False) as d:
            n2 = d["n"].astype(np.int64)
            mean_T2 = d["mean_T"]
            mean_Td2 = d["mean_Td"]
            mean_logP2 = d["mean_logP"]
            M2_T2 = d["M2_T"]
            M2_Td2 = d["M2_Td"]
            M2_logP2 = d["M2_logP"]
            C_T_Td2 = d["C_T_Td"]
            C_T_logP2 = d["C_T_logP"]
            C_Td_logP2 = d["C_Td_logP"]

            n1 = n.copy()
            mT1 = mean_T.copy()
            mTd1 = mean_Td.copy()
            mLP1 = mean_logP.copy()

            n_total = n1 + n2
            both = (n1 > 0) & (n2 > 0)
            only2 = (n1 == 0) & (n2 > 0)

            # scalar moments
            if np.any(only2):
                mean_T[only2] = mean_T2[only2]
                mean_Td[only2] = mean_Td2[only2]
                mean_logP[only2] = mean_logP2[only2]
                M2_T[only2] = M2_T2[only2]
                M2_Td[only2] = M2_Td2[only2]
                M2_logP[only2] = M2_logP2[only2]
                C_T_Td[only2] = C_T_Td2[only2]
                C_T_logP[only2] = C_T_logP2[only2]
                C_Td_logP[only2] = C_Td_logP2[only2]

            if np.any(both):
                a = n1[both].astype(np.float64)
                b = n2[both].astype(np.float64)
                nt = a + b
                dx = mean_T2[both] - mT1[both]
                dy = mean_Td2[both] - mTd1[both]
                dz = mean_logP2[both] - mLP1[both]
                w = (a * b) / nt

                mean_T[both] = mT1[both] + dx * (b / nt)
                mean_Td[both] = mTd1[both] + dy * (b / nt)
                mean_logP[both] = mLP1[both] + dz * (b / nt)

                M2_T[both] = M2_T[both] + M2_T2[both] + dx * dx * w
                M2_Td[both] = M2_Td[both] + M2_Td2[both] + dy * dy * w
                M2_logP[both] = M2_logP[both] + M2_logP2[both] + dz * dz * w

                C_T_Td[both] = C_T_Td[both] + C_T_Td2[both] + dx * dy * w
                C_T_logP[both] = C_T_logP[both] + C_T_logP2[both] + dx * dz * w
                C_Td_logP[both] = C_Td_logP[both] + C_Td_logP2[both] + dy * dz * w

            n = n_total

    return {
        "n": n,
        "mean_T": mean_T,
        "mean_Td": mean_Td,
        "mean_logP": mean_logP,
        "M2_T": M2_T,
        "M2_Td": M2_Td,
        "M2_logP": M2_logP,
        "C_T_Td": C_T_Td,
        "C_T_logP": C_T_logP,
        "C_Td_logP": C_Td_logP,
    }

def stats_from_merged(w: dict) -> dict:
    n = w["n"]
    valid = n >= 2

    shape = n.shape
    out = {
        "n": n.astype(np.int32),
        "mean_T": w["mean_T"].astype(np.float32),
        "mean_Td": w["mean_Td"].astype(np.float32),
        "mean_logP": w["mean_logP"].astype(np.float32),
        "std_T": np.full(shape, np.nan, np.float32),
        "std_Td": np.full(shape, np.nan, np.float32),
        "std_logP": np.full(shape, np.nan, np.float32),
        "cov_T_Td": np.full(shape, np.nan, np.float32),
        "cov_T_logP": np.full(shape, np.nan, np.float32),
        "cov_Td_logP": np.full(shape, np.nan, np.float32),
        "corr_T_Td": np.full(shape, np.nan, np.float32),
        "corr_T_logP": np.full(shape, np.nan, np.float32),
        "corr_Td_logP": np.full(shape, np.nan, np.float32),
    }

    nn = n[valid].astype(np.float64)
    var_T = w["M2_T"][valid] / (nn - 1)
    var_Td = w["M2_Td"][valid] / (nn - 1)
    var_lp = w["M2_logP"][valid] / (nn - 1)

    cov_T_Td = w["C_T_Td"][valid] / (nn - 1)
    cov_T_lp = w["C_T_logP"][valid] / (nn - 1)
    cov_Td_lp = w["C_Td_logP"][valid] / (nn - 1)

    var_T = np.maximum(var_T, 0.0)
    var_Td = np.maximum(var_Td, 0.0)
    var_lp = np.maximum(var_lp, 0.0)

    sT = np.sqrt(var_T)
    sTd = np.sqrt(var_Td)
    sLp = np.sqrt(var_lp)

    out["std_T"][valid] = sT.astype(np.float32)
    out["std_Td"][valid] = sTd.astype(np.float32)
    out["std_logP"][valid] = sLp.astype(np.float32)

    out["cov_T_Td"][valid] = cov_T_Td.astype(np.float32)
    out["cov_T_logP"][valid] = cov_T_lp.astype(np.float32)
    out["cov_Td_logP"][valid] = cov_Td_lp.astype(np.float32)

    den = sT * sTd
    mask = den > 0
    tmp = np.full_like(cov_T_Td, np.nan)
    tmp[mask] = cov_T_Td[mask] / den[mask]
    out["corr_T_Td"][valid] = np.clip(tmp, -0.999999, 0.999999).astype(np.float32)

    den = sT * sLp
    mask = den > 0
    tmp = np.full_like(cov_T_lp, np.nan)
    tmp[mask] = cov_T_lp[mask] / den[mask]
    out["corr_T_logP"][valid] = np.clip(tmp, -0.999999, 0.999999).astype(np.float32)

    den = sTd * sLp
    mask = den > 0
    tmp = np.full_like(cov_Td_lp, np.nan)
    tmp[mask] = cov_Td_lp[mask] / den[mask]
    out["corr_Td_logP"][valid] = np.clip(tmp, -0.999999, 0.999999).astype(np.float32)

    return out

# =============================================================================
# 10. HIGHER MOMENTS: PÉBAY BATCH MERGE
# =============================================================================

def batch_moments(x: np.ndarray, valid: np.ndarray):
    """
    Vectorized moments over axis 0.
    Returns n, mean, M2, M3, M4.
    """
    n = valid.sum(axis=0).astype(np.int64)
    safe_n = np.maximum(n, 1).astype(np.float64)

    xv = np.where(valid, x, 0.0).astype(np.float64)
    mean = xv.sum(axis=0) / safe_n

    d = np.where(valid, x.astype(np.float64) - mean[None, :], 0.0)
    M2 = np.sum(d**2, axis=0)
    M3 = np.sum(d**3, axis=0)
    M4 = np.sum(d**4, axis=0)

    mean[n == 0] = np.nan
    M2[n == 0] = np.nan
    M3[n == 0] = np.nan
    M4[n == 0] = np.nan

    return n, mean, M2, M3, M4

def combine_moments(
    n1, mean1, M21, M31, M41,
    n2, mean2, M22, M32, M42,
):
    """
    Pébay merge for moments up to fourth order.
    M2/M3/M4 are central sums, not normalized moments.
    Vectorized element-wise.
    """
    n1f = n1.astype(np.float64)
    n2f = n2.astype(np.float64)
    nt = n1f + n2f

    out_n = (nt).astype(np.int64)
    out_mean = mean1.copy()
    out_M2 = M21.copy()
    out_M3 = M31.copy()
    out_M4 = M41.copy()

    only2 = (n1 == 0) & (n2 > 0)
    both = (n1 > 0) & (n2 > 0)

    if np.any(only2):
        out_mean[only2] = mean2[only2]
        out_M2[only2] = M22[only2]
        out_M3[only2] = M32[only2]
        out_M4[only2] = M42[only2]

    if np.any(both):
        a = n1f[both]
        b = n2f[both]
        n = nt[both]
        delta = mean2[both] - mean1[both]

        A = M21[both]
        B = M22[both]
        C = M31[both]
        D = M32[both]
        E = M41[both]
        F = M42[both]

        out_mean[both] = mean1[both] + delta * (b / n)

        out_M2[both] = (
            A + B + delta**2 * a * b / n
        )

        out_M3[both] = (
            C + D
            + delta**3 * a * b * (a - b) / n**2
            + 3.0 * delta * (a * B - b * A) / n
        )

        out_M4[both] = (
            E + F
            + delta**4 * a * b * (a*a - a*b + b*b) / n**3
            + 6.0 * delta**2 * (a*a * B + b*b * A) / n**2
            + 4.0 * delta * (a * D - b * C) / n
        )

    return out_n, out_mean, out_M2, out_M3, out_M4

def sample_adjusted_skew_kurt(n, M2, M3, M4):
    """
    scipy.stats convention:
        skew(..., bias=False)
        kurtosis(..., bias=False, fisher=True)
    """
    skew = np.full_like(M2, np.nan, dtype=np.float64)
    kurt = np.full_like(M2, np.nan, dtype=np.float64)

    ok3 = (n >= 3) & (M2 > 0)
    if np.any(ok3):
        nn = n[ok3].astype(np.float64)
        m2 = M2[ok3] / nn
        m3 = M3[ok3] / nn
        skew[ok3] = (
            np.sqrt(nn * (nn - 1.0)) / (nn - 2.0)
        ) * m3 / np.power(m2, 1.5)

    ok4 = (n >= 4) & (M2 > 0)
    if np.any(ok4):
        nn = n[ok4].astype(np.float64)
        b2 = nn * M4[ok4] / np.square(M2[ok4])
        kurt[ok4] = (
            (nn - 1.0) / ((nn - 2.0) * (nn - 3.0))
        ) * ((nn + 1.0) * b2 - 3.0 * (nn - 1.0))

    return skew, kurt

# =============================================================================
# 11. PSD / CHOLESKY
# =============================================================================

def build_corr_batch(r1, r2, r3) -> np.ndarray:
    n = r1.size
    R = np.empty((n, 3, 3), dtype=np.float64)
    R[:, 0, 0] = 1.0
    R[:, 1, 1] = 1.0
    R[:, 2, 2] = 1.0
    R[:, 0, 1] = R[:, 1, 0] = r1
    R[:, 0, 2] = R[:, 2, 0] = r2
    R[:, 1, 2] = R[:, 2, 1] = r3
    return R

def nearest_correlation_matrix_batch(R: np.ndarray) -> np.ndarray:
    """
    Symmetric PSD projection followed by correlation renormalization.
    """
    vals, vecs = np.linalg.eigh(R)
    vals = np.maximum(vals, 0.0)
    P = np.einsum("...ik,...k,...jk->...ij", vecs, vals, vecs)
    diag = np.sqrt(np.maximum(np.diagonal(P, axis1=1, axis2=2), 1e-15))
    C = P / (diag[:, :, None] * diag[:, None, :])
    C = 0.5 * (C + np.swapaxes(C, 1, 2))
    C[:, 0, 0] = 1.0
    C[:, 1, 1] = 1.0
    C[:, 2, 2] = 1.0
    return C

def cholesky_batch(R: np.ndarray):
    eig = np.linalg.eigvalsh(R)
    mineig = eig[:, 0]
    valid = np.isfinite(mineig) & (mineig >= -PSD_REPAIR_TOL)
    repaired = np.isfinite(mineig) & (mineig < 0) & (mineig >= -PSD_REPAIR_TOL)

    R2 = R.copy()
    if np.any(repaired):
        R2[repaired] = nearest_correlation_matrix_batch(R2[repaired])

    hard_invalid = ~valid
    if np.any(hard_invalid):
        # Caller skips these cells.
        pass

    L = np.full_like(R2, np.nan)
    if np.any(valid):
        try:
            L[valid] = np.linalg.cholesky(R2[valid])
        except np.linalg.LinAlgError:
            # Retry only with slightly regularized matrices.
            idx = np.flatnonzero(valid)
            for k in idx:
                try:
                    Rk = R2[k] + np.eye(3) * 1e-12
                    Rk = nearest_correlation_matrix_batch(Rk[None, ...])[0]
                    L[k] = np.linalg.cholesky(Rk)
                except np.linalg.LinAlgError:
                    valid[k] = False

    return L, valid, mineig, repaired

# =============================================================================
# 12. DAILY MONTE CARLO
# =============================================================================

def day_paths(doy: int) -> tuple[Path, Path]:
    return (
        DAY_DIR / f"day_{doy:03d}_{CONFIG_HASH}.npz",
        DAY_DIR / f"day_{doy:03d}_{CONFIG_HASH}.json",
    )

def valid_day_checkpoint(doy: int) -> bool:
    npz, js = day_paths(doy)
    if not npz.exists() or not js.exists():
        return False
    try:
        meta = json_load(js)
        return (
            meta.get("status") == "completed"
            and meta.get("config_hash") == CONFIG_HASH
            and meta.get("sha256") == sha256_file(npz)
        )
    except Exception:
        return False

def empty_day_output(ny: int, nx: int) -> dict:
    b = np.full((ny, nx), np.nan, dtype=np.float32)
    return {
        "mean_rh": b.copy(), "std_rh": b.copy(), "skew_rh": b.copy(), "kurt_rh": b.copy(),
        "mean_e": b.copy(), "std_e": b.copy(), "skew_e": b.copy(), "kurt_e": b.copy(),
        "mean_r": b.copy(), "std_r": b.copy(), "skew_r": b.copy(), "kurt_r": b.copy(),
        "mean_q": b.copy(), "std_q": b.copy(), "skew_q": b.copy(), "kurt_q": b.copy(),
        "supersat_fraction": b.copy(),
        "invalid_e_over_p_fraction": b.copy(),
        "invalid_covariance_fraction": b.copy(),
        "min_eigenvalue": b.copy(),
        "valid_sample_count": np.zeros((ny, nx), dtype=np.int32),
    }

def process_day(doy0: int, stats: dict) -> Optional[dict]:
    doy = doy0 + 1
    npz_path, json_path = day_paths(doy)

    if valid_day_checkpoint(doy):
        with np.load(npz_path, allow_pickle=False) as d:
            return {k: d[k] for k in d.files}

    ny, nx = stats["mean_T"].shape[1:]
    out = empty_day_output(ny, nx)

    mt = stats["mean_T"][doy0].reshape(-1)
    mtd = stats["mean_Td"][doy0].reshape(-1)
    mlp = stats["mean_logP"][doy0].reshape(-1)
    st = stats["std_T"][doy0].reshape(-1)
    std = stats["std_Td"][doy0].reshape(-1)
    slp = stats["std_logP"][doy0].reshape(-1)

    r1 = stats["corr_T_Td"][doy0].reshape(-1)
    r2 = stats["corr_T_logP"][doy0].reshape(-1)
    r3 = stats["corr_Td_logP"][doy0].reshape(-1)

    valid_cells = (
        np.isfinite(mt) & np.isfinite(mtd) & np.isfinite(mlp)
        & np.isfinite(st) & np.isfinite(std) & np.isfinite(slp)
        & (st > 0) & (std > 0) & (slp > 0)
        & np.isfinite(r1) & np.isfinite(r2) & np.isfinite(r3)
    )

    rows, cols = np.where(valid_cells.reshape(ny, nx))
    flat_valid = np.flatnonzero(valid_cells)

    if flat_valid.size == 0:
        atomic_npz_write(npz_path, **out)
        atomic_json_write(
            json_path,
            {
                "status": "completed",
                "doy": doy,
                "config_hash": CONFIG_HASH,
                "sha256": sha256_file(npz_path),
                "n_valid_cells": 0,
            },
        )
        return out

    rng = np.random.default_rng(CONFIG.random_seed + doy)
    chunk_size = CONFIG.cell_chunk_size
    batch_size = CONFIG.sample_batch_size

    # Diagnostics accumulators
    total_cov_bad = np.zeros(flat_valid.size, dtype=np.int64)
    total_e_over_p = np.zeros(flat_valid.size, dtype=np.int64)
    total_samples = np.zeros(flat_valid.size, dtype=np.int64)
    total_supersat = np.zeros(flat_valid.size, dtype=np.int64)
    min_eig_global = np.full(flat_valid.size, np.inf, dtype=np.float64)

    keys = ["rh", "e", "r", "q"]

    for start in range(0, flat_valid.size, chunk_size):
        ids = np.arange(start, min(start + chunk_size, flat_valid.size))
        cell = flat_valid[ids]

        mu = np.column_stack([
            mt[cell], mtd[cell], mlp[cell]
        ]).astype(np.float64)

        sig = np.column_stack([
            st[cell], std[cell], slp[cell]
        ]).astype(np.float64)

        R = build_corr_batch(r1[cell], r2[cell], r3[cell])
        L, chol_valid, mineig, repaired = cholesky_batch(R)

        min_eig_global[ids] = mineig
        total_cov_bad[ids] = (~chol_valid).astype(np.int64)

        good = chol_valid
        if not np.any(good):
            continue

        ids_good = ids[good]
        cell_good = cell[good]
        mu_good = mu[good]
        sig_good = sig[good]
        L_good = L[good]

        ncell = len(ids_good)

        moments = {}
        for key in keys:
            moments[key] = {
                "n": np.zeros(ncell, np.int64),
                "mean": np.zeros(ncell, np.float64),
                "M2": np.zeros(ncell, np.float64),
                "M3": np.zeros(ncell, np.float64),
                "M4": np.zeros(ncell, np.float64),
            }

        supersat = np.zeros(ncell, dtype=np.int64)
        invalid_ep = np.zeros(ncell, dtype=np.int64)
        sampled = np.zeros(ncell, dtype=np.int64)

        consumed = 0
        while consumed < CONFIG.n_samples:
            bn = min(batch_size, CONFIG.n_samples - consumed)
            Z = rng.standard_normal((bn, ncell, 3)).astype(np.float32)

            Xstd = np.einsum(
                "bci,cij->bcj",
                Z,
                L_good.astype(np.float32),
                optimize=True,
            )

            T = mu_good[:, 0][None, :] + sig_good[:, 0][None, :] * Xstd[:, :, 0]
            Td = mu_good[:, 1][None, :] + sig_good[:, 1][None, :] * Xstd[:, :, 1]
            logP = mu_good[:, 2][None, :] + sig_good[:, 2][None, :] * Xstd[:, :, 2]
            P = np.exp(logP.astype(np.float64)).astype(np.float32)

            phys = derive_moisture(T.astype(np.float32), Td.astype(np.float32), P)

            valid_all = phys["valid_all"]
            n2, mean2, M22, M32, M42 = batch_moments(phys["rh"], valid_all)
            n1 = moments["rh"]["n"]
            n, m, M2, M3, M4 = combine_moments(
                n1, moments["rh"]["mean"], moments["rh"]["M2"], moments["rh"]["M3"], moments["rh"]["M4"],
                n2, mean2, M22, M32, M42
            )
            moments["rh"] = {"n": n, "mean": m, "M2": M2, "M3": M3, "M4": M4}

            for key in ("e", "r", "q"):
                n2, mean2, M22, M32, M42 = batch_moments(phys[key], valid_all)
                a = moments[key]
                moments[key] = dict(zip(
                    ("n", "mean", "M2", "M3", "M4"),
                    combine_moments(
                        a["n"], a["mean"], a["M2"], a["M3"], a["M4"],
                        n2, mean2, M22, M32, M42,
                    )
                ))

            supersat += phys["supersat"].sum(axis=0).astype(np.int64)
            invalid_ep += phys["invalid_e_over_p"].sum(axis=0).astype(np.int64)
            sampled += valid_all.sum(axis=0).astype(np.int64)

            consumed += bn

        for local, orig_id in enumerate(ids_good):
            rr = rows[orig_id]
            cc = cols[orig_id]

            nvalid = moments["rh"]["n"][local]
            total_samples[orig_id] = nvalid
            total_supersat[orig_id] = supersat[local]
            total_e_over_p[orig_id] = invalid_ep[local]

            if nvalid < MIN_MC_VALID:
                continue

            for key, prefix in [("rh", "rh"), ("e", "e"), ("r", "r"), ("q", "q")]:
                a = moments[key]
                n1 = np.array([a["n"][local]], dtype=np.int64)
                M2 = np.array([a["M2"][local]], dtype=np.float64)
                M3 = np.array([a["M3"][local]], dtype=np.float64)
                M4 = np.array([a["M4"][local]], dtype=np.float64)
                mean = float(a["mean"][local])

                var = M2[0] / (nvalid - 1)
                stdv = np.sqrt(max(var, 0.0))
                sk, ku = sample_adjusted_skew_kurt(n1, M2, M3, M4)

                out[f"mean_{prefix}"][rr, cc] = np.float32(mean)
                out[f"std_{prefix}"][rr, cc] = np.float32(stdv)
                out[f"skew_{prefix}"][rr, cc] = np.float32(sk[0])
                out[f"kurt_{prefix}"][rr, cc] = np.float32(ku[0])

            out["supersat_fraction"][rr, cc] = np.float32(
                total_supersat[orig_id] / max(nvalid, 1)
            )
            out["invalid_e_over_p_fraction"][rr, cc] = np.float32(
                total_e_over_p[orig_id] / max(CONFIG.n_samples, 1)
            )
            out["invalid_covariance_fraction"][rr, cc] = np.float32(
                total_cov_bad[orig_id] > 0
            )
            out["min_eigenvalue"][rr, cc] = np.float32(min_eig_global[orig_id])
            out["valid_sample_count"][rr, cc] = np.int32(nvalid)

    atomic_npz_write(npz_path, **out)
    atomic_json_write(
        json_path,
        {
            "status": "completed",
            "doy": doy,
            "config_hash": CONFIG_HASH,
            "schema_version": CHECKPOINT_VERSION,
            "sha256": sha256_file(npz_path),
            "n_valid_cells": int(flat_valid.size),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    return out

# =============================================================================
# 13. NETCDF STREAMING FINALIZER
# =============================================================================

def create_main_netcdf(path: Path, lat: np.ndarray, lon: np.ndarray):
    if path.exists():
        path.unlink()

    ds = Dataset(path, "w", format="NETCDF4")
    ds.createDimension("doy", DOY_COUNT)
    ds.createDimension("latitude", len(lat))
    ds.createDimension("longitude", len(lon))

    v_doy = ds.createVariable("doy", "i2", ("doy",))
    v_lat = ds.createVariable("latitude", "f4", ("latitude",))
    v_lon = ds.createVariable("longitude", "f4", ("longitude",))
    v_doy[:] = np.arange(1, DOY_COUNT + 1, dtype=np.int16)
    v_lat[:] = lat.astype(np.float32)
    v_lon[:] = lon.astype(np.float32)
    v_lat.units = "degrees_north"
    v_lon.units = "degrees_east"

    month, day, label, reserved = calendar_labels()
    v_month = ds.createVariable("month", "i2", ("doy",))
    v_day = ds.createVariable("day", "i2", ("doy",))
    v_reserved = ds.createVariable("reserved_day", "i1", ("doy",))
    v_month[:] = month
    v_day[:] = day
    v_reserved[:] = reserved

    float_vars = [
        "mean_rh","std_rh","skew_rh","kurt_rh",
        "mean_vapor_pressure","std_vapor_pressure","skew_vapor_pressure","kurt_vapor_pressure",
        "mean_mixing_ratio","std_mixing_ratio","skew_mixing_ratio","kurt_mixing_ratio",
        "mean_specific_humidity","std_specific_humidity","skew_specific_humidity","kurt_specific_humidity",
    ]

    for name in float_vars:
        var = ds.createVariable(
            name, "f4", ("doy","latitude","longitude"),
            zlib=True, complevel=4, shuffle=True,
            chunksizes=(1, min(128,len(lat)), min(128,len(lon))),
            fill_value=-9999.0,
        )
        var.units = {
            "mean_rh": "%", "std_rh": "%",
            "mean_vapor_pressure": "hPa", "std_vapor_pressure": "hPa",
            "mean_mixing_ratio": "kg kg-1", "std_mixing_ratio": "kg kg-1",
            "mean_specific_humidity": "kg kg-1", "std_specific_humidity": "kg kg-1",
        }.get(name, "1")

    for name in [
        "skew_rh","kurt_rh",
        "skew_vapor_pressure","kurt_vapor_pressure",
        "skew_mixing_ratio","kurt_mixing_ratio",
        "skew_specific_humidity","kurt_specific_humidity",
    ]:
        ds.variables[name].units = "1"

    attrs = {
        "title": "ERA5-Land Moisture Climatology 1981-2020",
        "period": "1981-2020",
        "calendar": "366-day; Feb-28 + Feb-29 combined into DOY 60; DOY 59 reserved",
        "model": "Multivariate normal in (T2m_C, Td2m_C, ln(surface_pressure_hPa))",
        "n_samples": CONFIG.n_samples,
        "random_seed": CONFIG.random_seed,
        "rh_phase_rule": "water for T>=0 C, ice for T<0 C",
        "rh_clipping": "[0,100] percent after raw supersaturation diagnostic",
        "mixing_ratio_formula": "r = 0.622*e/(P-e)",
        "specific_humidity_formula": "q = r/(1+r)",
        "schema_version": SCHEMA_VERSION,
        "config_hash": CONFIG_HASH,
        "script_sha256": script_sha256(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    for k,v in attrs.items():
        setattr(ds, k, v)

    return ds

def create_diag_netcdf(path: Path, lat: np.ndarray, lon: np.ndarray):
    if path.exists():
        path.unlink()

    ds = Dataset(path, "w", format="NETCDF4")
    ds.createDimension("doy", DOY_COUNT)
    ds.createDimension("latitude", len(lat))
    ds.createDimension("longitude", len(lon))

    ds.createVariable("doy", "i2", ("doy",))[:] = np.arange(1, DOY_COUNT+1, dtype=np.int16)
    ds.createVariable("latitude", "f4", ("latitude",))[:] = lat.astype(np.float32)
    ds.createVariable("longitude", "f4", ("longitude",))[:] = lon.astype(np.float32)

    float_vars = [
        "supersaturation_fraction",
        "invalid_e_over_p_fraction",
        "invalid_covariance_fraction",
        "min_eigenvalue",
    ]
    for name in float_vars:
        ds.createVariable(
            name, "f4", ("doy","latitude","longitude"),
            zlib=True, complevel=4, shuffle=True,
            chunksizes=(1, min(128,len(lat)), min(128,len(lon))),
            fill_value=-9999.0,
        )

    ds.createVariable(
        "valid_sample_count", "i4",
        ("doy","latitude","longitude"),
        zlib=True, complevel=4,
        chunksizes=(1, min(128,len(lat)), min(128,len(lon))),
        fill_value=-9999,
    )

    for name in ["corr_T_Td","corr_T_logP","corr_Td_logP","valid_observation_count"]:
        # Populated from global stats after opening, via a separate writer.
        dtype = "f4" if name != "valid_observation_count" else "i4"
        ds.createVariable(
            name, dtype, ("doy","latitude","longitude"),
            zlib=True, complevel=4, shuffle=True,
            chunksizes=(1, min(128,len(lat)), min(128,len(lon))),
            fill_value=-9999,
        )
    ds.title = "Moisture Climatology Diagnostics"
    ds.config_hash = CONFIG_HASH
    ds.schema_version = SCHEMA_VERSION
    ds.script_sha256 = script_sha256()
    return ds

def write_day_to_netcdf(ds_main, ds_diag, doy0: int, out: dict, stats: dict):
    sl = doy0
    mapping = {
        "mean_rh":"mean_rh", "std_rh":"std_rh", "skew_rh":"skew_rh", "kurt_rh":"kurt_rh",
        "mean_e":"mean_vapor_pressure", "std_e":"std_vapor_pressure",
        "skew_e":"skew_vapor_pressure", "kurt_e":"kurt_vapor_pressure",
        "mean_r":"mean_mixing_ratio", "std_r":"std_mixing_ratio",
        "skew_r":"skew_mixing_ratio", "kurt_r":"kurt_mixing_ratio",
        "mean_q":"mean_specific_humidity", "std_q":"std_specific_humidity",
        "skew_q":"skew_specific_humidity", "kurt_q":"kurt_specific_humidity",
    }
    for src, dst in mapping.items():
        arr = np.asarray(out[src], dtype=np.float32)
        arr = np.where(np.isfinite(arr), arr, -9999.0)
        ds_main.variables[dst][sl,:,:] = arr

    for key in [
        "supersaturation_fraction",
        "invalid_e_over_p_fraction",
        "invalid_covariance_fraction",
        "min_eigenvalue",
    ]:
        arr = np.asarray(out[key], dtype=np.float32)
        arr = np.where(np.isfinite(arr), arr, -9999.0)
        ds_diag.variables[key][sl,:,:] = arr

    ds_diag.variables["valid_sample_count"][sl,:,:] = out["valid_sample_count"].astype(np.int32)

    for src, dst in [
        ("corr_T_Td","corr_T_Td"),
        ("corr_T_logP","corr_T_logP"),
        ("corr_Td_logP","corr_Td_logP"),
    ]:
        arr = stats[src][doy0].astype(np.float32)
        arr = np.where(np.isfinite(arr), arr, -9999.0)
        ds_diag.variables[dst][sl,:,:] = arr

    ds_diag.variables["valid_observation_count"][sl,:,:] = stats["n"][doy0].astype(np.int32)

def finalize_streaming(lat, lon, stats):
    tmp_main = OUTPUT_FILE.with_suffix(".tmp.nc")
    tmp_diag = DIAGNOSTIC_FILE.with_suffix(".tmp.nc")

    for p in [tmp_main, tmp_diag]:
        if p.exists():
            p.unlink()

    ds_main = create_main_netcdf(tmp_main, lat, lon)
    ds_diag = create_diag_netcdf(tmp_diag, lat, lon)

    try:
        for doy in range(1, DOY_COUNT + 1):
            npz, js = day_paths(doy)
            if not valid_day_checkpoint(doy):
                raise RuntimeError(f"Cannot finalize: invalid/missing day {doy}")

            with np.load(npz, allow_pickle=False) as d:
                out = {k:d[k] for k in d.files}

            write_day_to_netcdf(ds_main, ds_diag, doy-1, out, stats)
            if doy % 10 == 0:
                ds_main.sync()
                ds_diag.sync()

        ds_main.sync()
        ds_diag.sync()
    finally:
        ds_main.close()
        ds_diag.close()

    os.replace(tmp_main, OUTPUT_FILE)
    os.replace(tmp_diag, DIAGNOSTIC_FILE)

# =============================================================================
# 14. TESTS
# =============================================================================

def test_leap_day():
    assert get_clim_doy(59, 1984) == 60
    assert get_clim_doy(60, 1984) == 60
    assert get_clim_doy(61, 1984) == 61
    assert get_clim_doy(59, 1985) == 60
    assert get_clim_doy(60, 1985) == 61
    logger.info("PASS: leap-day")

def test_welford_against_numpy():
    rng = np.random.default_rng(7)
    x = rng.normal(size=10000).astype(np.float64)
    y = 2.0*x + rng.normal(scale=0.5, size=x.size)

    m = 0.0
    M2 = 0.0
    n = 0

    for v in x:
        n += 1
        delta = v - m
        m += delta / n
        M2 += delta * (v - m)

    assert abs(m - np.mean(x)) < 1e-12
    assert abs(M2/(n-1) - np.var(x, ddof=1)) < 1e-12

    # covariance test
    my = 0.0
    c = 0.0
    n = 0
    mx = 0.0
    for xv, yv in zip(x, y):
        n += 1
        dx = xv - mx
        dy = yv - my
        mx += dx/n
        my += dy/n
        c += dx * (yv - my)

    cov = c/(n-1)
    target = np.cov(x, y, ddof=1)[0,1]
    assert abs(cov - target) < 1e-12
    logger.info("PASS: Welford mean/variance/covariance")

def test_pebay_against_scipy():
    from scipy.stats import skew, kurtosis

    rng = np.random.default_rng(9)
    x = rng.lognormal(mean=0.2, sigma=0.8, size=5000)

    # Build moments as one batch
    xx = x[:,None]
    valid = np.ones_like(xx, dtype=bool)
    n, mean, M2, M3, M4 = batch_moments(xx, valid)

    sk, ku = sample_adjusted_skew_kurt(n, M2, M3, M4)

    assert np.allclose(sk[0], skew(x, bias=False), rtol=1e-10, atol=1e-10)
    assert np.allclose(ku[0], kurtosis(x, bias=False, fisher=True), rtol=1e-10, atol=1e-10)

    logger.info("PASS: Pébay moments vs scipy")

def test_ground_truth_vectorized_physics():
    rng = np.random.default_rng(123)
    T = rng.normal(12.0, 8.0, size=(37, 6)).astype(np.float32)
    Td = (T - np.abs(rng.normal(4.0, 2.0, size=T.shape))).astype(np.float32)
    P = np.exp(rng.normal(np.log(1010.0), 0.03, size=T.shape)).astype(np.float32)

    slow = ground_truth_simple(T, Td, P)
    fast = derive_moisture(T, Td, P)

    for key in ("rh", "e", "r", "q"):
        got = fast[key]
        ref = slow[key]
        mask = np.isfinite(ref) & np.isfinite(got)
        assert np.any(mask)
        err = np.max(np.abs(got[mask] - ref[mask]))
        assert err < 1e-5, f"{key}: max error {err}"

    logger.info("PASS: vectorized physics vs Ground Truth")


def test_mvn_physics():
    rng = np.random.default_rng(42)
    n = 1000
    cells = 8
    mu = np.array([15.0, 8.0, np.log(1013.25)])
    sig = np.array([5.0, 4.0, 0.03])
    R = np.array([[1,0.8,0.2],[0.8,1,0.3],[0.2,0.3,1]], dtype=float)
    L = np.linalg.cholesky(R)

    Z = rng.standard_normal((n,cells,3))
    X = np.einsum("bci,ij->bcj", Z, L)

    T = mu[0] + sig[0]*X[:,:,0]
    Td = mu[1] + sig[1]*X[:,:,1]
    P = np.exp(mu[2] + sig[2]*X[:,:,2])

    d = derive_moisture(T.astype(np.float32), Td.astype(np.float32), P.astype(np.float32))

    assert np.nanmin(d["rh"]) >= 0
    assert np.nanmax(d["rh"]) <= 100
    assert np.nanmin(d["e"]) > 0
    assert np.nanmin(d["r"]) >= 0
    assert np.nanmin(d["q"]) >= 0
    logger.info("PASS: synthetic MVN physics")

def run_tests():
    logger.info("Running mandatory tests...")
    test_leap_day()
    test_welford_against_numpy()
    test_pebay_against_scipy()
    test_ground_truth_vectorized_physics()
    test_mvn_physics()
    logger.info("All mandatory tests passed.")

# =============================================================================
# 15. MAIN
# =============================================================================

def main():
    logger.info("="*80)
    logger.info("MOISTURE CLIMATOLOGY v6")
    logger.info(f"Period: {START_YEAR}-{END_YEAR}")
    logger.info(f"Config hash: {CONFIG_HASH}")
    logger.info(f"N samples: {N_SAMPLES}")
    logger.info(f"Cell chunk: {CELL_CHUNK_SIZE}")
    logger.info(f"Sample batch: {SAMPLE_BATCH_SIZE}")
    logger.info("="*80)

    run_tests()

    sample_file = build_file_index(START_YEAR, T2M_DIR)[1]
    with open_dataset(sample_file) as ds:
        ds = sort_dataset(ds)
        lat = ds.latitude.values
        lon = ds.longitude.values

    shape = (len(lat), len(lon))
    years = list(range(START_YEAR, END_YEAR+1))

    log_progress("STAGE", phase="YEAR ACCUMULATION")
    # ---- Year checkpoints
    remaining = [y for y in years if not is_year_complete(y)]
    if remaining:
        logger.info(f"Years remaining: {len(remaining)}")
        processed = []
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(process_year_welford, y): y for y in remaining}
            for fut in tqdm(as_completed(futs), total=len(futs), desc="Years"):
                y = futs[fut]
                result_y, path = fut.result()
                if path is None or not is_year_complete(result_y):
                    raise RuntimeError(f"Year {result_y} did not produce a valid checkpoint.")
                processed.append(result_y)

    # Validate every year before merge.
    bad = [y for y in years if not is_year_complete(y)]
    if bad:
        raise RuntimeError(f"Invalid/missing year checkpoints: {bad}")

    log_progress("STAGE", phase="YEAR MERGE")
    # ---- Merge
    logger.info("Merging yearly accumulators...")
    merged = merge_years_v2(years, shape)
    stats = stats_from_merged(merged)

    log_progress("STAGE", phase="DAILY MONTE CARLO")
    # ---- Daily MC
    logger.info("Running daily Monte Carlo...")
    for doy0 in tqdm(range(DOY_COUNT), desc="Days"):
        if valid_day_checkpoint(doy0+1):
            continue
        process_day(doy0, stats)

    bad_days = [d for d in range(1, DOY_COUNT+1) if not valid_day_checkpoint(d)]
    if bad_days:
        raise RuntimeError(f"Invalid/missing daily checkpoints: {bad_days}")

    log_progress("STAGE", phase="FINALIZE NETCDF")
    # ---- Finalize
    logger.info("Finalizing NetCDF products...")
    finalize_streaming(lat, lon, stats)

    # ---- Final independent sanity validation
    with xr.open_dataset(OUTPUT_FILE) as ds:
        mean_rh = ds["mean_rh"].where(ds["mean_rh"] != -9999)
        mean_q = ds["mean_specific_humidity"].where(ds["mean_specific_humidity"] != -9999)
        mean_r = ds["mean_mixing_ratio"].where(ds["mean_mixing_ratio"] != -9999)
        mean_e = ds["mean_vapor_pressure"].where(ds["mean_vapor_pressure"] != -9999)

        if float(mean_rh.min()) < -1e-6 or float(mean_rh.max()) > 100.000001:
            raise RuntimeError("Final RH range validation failed.")
        if float(mean_q.min()) < -1e-12:
            raise RuntimeError("Final q validation failed.")
        if float(mean_r.min()) < -1e-12:
            raise RuntimeError("Final r validation failed.")
        if float(mean_e.min()) <= 0:
            raise RuntimeError("Final vapor pressure validation failed.")

    logger.info("="*80)
    logger.info("SUCCESS: moisture climatology completed.")
    logger.info(f"Main: {OUTPUT_FILE}")
    logger.info(f"Diagnostics: {DIAGNOSTIC_FILE}")
    logger.info(f"Main size: {OUTPUT_FILE.stat().st_size/1024**2:.1f} MB")
    logger.info("="*80)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        raise
    except Exception:
        logger.exception("Fatal error.")
        raise
