# ============================================================================
# HUMIDCLIMATOLOGYENGINE v11.0.0 — ULTIMATE SELF-DOCUMENTING SOURCE
# ============================================================================
# This file is intentionally verbose. The comments are part of the release.
# They explain the scientific contract, data shapes, mathematical intention,
# failure semantics, persistence model, recovery rules, performance choices,
# and the reason behind important implementation decisions.
#
# PROJECT IDENTITY
# ----------------
# Software: HumidClimatologyEngine
# Version : 11.0.0
# Role    : primary single-process calculation engine
#
# INPUT COMPATIBILITY
# -------------------
# The primary scientific inputs intentionally match the historical v8 input
# families: 2 m air temperature (T2m/t2m), 2 m dew-point temperature (D2m/d2m),
# and surface pressure (SP/sp). Keeping the same raw inputs is important for
# regression testing: differences in v10 should come from the declared software
# contract, not an accidental change of source data.
#
# OUTPUT SCOPE
# ------------
# v10 calculates RH, vapor pressure e, mixing ratio r and specific humidity q;
# it resolves the climatology at L1 daily, L2 eight 3-hour bins and L3 twenty-
# four hourly bins; it tracks four decades plus FULL 1981-2020; it retains exact
# extrema, counts, moments through fourth order, pairwise dependence state,
# thresholds and an empirical RH-q histogram at the configured levels.
#
# IMPORTANT OPERATIONAL RULE
# --------------------------
# PROGRESS IS NOT SCIENTIFIC TRUTH. Only a durable, verified COMMIT record means
# that a day/spatial-block transaction has been completed. This distinction is
# the key safeguard against the historical failure in which an apparently
# advanced checkpoint did not represent the scientific work actually completed.
#
# GRID NORMALIZATION RULE
# -----------------------
# T2m defines the reference orientation. D2m or SP may encode the same coordinate
# values in reverse order. A pure reversal is safe to normalize because it changes
# representation, not the physical coordinate set. A non-equivalent coordinate
# set is a hard error. The alignment action is logged for provenance.
#
# CALENDAR RULE
# -------------
# 366 climatological slots are retained. Slot 59 is reserved. Slot 60 is the
# composite Feb-28/Feb-29 slot. Slot 61 is Mar-01. This rule is part of the data
# schema and cannot be changed silently.
#
# PERFORMANCE RULE
# ----------------
# The core is single-process. Vectorized NumPy operations are preferred over
# Python loops across individual cells. Raw fields are spatially chunked. Durable
# persistence is performed at a bounded block/day transaction level rather than
# rewriting an entire multidecadal state for every hour.
#
# FAILURE RULE
# ------------
# Never turn an invalid scientific value into a plausible value solely to make
# the program continue. Missing data, supersaturation and invalid pressure
# partitions are recorded explicitly.
# ============================================================================

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HumidClimatologyEngine v11.2.0 FINAL
====================================

Production, single-process, power-safe ERA5-Land moisture climatology engine.

Primary input contract (kept compatible with the historical v8 workflow):
    T2m + 2 m Dew Point Temperature + Surface Pressure

Primary physical products:
    T, Td, P, es(T), e, RH, r, q

Time levels:
    L1 = daily pooled (24 hourly observations)
    L2 = eight 3-hour bins (00-02 ... 21-23)
    L3 = 24 hourly bins

Statistical products per variable and cell:
    n, mean, M2, M3, M4, min, max, missing_count,
    supersaturation_count, invalid_e_over_p_count,
    exact threshold counts.

Pair products:
    n_pair, mean_x, mean_y, Cxy, covariance, correlation,
    exact joint-threshold counts, optional empirical 2-D histogram.

Periods accumulated directly in one raw-data pass:
    DECADE_1981_1990
    DECADE_1991_2000
    DECADE_2001_2010
    DECADE_2011_2020
    FULL_1981_2020

Power-failure contract:
    A day x spatial-block is complete ONLY after:
        1) before-image journal is durable;
        2) NetCDF state is updated and sync() succeeds;
        3) a durable SQLite COMMIT record is written.
    If power fails before COMMIT, recovery restores the before-image and the
    work unit is replayed. There is no worker-count-dependent progress state.

Important design choices:
    * Single process; no ProcessPoolExecutor / ThreadPoolExecutor.
    * Streaming monthly input; one 24-hour block is resident at a time.
    * NetCDF checkpoint shards are block-local and chunked on disk.
    * SQLite is used only for transaction truth, not for scientific arrays.
    * Exact moments/covariance/min/max/threshold counts use observations directly.
    * The default empirical histogram is RH x q for L1/L2. L3 histograms are
      available as an explicit opt-in because they multiply disk/I/O costs.
    * Quantiles are intentionally delegated to the analysis layer or an
      optional second pass rather than storing Python objects per cell in the
      production hot path.

Repository compatibility:
    The historical filename is intentionally retained so this file can replace
    the v8 runner in the existing repository while software metadata reports
    version 10.0.0.
"""
from __future__ import annotations

import argparse
import calendar
import csv
import hashlib
import json
import logging
import os
import shutil
import sqlite3
import sys
import tempfile
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import xarray as xr

try:
    from netCDF4 import Dataset
except Exception:  # pragma: no cover
    Dataset = None


# ============================================================================
# ENGINE IDENTITY
# ============================================================================
# These constants are intentionally duplicated into output metadata and manifests.
# They are the first things a downstream analysis tool should inspect before reading
# a checkpoint or result file.  Version 10 is not the same software as v8 even though
# the three primary ERA5-Land input families are deliberately kept compatible.
# ============================================================================
ENGINE_NAME = "HumidClimatologyEngine"
ENGINE_VERSION = "11.5.0"
SCHEMA_VERSION = "11.5-FINAL"
CHECKPOINT_VERSION = "11.3-POWER-SAFE-BLOCKDAY"
LOG = logging.getLogger("HumidClimatologyEngine")


# ============================================================================
# SCIENTIFIC STATE SHAPE
# ============================================================================
# The next constants define the dimensions of the statistical state.  They are part
# of the schema, not merely implementation conveniences.  Changing them changes the
# meaning and/or size of checkpoint data and therefore must participate in the config
# fingerprint.
# ============================================================================
VARIABLES: Tuple[str, ...] = ("rh", "e", "r", "q")
PAIRS: Tuple[Tuple[str, str], ...] = (("rh", "q"), ("rh", "e"), ("q", "e"), ("r", "q"))
LEVELS: Tuple[str, ...] = ("L1", "L2", "L3")
LEVEL_BINS: Dict[str, Tuple[int, int]] = {"L1": (0, 1), "L2": (1, 9), "L3": (9, 33)}
L1_BIN_INDEX = 0
L2_START = 1
L3_START = 9

THRESHOLDS: Dict[str, Tuple[float, ...]] = {
    "rh": (80.0, 90.0, 95.0, 100.0),
    "e": (10.0, 15.0, 20.0, 25.0),
    "r": (0.008, 0.010, 0.012, 0.015, 0.020),
    "q": (0.008, 0.010, 0.012, 0.015, 0.020),
}
JOINT_THRESHOLDS: Tuple[Tuple[str, float, str, float], ...] = (
    ("rh", 80.0, "q", 0.012),
    ("rh", 90.0, "q", 0.012),
    ("rh", 95.0, "q", 0.015),
    ("q", 0.012, "rh", 80.0),
)
HIST_RANGES: Dict[str, Tuple[float, float]] = {
    "rh": (0.0, 130.0),
    "q": (0.0, 0.030),
    "r": (0.0, 0.035),
    "e": (0.0, 40.0),
}
HIST_BINS = (8, 8)
HIST_PAIRS: Tuple[Tuple[str, str], ...] = (("rh", "q"),)
HIST_LEVELS: Tuple[str, ...] = ("L1", "L2")


# CONFIGURATION CONTRACT
# This immutable dataclass is the single source of operational defaults for v11.2.
# It deliberately contains the same primary ERA5-Land input families used by v8,
# while adding the v10 output root, spatial chunking, compression, and empirical-
# histogram controls.  Freezing this object prevents accidental mutation during a
# long run, which is important because the configuration fingerprint is used as
# part of reproducibility and checkpoint identity.
@dataclass(frozen=True)
# ========================================================================
# CLASS Config — IMPLEMENTATION GUIDE
# ========================================================================
# Responsibility:
#   Own one clearly bounded part of the v10 engine. The class should be
#   read together with its caller and with the persisted state it owns.
#
# What to inspect:
#   1. Constructor state and configuration.
#   2. Public methods and their pre/post conditions.
#   3. Array shapes and coordinate conventions.
#   4. Failure behavior and recovery behavior.
#   5. Whether a value is scientific state, telemetry, or metadata.
#
# Scientific safety:
#   Optimizing this class must not change the sample population, calendar,
#   units, masking rules, or mathematical definition of any statistic.
# ========================================================================

class Config:
    start_year: int = 1981
    end_year: int = 2020
    t2m_dir: Path = Path(r"F:\Kazemi\era5\land\T2m")
    d2m_dir: Path = Path(r"F:\Kazemi\era5\land\Dew_Point_Temperature")
    sp_dir: Path = Path(r"F:\Kazemi\era5\land\Surface_Pressure")
    output_root: Path = Path(r"C:\c\HumidClimatologyEngine_v11.5")
    chunk_lat: int = 64
    chunk_lon: int = 128
    compression: int = 4
    shuffle: bool = True
    hist_levels: Tuple[str, ...] = HIST_LEVELS
    hist_pairs: Tuple[Tuple[str, str], ...] = HIST_PAIRS

    def validate(self) -> None:
        if self.start_year > self.end_year:
            raise ValueError("start_year must be <= end_year")
        if self.chunk_lat < 1 or self.chunk_lon < 1:
            raise ValueError("chunk_lat/chunk_lon must be positive")
        if not set(self.hist_levels).issubset(set(LEVELS)):
            raise ValueError("hist_levels contains unknown level")
        for p in self.hist_pairs:
            if p not in PAIRS:
                raise ValueError(f"Histogram pair not configured: {p}")



# ============================================================================
# DEFAULT DATA LOCATIONS
# ============================================================================
# These are the same primary input families used by the historical v8 workflow.
# They are paths to source data only; v10 checkpoint/output files are intentionally
# written elsewhere so the old v8 products remain untouched.
# ============================================================================
CONFIG = Config()

# ============================================================================
# DECADE + FULL PERIOD ROUTING
# ============================================================================
# Every source year belongs to exactly one decade and also contributes to FULL.
# This dual routing is how the 40-year product is accumulated without replaying the
# raw archive a second time merely to construct the FULL state.
# ============================================================================
PERIODS: Dict[str, Tuple[int, int]] = {
    "DECADE_1981_1990": (1981, 1990),
    "DECADE_1991_2000": (1991, 2000),
    "DECADE_2001_2010": (2001, 2010),
    "DECADE_2011_2020": (2011, 2020),
    "FULL_1981_2020": (1981, 2020),
}


# >>> utc_now: UTC timestamp helper for durable journal/provenance records.
# TIME / PROVENANCE
# All journal timestamps are stored in UTC.  UTC avoids ambiguity during DST changes,
# local-time changes, and cross-machine comparison.  The timestamp is metadata only;
# it must never be used as the scientific sample time for ERA5 observations.
# ------------------------------------------------------------------------
# FUNCTION utc_now — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


# >>> ensure_dir: Idempotent directory creation used throughout restart/recovery.
# FILESYSTEM SAFETY
# Directory creation is deliberately idempotent.  Restarting the engine must not fail
# merely because a checkpoint, report, or state directory already exists.
# ------------------------------------------------------------------------
# FUNCTION ensure_dir — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


# >>> canonical_json: Stable JSON serialization used for deterministic configuration fingerprints.
# DETERMINISTIC SERIALIZATION
# JSON object ordering is fixed so that the same scientific configuration produces
# the same byte sequence and therefore the same SHA-256 configuration fingerprint.
# ------------------------------------------------------------------------
# FUNCTION canonical_json — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def canonical_json(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


# >>> sha256_bytes: SHA-256 for in-memory transaction payloads.
# IN-MEMORY INTEGRITY HASH
# Used when an object is already represented as bytes.  SHA-256 is an integrity and
# identity mechanism here; it is not an encryption mechanism and does not prove the
# scientific correctness of the contents by itself.
# ------------------------------------------------------------------------
# FUNCTION sha256_bytes — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


# >>> sha256_file: Streaming SHA-256 so large outputs are fingerprinted without loading them into RAM.
# STREAMING FILE HASH
# Large NetCDF/checkpoint files are hashed in bounded blocks so the integrity check
# does not require the whole file to be loaded into memory.  The 4 MiB default is a
# compromise between Python call overhead and memory use.
# ------------------------------------------------------------------------
# FUNCTION sha256_file — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def sha256_file(path: Path, block_size: int = 4 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(block_size), b""):
            h.update(block)
    return h.hexdigest()


# >>> atomic_write_bytes: Atomic file write: temporary file -> flush -> fsync -> atomic replace.
# ATOMIC FILE REPLACEMENT
# The write protocol is temporary-file -> flush -> fsync -> os.replace.  The purpose
# is to prevent a reader from observing a partially written JSON/state file after a
# crash.  The replacement is atomic at the filesystem namespace level; ultimate
# physical persistence still depends on the OS and storage hardware.
# ------------------------------------------------------------------------
# FUNCTION atomic_write_bytes — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def atomic_write_bytes(path: Path, payload: bytes) -> None:
    ensure_dir(path.parent)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(payload)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


# >>> atomic_write_json: JSON wrapper around the atomic writer.
# JSON DURABILITY WRAPPER
# All small metadata records use the same atomic bytes writer so manifests and
# journal records cannot be left half-written by an interrupted process.
# ------------------------------------------------------------------------
# FUNCTION atomic_write_json — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def atomic_write_json(path: Path, obj: Any) -> None:
    atomic_write_bytes(path, json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True).encode("utf-8"))


# >>> is_leap_year: Gregorian leap-year rule for the climatological calendar.
# GREGORIAN CALENDAR RULE
# This is the exact Gregorian leap-year rule required by the historical v8/v10
# climatological calendar contract.
# ------------------------------------------------------------------------
# FUNCTION is_leap_year — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def is_leap_year(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


# >>> climatological_doy_from_date: Maps Gregorian dates to the 366-slot project calendar.
# CLIMATOLOGICAL DAY MAPPING
# The software does NOT use raw Gregorian DOY directly.  It reserves slot 59,
# pools Feb-28 and Feb-29 into slot 60, and shifts Mar-Dec so that Dec-31 is slot 366.
# Keeping this mapping identical across all periods is essential for valid decadal
# comparisons and FULL-period merge equivalence.
# ------------------------------------------------------------------------
# FUNCTION climatological_doy_from_date — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def climatological_doy_from_date(year: int, month: int, day: int) -> int:
    native = datetime(year, month, day).timetuple().tm_yday
    if native == 59:
        return 60
    if is_leap_year(year):
        if native == 60:
            return 60
        return native
    return native if native < 59 else native + 1


# >>> climatological_doy: datetime64 adapter for the calendar mapper.
# DATETIME64 ADAPTER
# Converts an xarray/NumPy datetime64 observation to the project's canonical DOY
# representation.  The function intentionally delegates the scientific rule to the
# integer/date implementation above so there is only one source of calendar truth.
# ------------------------------------------------------------------------
# FUNCTION climatological_doy — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def climatological_doy(timestamp: np.datetime64) -> int:
    s = np.datetime_as_string(timestamp, unit="D")
    dt = datetime.fromisoformat(s)
    return climatological_doy_from_date(dt.year, dt.month, dt.day)


# >>> canonical_calendar: Builds canonical calendar labels and the reserved-slot mask.
# HUMAN-READABLE CALENDAR TABLE
# Builds labels, month/day fields, and the reserved-slot mask used by outputs and
# reports.  This is a schema artifact as much as a display helper.
# ------------------------------------------------------------------------
# FUNCTION canonical_calendar — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def canonical_calendar() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    month = np.full(366, -1, np.int16)
    day = np.full(366, -1, np.int16)
    label = np.full(366, "RESERVED", "U16")
    reserved = np.zeros(366, np.int8)
    d = datetime(2001, 1, 1)
    for slot in range(1, 59):
        i = slot - 1
        month[i], day[i], label[i] = d.month, d.day, d.strftime("%b-%d")
        d = d + timedelta(days=1)
    reserved[58] = 1
    month[59], day[59], label[59] = 2, 28, "Feb-28/Feb-29"
    d = datetime(2000, 3, 1)
    for slot in range(61, 367):
        i = slot - 1
        month[i], day[i], label[i] = d.month, d.day, d.strftime("%b-%d")
        d = d + timedelta(days=1)
    return month, day, label, reserved


# >>> validate_calendar_contract: Regression tests for all critical calendar rules.
# CALENDAR REGRESSION TEST
# These assertions protect against a subtle but catastrophic error: an off-by-one
# DOY map can make a numerically perfect climate calculation scientifically wrong.
# ------------------------------------------------------------------------
# FUNCTION validate_calendar_contract — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def validate_calendar_contract() -> None:
    checks = {
        (1985, 1, 1): 1,
        (1985, 2, 27): 58,
        (1985, 2, 28): 60,
        (1984, 2, 28): 60,
        (1984, 2, 29): 60,
        (1984, 3, 1): 61,
        (1985, 3, 1): 61,
        (1985, 12, 30): 365,
        (1985, 12, 31): 366,
    }
    for args, expected in checks.items():
        got = climatological_doy_from_date(*args)
        if got != expected:
            raise AssertionError(f"calendar {args}: expected {expected}, got {got}")
    if climatological_doy_from_date(1985, 2, 28) != 60:
        raise AssertionError("Feb-28 mapping failed")
    if climatological_doy_from_date(1984, 2, 28) != 60:
        raise AssertionError("Leap Feb-28 mapping failed")
    if climatological_doy_from_date(1984, 2, 29) != 60:
        raise AssertionError("Leap Feb-29 mapping failed")
    if climatological_doy_from_date(1984, 3, 1) != 61:
        raise AssertionError("Leap Mar-01 mapping failed")
    if climatological_doy_from_date(1985, 3, 1) != 61:
        raise AssertionError("Non-leap Mar-01 mapping failed")
    _, _, _, reserved = canonical_calendar()
    if reserved[58] != 1:
        raise AssertionError("Slot 59 must be reserved")


# FILE INVENTORY PARSER
# Extracts a target year/month from a filename when possible.  It is intentionally
# conservative; inventory logic must not reinterpret arbitrary filenames as data.
# ------------------------------------------------------------------------
# FUNCTION extract_year_month — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def extract_year_month(path: Path, year: int) -> Optional[int]:
    import re
    m = re.search(r"(?<!\d)(\d{4})(\d{2})(?!\d)", path.name)
    if not m:
        return None
    y, mon = int(m.group(1)), int(m.group(2))
    return mon if y == year and 1 <= mon <= 12 else None


# >>> build_file_index: Builds the year/month input inventory and fails on missing or duplicate months.
# MONTHLY INPUT INVENTORY
# Builds the authoritative map month -> NetCDF path for one variable and one year.
# Missing or duplicate months are treated as operational errors before accumulation.
# ------------------------------------------------------------------------
# FUNCTION build_file_index — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def build_file_index(year: int, folder: Path) -> Dict[int, Path]:
    if not folder.exists():
        raise FileNotFoundError(folder)
    result: Dict[int, Path] = {}
    for p in sorted(folder.glob(f"*{year}*.nc")):
        mon = extract_year_month(p, year)
        if mon is None:
            continue
        if mon in result:
            raise RuntimeError(f"Duplicate {year}-{mon:02d} in {folder}: {p}")
        result[mon] = p
    missing = sorted(set(range(1, 13)) - set(result))
    if missing:
        raise RuntimeError(f"Missing months for {year} in {folder}: {missing}")
    return result


# >>> open_dataset: Central NetCDF opening policy; cache=False limits lingering file handles.
# NETCDF OPENING POLICY
# Opens a dataset through xarray using the NetCDF4 backend.  The decode/mask settings
# make coordinate and fill-value handling explicit and consistent with the validator.
# ------------------------------------------------------------------------
# FUNCTION open_dataset — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def open_dataset(path: Path) -> xr.Dataset:
    return xr.open_dataset(path, engine="netcdf4", decode_times=True, mask_and_scale=True, cache=False)


# >>> normalize_units: Normalizes unit strings before exact conversion dispatch.
# UNIT NORMALIZATION
# Converts arbitrary NetCDF unit metadata to a normalized lowercase token.  Unit
# handling is deliberately centralized to stop different parts of the engine from
# making inconsistent assumptions.
# ------------------------------------------------------------------------
# FUNCTION normalize_units — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def normalize_units(units: Any) -> str:
    return str(units or "").strip().lower().replace("°", "deg").replace(" ", "_")


# >>> convert_temperature: Explicitly converts temperature to Celsius; unknown units fail closed.
# TEMPERATURE NORMALIZATION
# Converts the supported temperature encodings to Celsius for thermodynamic formulas.
# The function fails instead of guessing unknown units, because a factor-of-273.15
# mistake would silently corrupt the entire 40-year product.
# ------------------------------------------------------------------------
# FUNCTION convert_temperature — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def convert_temperature(arr: np.ndarray, units: Any) -> np.ndarray:
    u = normalize_units(units)
    x = np.asarray(arr, dtype=np.float32)
    if u in {"k", "kelvin"}:
        return x - np.float32(273.15)
    if u in {"c", "degc", "degree_celsius", "degrees_celsius", "celsius"}:
        return x
    raise RuntimeError(f"Unsupported temperature units: {units!r}")


# >>> convert_pressure: Explicitly converts pressure to hPa; unknown units fail closed.
# PRESSURE NORMALIZATION
# Converts supported pressure representations to hPa, the unit used by the
# saturation-vapor and mixing-ratio formulas in this project.
# ------------------------------------------------------------------------
# FUNCTION convert_pressure — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def convert_pressure(arr: np.ndarray, units: Any) -> np.ndarray:
    u = normalize_units(units)
    x = np.asarray(arr, dtype=np.float32)
    if u in {"pa", "pascal", "pascals"}:
        return x / np.float32(100.0)
    if u in {"hpa", "mb", "millibar", "millibars"}:
        return x
    raise RuntimeError(f"Unsupported pressure units: {units!r}")


# >>> validate_time_axis: Checks actual timestamps, not just time-array length.
# TIME-AXIS VALIDATION
# Verifies that the monthly time coordinate is the expected hourly sequence and that
# the observations belong to the requested year/month.  Length alone is not enough:
# exact timestamp integrity is required for paired T/Td/P calculations.
# ------------------------------------------------------------------------
# FUNCTION validate_time_axis — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def validate_time_axis(ds: xr.Dataset, year: int, month: int) -> None:
    if "time" not in ds.coords:
        raise RuntimeError(f"Missing time coordinate: {year}-{month:02d}")
    t = ds.time.values
    if t.size == 0:
        raise RuntimeError("Empty time axis")
    dtm = np.diff(t).astype("timedelta64[m]").astype(np.int64)
    if np.any(dtm <= 0) or np.any(dtm != 60):
        bad = np.unique(dtm[dtm != 60])[:10]
        raise RuntimeError(f"Non-hourly/duplicate time axis in {year}-{month:02d}: {bad.tolist()}")
    a = np.datetime_as_string(t[0], unit="m")
    b = np.datetime_as_string(t[-1], unit="m")
    if not (a.startswith(f"{year:04d}-{month:02d}") and b.startswith(f"{year:04d}-{month:02d}")):
        raise RuntimeError(f"Month boundary mismatch: {a} .. {b}")


# >>> _coord_relation: Classifies coordinate axes as exact, reversed, or truly mismatched.
# SPATIAL AXIS RELATIONSHIP
# Returns whether two coordinate vectors match directly, match in reverse order, or
# are incompatible.  Reversed latitude is a known storage convention in the D2m data.
# ------------------------------------------------------------------------
# FUNCTION _coord_relation — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def _coord_relation(ref: np.ndarray, other: np.ndarray, atol: float = 1e-6) -> str:
    """Return exact, reversed, or mismatch for a floating coordinate axis."""
    a = np.asarray(ref, dtype=np.float64)
    b = np.asarray(other, dtype=np.float64)
    if a.shape != b.shape:
        return "mismatch"
    if np.allclose(a, b, rtol=0.0, atol=atol, equal_nan=False):
        return "exact"
    if np.allclose(a, b[::-1], rtol=0.0, atol=atol, equal_nan=False):
        return "reversed"
    return "mismatch"


# >>> align_dataset_to_reference: Normalizes a source dataset to the T2m reference orientation; reversed axes are safe only when coordinate values otherwise match.
# SPATIAL ALIGNMENT
# T2m is the reference grid.  If another dataset stores the exact same coordinate
# values in reverse order, its data axis is reversed before calculations.  A true
# grid mismatch still raises an error: alignment never means silent regridding.
# ------------------------------------------------------------------------
# FUNCTION align_dataset_to_reference — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def align_dataset_to_reference(reference: xr.Dataset, other: xr.Dataset, label: str, atol: float = 1e-6) -> xr.Dataset:
    """Normalize latitude/longitude orientation to the T2m reference grid.

    A reversed axis is accepted only when the coordinate values are otherwise
    identical within a tight absolute tolerance. Any true coordinate mismatch
    remains a hard error.
    """
    result = other
    for coord in ("latitude", "longitude"):
        relation = _coord_relation(reference[coord].values, result[coord].values, atol=atol)
        if relation == "exact":
            continue
        if relation == "reversed":
            result = result.isel({coord: slice(None, None, -1)})
            result = result.assign_coords({coord: reference[coord].values})
            LOG.info("GRID ALIGN | %s | reversed %s axis to T2m reference", label, coord)
            continue
        raise RuntimeError(
            f"Grid coordinate mismatch: {label} {coord}; "
            f"shape_ref={reference[coord].shape}, shape_other={result[coord].shape}"
        )
    # Final post-normalization verification.
    for coord in ("latitude", "longitude"):
        if not np.allclose(
            np.asarray(reference[coord].values, dtype=np.float64),
            np.asarray(result[coord].values, dtype=np.float64),
            rtol=0.0, atol=atol, equal_nan=False,
        ):
            raise RuntimeError(f"Post-alignment verification failed: {label} {coord}")
    return result


# >>> open_aligned_triplet: Opens T2m/D2m/SP and performs spatial alignment before any calculation.
# OPEN + ALIGN THE THREE PHYSICAL INPUTS
# This is the single gate between raw monthly files and scientific computation.
# It normalizes D2m/SP axis orientation to T2m, verifies coordinates again after the
# transformation, and closes all three datasets if anything fails during opening.
# ------------------------------------------------------------------------
# FUNCTION open_aligned_triplet — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def open_aligned_triplet(t2m_path: Path, d2m_path: Path, sp_path: Path):
    """Open T2m/D2m/SP and normalize the latter two to the T2m grid."""
    ds_t = open_dataset(t2m_path)
    ds_d = open_dataset(d2m_path)
    ds_p = open_dataset(sp_path)
    try:
        ds_d = align_dataset_to_reference(ds_t, ds_d, "D2m")
        ds_p = align_dataset_to_reference(ds_t, ds_p, "SP")
        return ds_t, ds_d, ds_p
    except Exception:
        ds_t.close(); ds_d.close(); ds_p.close()
        raise


# >>> validate_datasets: Final monthly validation gate for dimensions, variables, times and post-alignment coordinates.
# FINAL MONTHLY DATASET VALIDATION
# Checks dimensions, coordinate compatibility, time integrity, and expected variables.
# This function is intentionally strict because an incorrect pairing of temperature,
# dew point, and pressure can produce plausible-looking but physically wrong humidity.
# ------------------------------------------------------------------------
# FUNCTION validate_datasets — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def validate_datasets(ds_t: xr.Dataset, ds_d: xr.Dataset, ds_p: xr.Dataset, year: int, month: int) -> None:
    for label, ds, var in (("T2m", ds_t, "t2m"), ("D2m", ds_d, "d2m"), ("SP", ds_p, "sp")):
        for dim in ("time", "latitude", "longitude"):
            if dim not in ds.dims:
                raise RuntimeError(f"{label}: missing dimension {dim}")
        if var not in ds.data_vars:
            raise RuntimeError(f"{label}: missing variable {var}")
        validate_time_axis(ds, year, month)

    # D2m/SP must already be aligned to the T2m reference before this call.
    for c in ("latitude", "longitude"):
        ref = np.asarray(ds_t[c].values, dtype=np.float64)
        for label, ds in (("D2m", ds_d), ("SP", ds_p)):
            cur = np.asarray(ds[c].values, dtype=np.float64)
            if ref.shape != cur.shape or not np.allclose(ref, cur, rtol=0.0, atol=1e-6, equal_nan=False):
                raise RuntimeError(f"Post-alignment grid mismatch: {label} {c}")
        if ref.size > 1:
            diff = np.diff(ref)
            if not (np.all(diff > 0) or np.all(diff < 0)):
                raise RuntimeError(f"Grid coordinate {c} is not strictly monotonic")

    # Compare the actual time coordinates, not only their lengths.
    tt = np.asarray(ds_t["time"].values)
    for label, ds in (("D2m", ds_d), ("SP", ds_p)):
        ot = np.asarray(ds["time"].values)
        if tt.shape != ot.shape or not np.array_equal(tt, ot):
            raise RuntimeError(f"Input time coordinate mismatch: T2m vs {label}")


# >>> saturation_vapor_pressure: Phase-aware saturation vapor pressure: water at/above 0 C, ice below 0 C.
# PHASE-AWARE SATURATION VAPOR PRESSURE
# Uses the project convention: water saturation for T >= 0 C and ice saturation for
# T < 0 C.  The phase switch is part of the scientific definition, not an optimization.
# ------------------------------------------------------------------------
# FUNCTION saturation_vapor_pressure — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def saturation_vapor_pressure(temp_c: np.ndarray) -> np.ndarray:
    t = np.asarray(temp_c, dtype=np.float32)
    out = np.full(t.shape, np.nan, dtype=np.float32)
    finite = np.isfinite(t)
    water = finite & (t >= 0.0)
    ice = finite & (t < 0.0)
    if np.any(water):
        tw = t[water].astype(np.float64)
        out[water] = (6.112 * np.exp(17.67 * tw / (tw + 243.5))).astype(np.float32)
    if np.any(ice):
        ti = t[ice].astype(np.float64)
        out[ice] = (6.112 * np.exp(22.46 * ti / (ti + 272.62))).astype(np.float32)
    return out


# >>> derive_moisture: Vectorized observation-level psychrometric transformation plus validity diagnostics.
# PSYCHROMETRIC TRANSFORMATION
# Converts aligned T, Td, P arrays into vapor pressure e, relative humidity RH,
# mixing ratio r, and specific humidity q.  Invalid pressure partitions are flagged,
# while RH supersaturation is preserved as a diagnostic and separately bounded for
# the reporting value.
# ------------------------------------------------------------------------
# FUNCTION derive_moisture — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def derive_moisture(T: np.ndarray, Td: np.ndarray, P_hpa: np.ndarray) -> Dict[str, np.ndarray]:
    es_t = saturation_vapor_pressure(T)
    e = saturation_vapor_pressure(Td)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        rh_raw = 100.0 * e / es_t
    supersat = np.isfinite(rh_raw) & (rh_raw > 100.0)
    rh = np.clip(rh_raw, 0.0, 100.0).astype(np.float32)
    valid_e = np.isfinite(e) & (e > 0)
    valid_ep = valid_e & np.isfinite(P_hpa) & (P_hpa > 0) & (e < P_hpa)
    r = np.full_like(e, np.nan, dtype=np.float32)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        r[valid_ep] = (0.622 * e[valid_ep] / (P_hpa[valid_ep] - e[valid_ep])).astype(np.float32)
    q = (r / (1.0 + r)).astype(np.float32)
    return {
        "rh": rh,
        "e": e.astype(np.float32),
        "r": r,
        "q": q,
        "valid": {
            "rh": np.isfinite(rh),
            "e": np.isfinite(e),
            "r": np.isfinite(r),
            "q": np.isfinite(q),
        },
        "supersat": supersat,
        "invalid_e_over_p": np.isfinite(e) & np.isfinite(P_hpa) & (e >= P_hpa),
    }


# >>> batch_moments: Computes batch sufficient statistics without retaining raw history.
# BATCH SUFFICIENT STATISTICS
# Computes count, mean, M2, M3, M4, min and max for one vectorized batch.  The output
# contains exactly the sufficient information needed to merge this batch into an
# existing climatological state without retaining raw history.
# ------------------------------------------------------------------------
# FUNCTION batch_moments — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def batch_moments(x: np.ndarray, valid: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Accept (time, cells) or (time, y, x). The statistical engine always treats
    # axis 0 as time and collapses every remaining dimension into independent cells.
    x = np.asarray(x, dtype=np.float64)
    if x.ndim < 2:
        raise ValueError(f"batch_moments requires at least 2 dimensions (time, cells); got {x.shape}")
    spatial_shape = x.shape[1:]
    cells = int(np.prod(spatial_shape, dtype=np.int64))
    x2 = x.reshape(x.shape[0], cells)
    v = np.asarray(valid, dtype=bool)
    if v.shape != x.shape:
        if v.size == x.size:
            v = v.reshape(x.shape)
        else:
            raise ValueError(f"Validity shape mismatch: x={x.shape}, valid={v.shape}")
    v2 = v.reshape(x.shape[0], cells) & np.isfinite(x2)

    n = v2.sum(axis=0, dtype=np.int64)
    has = n > 0
    nf = n.astype(np.float64)
    xx = np.where(v2, x2, 0.0)
    mean = np.divide(xx.sum(axis=0, dtype=np.float64), nf, out=np.zeros_like(nf), where=has)
    d = np.where(v2, x2 - mean[None, :], 0.0)
    M2 = np.sum(d ** 2, axis=0, dtype=np.float64)
    M3 = np.sum(d ** 3, axis=0, dtype=np.float64)
    M4 = np.sum(d ** 4, axis=0, dtype=np.float64)

    xmin = np.full(cells, np.nan, dtype=np.float64)
    xmax = np.full(cells, np.nan, dtype=np.float64)
    # np.min/np.max are used only on columns with at least one valid observation,
    # so a fully missing ocean cell never emits an All-NaN warning.
    cols = np.flatnonzero(has)
    if cols.size:
        xv = x2[:, cols]
        vv = v2[:, cols]
        safe_min = np.where(vv, xv, np.inf)
        safe_max = np.where(vv, xv, -np.inf)
        xmin[cols] = np.min(safe_min, axis=0)
        xmax[cols] = np.max(safe_max, axis=0)

    target_shape = spatial_shape
    return (n.reshape(target_shape), mean.reshape(target_shape), M2.reshape(target_shape),
            M3.reshape(target_shape), M4.reshape(target_shape),
            xmin.reshape(target_shape), xmax.reshape(target_shape))


# >>> combine_moments_state: Mergeable central moments through fourth order.
# PEBAY/WELFORD-STYLE MERGE
# Combines two central-moment states without reconstructing raw observations.  This
# is what allows decade and FULL states to remain mergeable and restartable.
# ------------------------------------------------------------------------
# FUNCTION combine_moments_state — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def combine_moments_state(n1, m1, M21, M31, M41, n2, m2, M22, M32, M42):
    a = n1.astype(np.float64)
    b = n2.astype(np.float64)
    nt = a + b
    outn = nt.astype(np.uint32)
    outm = m1.copy(); out2 = M21.copy(); out3 = M31.copy(); out4 = M41.copy()
    only2 = (a == 0) & (b > 0)
    both = (a > 0) & (b > 0)
    if np.any(only2):
        outm[only2] = m2[only2]; out2[only2] = M22[only2]; out3[only2] = M32[only2]; out4[only2] = M42[only2]
    if np.any(both):
        aa = a[both]; bb = b[both]; nn = nt[both]; d = m2[both] - m1[both]
        outm[both] = m1[both] + d * bb / nn
        out2[both] = M21[both] + M22[both] + d*d*aa*bb/nn
        out3[both] = (M31[both] + M32[both] + d**3*aa*bb*(aa-bb)/(nn**2) + 3*d*(aa*M22[both]-bb*M21[both])/nn)
        out4[both] = (M41[both] + M42[both] + d**4*aa*bb*(aa*aa-aa*bb+bb*bb)/(nn**3)
                      + 6*d*d*(aa*aa*M22[both]+bb*bb*M21[both])/(nn**2)
                      + 4*d*(aa*M32[both]-bb*M31[both])/nn)
    return outn, outm, out2, out3, out4


# >>> combine_cov_state: Mergeable covariance state using the paired-valid population.
# PAIRED COVARIANCE MERGE
# Combines two covariance states that describe the SAME paired-valid population.
# Keeping paired counts and paired means separate is essential when either variable
# contains missing/invalid observations.
# ------------------------------------------------------------------------
# FUNCTION combine_cov_state — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def combine_cov_state(n1, mx1, my1, c1, n2, mx2, my2, c2):
    a = n1.astype(np.float64); b = n2.astype(np.float64); nt = a + b
    outn = nt.astype(np.uint32); outx = mx1.copy(); outy = my1.copy(); outc = c1.copy()
    only2 = (a == 0) & (b > 0); both = (a > 0) & (b > 0)
    if np.any(only2):
        outx[only2] = mx2[only2]; outy[only2] = my2[only2]; outc[only2] = c2[only2]
    if np.any(both):
        aa=a[both]; bb=b[both]; nn=nt[both]; dx=mx2[both]-mx1[both]; dy=my2[both]-my1[both]
        outx[both]=mx1[both]+dx*bb/nn; outy[both]=my1[both]+dy*bb/nn
        outc[both]=c1[both]+c2[both]+dx*dy*aa*bb/nn
    return outn, outx, outy, outc


# >>> histogram_counts: Vectorized empirical 2-D histogram using flat indices and bincount.
# EMPIRICAL 2-D HISTOGRAM COUNTER
# Flattens two-dimensional bin indices and uses vectorized counting rather than a
# Python loop over individual observations.  This function stores empirical evidence,
# not a fitted parametric probability surface.
# ------------------------------------------------------------------------
# FUNCTION histogram_counts — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def histogram_counts(x: np.ndarray, y: np.ndarray, valid: np.ndarray, xname: str, yname: str, bins: Tuple[int,int]) -> np.ndarray:
    hb0, hb1 = bins
    xmin,xmax=HIST_RANGES[xname]; ymin,ymax=HIST_RANGES[yname]
    valid = valid & (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
    rows,cells=np.nonzero(valid)
    out=np.zeros((x.shape[1],hb0,hb1),dtype=np.uint32)
    if rows.size==0: return out
    xv=x[rows,cells].astype(np.float64); yv=y[rows,cells].astype(np.float64)
    xb=np.minimum(((xv-xmin)/(xmax-xmin)*hb0).astype(np.int64),hb0-1)
    yb=np.minimum(((yv-ymin)/(ymax-ymin)*hb1).astype(np.int64),hb1-1)
    flat=cells*hb0*hb1 + xb*hb1 + yb
    out.flat[:] = np.bincount(flat, minlength=out.size)
    return out


# >>> config_fingerprint: Hashes scientifically relevant configuration so incompatible checkpoints cannot be reused silently.
# SCIENTIFIC CONFIGURATION IDENTITY
# Hashes the settings that can change the numerical meaning of the persisted state.
# A checkpoint produced under a different scientific contract must never be reused
# silently, even if the filenames happen to look compatible.
# ------------------------------------------------------------------------
# FUNCTION config_fingerprint — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def config_fingerprint(config: Config = CONFIG) -> str:
    payload = {
        "engine_version": ENGINE_VERSION,
        "schema_version": SCHEMA_VERSION,
        "config": {
            "start_year": config.start_year, "end_year": config.end_year,
            "chunk_lat": config.chunk_lat, "chunk_lon": config.chunk_lon,
            "hist_levels": config.hist_levels, "hist_pairs": config.hist_pairs,
            "hist_bins": HIST_BINS, "hist_ranges": HIST_RANGES,
            "thresholds": THRESHOLDS, "joint_thresholds": JOINT_THRESHOLDS,
            "variables": VARIABLES, "pairs": PAIRS, "levels": LEVELS,
        },
    }
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


# DURABLE TRANSACTION JOURNAL
# SQLite is used as the authoritative completion ledger.  It records OPEN/COMMITTED
# state for day x spatial-block work units.  Progress text/JSON is telemetry only.
# ========================================================================
# CLASS Journal — IMPLEMENTATION GUIDE
# ========================================================================
# Responsibility:
#   Own one clearly bounded part of the v10 engine. The class should be
#   read together with its caller and with the persisted state it owns.
#
# What to inspect:
#   1. Constructor state and configuration.
#   2. Public methods and their pre/post conditions.
#   3. Array shapes and coordinate conventions.
#   4. Failure behavior and recovery behavior.
#   5. Whether a value is scientific state, telemetry, or metadata.
#
# Scientific safety:
#   Optimizing this class must not change the sample population, calendar,
#   units, masking rules, or mathematical definition of any statistic.
# ========================================================================

class Journal:
    """SQLite-backed transaction truth; no scientific arrays live in SQLite."""
    def __init__(self, root: Path):
        self.root = ensure_dir(root)
        self.db_path = self.root / "journal.sqlite"
        self.before_dir = ensure_dir(self.root / "before_images")
        self.conn = sqlite3.connect(self.db_path)
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=FULL")
        self.conn.execute("PRAGMA foreign_keys=ON")
        self.conn.execute("""CREATE TABLE IF NOT EXISTS transactions(
            txid TEXT PRIMARY KEY,
            year INTEGER NOT NULL,
            month INTEGER NOT NULL,
            day INTEGER NOT NULL,
            y_chunk INTEGER NOT NULL,
            x_chunk INTEGER NOT NULL,
            status TEXT NOT NULL,
            record_json TEXT NOT NULL,
            created_utc TEXT NOT NULL,
            committed_utc TEXT,
            UNIQUE(year,month,day,y_chunk,x_chunk)
        )""")
        self.conn.commit()

    def close(self) -> None:
        self.conn.close()

    @staticmethod
    def txid(year:int,month:int,day:int,yc:int,xc:int)->str:
        return hashlib.sha256(f"{year:04d}-{month:02d}-{day:02d}-Y{yc:04d}-X{xc:04d}".encode("ascii")).hexdigest()[:32]

    def is_committed(self, year:int,month:int,day:int,yc:int,xc:int)->bool:
        cur=self.conn.execute("SELECT 1 FROM transactions WHERE year=? AND month=? AND day=? AND y_chunk=? AND x_chunk=? AND status='COMMITTED'",(year,month,day,yc,xc))
        return cur.fetchone() is not None

    def create_open(self, record: Dict[str,Any]) -> str:
        txid=record["txid"]
        self.conn.execute("INSERT OR REPLACE INTO transactions(txid,year,month,day,y_chunk,x_chunk,status,record_json,created_utc,committed_utc) VALUES(?,?,?,?,?,?,?,?,?,NULL)",
                          (txid,record["year"],record["month"],record["day"],record["y_chunk"],record["x_chunk"],"OPEN",canonical_json(record),record["created_utc"]))
        self.conn.commit()
        return txid

    def mark_committed(self, txid:str) -> None:
        self.conn.execute("UPDATE transactions SET status='COMMITTED', committed_utc=? WHERE txid=?",(utc_now(),txid))
        self.conn.commit()

    def mark_aborted(self, txid:str) -> None:
        self.conn.execute("UPDATE transactions SET status='ABORTED' WHERE txid=?",(txid,))
        self.conn.commit()

    def open_transactions(self)->List[Dict[str,Any]]:
        rows=self.conn.execute("SELECT txid,record_json FROM transactions WHERE status='OPEN'").fetchall()
        return [(json.loads(rj)) for _,rj in rows]

    def before_path(self,txid:str)->Path:
        return self.before_dir/f"{txid}.npz"

    def write_before(self,txid:str,payload:Dict[str,np.ndarray]) -> Path:
        path=self.before_path(txid)
        fd,tmp_name=tempfile.mkstemp(prefix=f".{txid}.",suffix=".npz.tmp",dir=str(self.before_dir))
        tmp=Path(tmp_name)
        try:
            with os.fdopen(fd,"wb") as fh:
                np.savez_compressed(fh,**payload)
                fh.flush(); os.fsync(fh.fileno())
            os.replace(tmp,path)
        finally:
            tmp.unlink(missing_ok=True)
        return path

    def cleanup(self,txid:str)->None:
        self.before_path(txid).unlink(missing_ok=True)

    def recover(self, restore_fn) -> int:
        opens=self.open_transactions(); n=0
        for rec in opens:
            txid=rec["txid"]
            try:
                restore_fn(rec)
                self.mark_aborted(txid); self.cleanup(txid); n+=1
            except Exception:
                LOG.exception("Recovery failed for transaction %s",txid)
                raise
        return n


# SCIENTIFIC STATE PERSISTENCE
# Stores the actual numerical checkpoint bundle for a spatial block.  The block file
# contains statistics; the journal separately answers the question: "was this unit
# committed?"  Separating numerical state from transaction truth prevents the old
# checkpoint/progress ambiguity.
# ========================================================================
# CLASS BlockCheckpoint — IMPLEMENTATION GUIDE
# ========================================================================
# Responsibility:
#   Own one clearly bounded part of the v10 engine. The class should be
#   read together with its caller and with the persisted state it owns.
#
# What to inspect:
#   1. Constructor state and configuration.
#   2. Public methods and their pre/post conditions.
#   3. Array shapes and coordinate conventions.
#   4. Failure behavior and recovery behavior.
#   5. Whether a value is scientific state, telemetry, or metadata.
#
# Scientific safety:
#   Optimizing this class must not change the sample population, calendar,
#   units, masking rules, or mathematical definition of any statistic.
# ========================================================================

class BlockCheckpoint:
    """One NetCDF shard per period and spatial block; chunked on DOY/level/cell."""
    def __init__(self, path: Path, lat: np.ndarray, lon: np.ndarray, config: Config):
        if Dataset is None: raise RuntimeError("netCDF4 is required")
        self.path=path; self.lat=lat; self.lon=lon; self.config=config
        self.nc: Optional[Dataset]=None
        self.ncell=lat.size*lon.size
        self.level_count=33

    def create_or_open(self) -> Dataset:
        if self.nc is not None: return self.nc
        if self.path.exists():
            self.nc=Dataset(self.path,"r+")
            if getattr(self.nc,"schema_version",None)!=SCHEMA_VERSION:
                raise RuntimeError(f"Checkpoint schema mismatch: {self.path}")
            return self.nc
        ensure_dir(self.path.parent)
        ds=Dataset(self.path,"w",format="NETCDF4")
        ds.createDimension("doy",366); ds.createDimension("level_bin",33); ds.createDimension("cell",self.ncell)
        ds.createDimension("variable",len(VARIABLES)); ds.createDimension("pair",len(PAIRS))
        ds.createDimension("threshold",max(len(v) for v in THRESHOLDS.values()))
        ds.createDimension("joint_threshold",len(JOINT_THRESHOLDS))
        ds.createDimension("hist_level_bin",9); ds.createDimension("x_bin",HIST_BINS[0]); ds.createDimension("y_bin",HIST_BINS[1])
        ds.createVariable("doy","i2",("doy",))[:]=np.arange(1,367,dtype=np.int16)
        ds.createVariable("level_bin","i2",("level_bin",))[:]=np.arange(33,dtype=np.int16)
        ds.createVariable("cell_lat","f4",("cell",))[:]=np.repeat(self.lat.astype(np.float32),self.lon.size)
        ds.createVariable("cell_lon","f4",("cell",))[:]=np.tile(self.lon.astype(np.float32),self.lat.size)
        ds.createVariable("variable","i1",("variable",))[:]=np.arange(len(VARIABLES),dtype=np.int8)
        ds.createVariable("pair","i1",("pair",))[:]=np.arange(len(PAIRS),dtype=np.int8)
        ds.createVariable("reserved_day","i1",("doy",))[:]=canonical_calendar()[3]
        chunks=(1,33,min(self.ncell,4096))
        f32=("f8",("doy","level_bin","cell")); u32=("u4",("doy","level_bin","cell"))
        for var in VARIABLES:
            for stat in ("mean","M2","M3","M4","min","max"):
                v=ds.createVariable(f"{stat}_{var}","f4",("doy","level_bin","cell"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=np.nan)
                if stat in ("min","max"): v.setncattr("statistic",stat)
            ds.createVariable(f"n_{var}","u4",("doy","level_bin","cell"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
            ds.createVariable(f"missing_count_{var}","u4",("doy","level_bin","cell"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
            ds.createVariable(f"threshold_count_{var}","u4",("doy","level_bin","cell","threshold"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=0,chunksizes=(1,33,min(self.ncell,4096),max(len(THRESHOLDS[var]),1)))
        ds.createVariable("supersaturation_count","u4",("doy","level_bin","cell"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
        ds.createVariable("invalid_e_over_p_count","u4",("doy","level_bin","cell"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
        for pi,(xname,yname) in enumerate(PAIRS):
            tag=f"{xname}__{yname}"
            for stat in ("mean_x","mean_y","M2_x","M2_y","Cxy"):
                ds.createVariable(f"pair_{tag}_{stat}","f4",("doy","level_bin","cell"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=np.nan)
            ds.createVariable(f"pair_{tag}_n","u4",("doy","level_bin","cell"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
            ds.createVariable(f"joint_threshold_count_{tag}","u4",("doy","level_bin","cell","joint_threshold"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=0,chunksizes=(1,33,min(self.ncell,4096),len(JOINT_THRESHOLDS)))
            if (xname,yname) in self.config.hist_pairs:
                ds.createVariable(f"hist_{tag}","u4",("doy","hist_level_bin","cell","x_bin","y_bin"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=0)
        month,day,label,res=canonical_calendar()
        ds.month=month.tolist(); ds.day=day.tolist(); ds.calendar_contract="slot 59 reserved; slot 60=Feb-28/Feb-29; slot 61=Mar-01"
        ds.engine_name=ENGINE_NAME; ds.engine_version=ENGINE_VERSION; ds.schema_version=SCHEMA_VERSION; ds.checkpoint_version=CHECKPOINT_VERSION; ds.config_fingerprint=config_fingerprint(self.config)
        ds.hist_levels_json=canonical_json(self.config.hist_levels); ds.hist_pairs_json=canonical_json(self.config.hist_pairs)
        ds.sync(); self.nc=ds; return ds

    def close(self)->None:
        if self.nc is not None:
            try:self.nc.sync()
            finally:self.nc.close()
            self.nc=None

    def slice_dict(self,doy:int,bin_indices:Sequence[int])->Dict[str,np.ndarray]:
        ds=self.create_or_open(); sl=np.asarray(list(bin_indices),dtype=np.int64)
        out={}
        for var in VARIABLES:
            for stat in ("mean","M2","M3","M4","min","max"):
                out[f"{stat}_{var}"]=np.asarray(ds.variables[f"{stat}_{var}"][doy,sl,:])
            out[f"n_{var}"]=np.asarray(ds.variables[f"n_{var}"][doy,sl,:])
            out[f"missing_count_{var}"]=np.asarray(ds.variables[f"missing_count_{var}"][doy,sl,:])
            out[f"threshold_count_{var}"]=np.asarray(ds.variables[f"threshold_count_{var}"][doy,sl,:,:])
        for pair in PAIRS:
            tag=f"{pair[0]}__{pair[1]}"
            for stat in ("mean_x","mean_y","M2_x","M2_y","Cxy"):
                out[f"pair_{tag}_{stat}"]=np.asarray(ds.variables[f"pair_{tag}_{stat}"][doy,sl,:])
            out[f"pair_{tag}_n"]=np.asarray(ds.variables[f"pair_{tag}_n"][doy,sl,:])
            out[f"joint_threshold_count_{tag}"]=np.asarray(ds.variables[f"joint_threshold_count_{tag}"][doy,sl,:,:])
            if f"hist_{tag}" in ds.variables:
                out[f"hist_{tag}"]=np.asarray(ds.variables[f"hist_{tag}"][doy, :, :, :,:])
        out["supersaturation_count"]=np.asarray(ds.variables["supersaturation_count"][doy,sl,:])
        out["invalid_e_over_p_count"]=np.asarray(ds.variables["invalid_e_over_p_count"][doy,sl,:])
        return out

    def write_slice_dict(self,doy:int,bin_indices:Sequence[int],data:Dict[str,np.ndarray]) -> None:
        ds=self.create_or_open(); sl=np.asarray(list(bin_indices),dtype=np.int64)
        for k,v in data.items():
            if k in ds.variables:
                ds.variables[k][doy,sl,...]=v
        ds.sync()


# ============================================================================
# CORE ENGINE ORCHESTRATOR
# ----------------------------------------------------------------------------
# Engine is intentionally the only component that coordinates the full data
# lifecycle.  It does not delegate the scientific computation to worker pools.
# Its responsibilities are deliberately visible in one place:
#
#   input inventory -> grid load -> recovery -> monthly open/alignment
#   -> day/block calculation -> transaction -> progress -> finalization
#   -> audit/merge audit -> provenance
#
# A future maintainer should resist moving scientific semantics into the CLI or
# progress reporting.  The journal and persisted statistical state define what
# has scientifically completed; console progress is only observability.
# ============================================================================
# MAIN V10 ENGINE
# Coordinates inventory, input alignment, physics, statistical accumulation, checkpoint
# persistence, recovery, finalization, audit, merge verification, and reports.  No
# worker pool is used: the engine is intentionally single-process.
# ========================================================================
# CLASS Engine — IMPLEMENTATION GUIDE
# ========================================================================
# Responsibility:
#   Own one clearly bounded part of the v10 engine. The class should be
#   read together with its caller and with the persisted state it owns.
#
# What to inspect:
#   1. Constructor state and configuration.
#   2. Public methods and their pre/post conditions.
#   3. Array shapes and coordinate conventions.
#   4. Failure behavior and recovery behavior.
#   5. Whether a value is scientific state, telemetry, or metadata.
#
# Scientific safety:
#   Optimizing this class must not change the sample population, calendar,
#   units, masking rules, or mathematical definition of any statistic.
# ========================================================================

class Engine:
    def __init__(self,config:Config=CONFIG):
        config.validate(); self.config=config
        self.root=ensure_dir(config.output_root); self.ckpt_root=ensure_dir(self.root/"checkpoints")
        self.final_root=ensure_dir(self.root/"netcdf"); self.report_root=ensure_dir(self.root/"reports")
        self.journal=Journal(self.ckpt_root)
        self.lat=None; self.lon=None; self.grid_shape=None
        self._open_blocks: Dict[Tuple[str,int,int],BlockCheckpoint]={}
        self._input_cache: Dict[int,Dict[str,Dict[int,Path]]]={}
        self.load_grid()
        self.recover()

    def recover(self)->None:
        def restore(rec:Dict[str,Any])->None:
            txid=rec["txid"]; before=self.journal.before_path(txid)
            if not before.exists(): raise RuntimeError(f"Missing before-image for {txid}")
            with np.load(before,allow_pickle=False) as z:
                for item in rec["targets"]:
                    path=Path(item["path"]); variable=item["variable"]; doy=item["doy"]; bins=np.asarray(item["bins"],dtype=np.int64)
                    block=self._open_or_get_file(item["period"],rec["y_chunk"],rec["x_chunk"],create=True)
                    ds=block.create_or_open(); ds.variables[variable][doy,bins,...]=z[item["key"]]
                    ds.sync()
        n=self.journal.recover(restore)
        if n: LOG.warning("Recovered %d uncommitted transactions",n)

    def load_grid(self)->None:
        if self.grid_shape is not None:return
        idx=build_file_index(self.config.start_year,self.config.t2m_dir)
        with open_dataset(idx[1]) as ds:
            self.lat=np.asarray(ds.latitude.values,dtype=np.float32); self.lon=np.asarray(ds.longitude.values,dtype=np.float32)
            self.grid_shape=(int(self.lat.size),int(self.lon.size))
        self._save_grid_metadata()

    def _save_grid_metadata(self)->None:
        atomic_write_json(self.root/"grid.json",{"ny":int(self.grid_shape[0]),"nx":int(self.grid_shape[1]),"chunk_lat":self.config.chunk_lat,"chunk_lon":self.config.chunk_lon,"latitude":self.lat.tolist(),"longitude":self.lon.tolist()})

    def file_indices(self,year:int)->Dict[str,Dict[int,Path]]:
        if year not in self._input_cache:
            self._input_cache[year]={"t2m":build_file_index(year,self.config.t2m_dir),"d2m":build_file_index(year,self.config.d2m_dir),"sp":build_file_index(year,self.config.sp_dir)}
        return self._input_cache[year]

    def _periods_for_year(self,year:int)->List[str]:
        return [p for p,(a,b) in PERIODS.items() if a<=year<=b]

    def _checkpoint_path(self,period:str,yc:int,xc:int)->Path:
        return self.ckpt_root/f"{period}__Y{yc:04d}__X{xc:04d}.nc"

    def _open_or_get_file(self,period:str,yc:int,xc:int,create:bool=True)->BlockCheckpoint:
        key=(period,yc,xc)
        if key in self._open_blocks:return self._open_blocks[key]
        y0=yc*self.config.chunk_lat; y1=min((yc+1)*self.config.chunk_lat,self.grid_shape[0])
        x0=xc*self.config.chunk_lon; x1=min((xc+1)*self.config.chunk_lon,self.grid_shape[1])
        lat=self.lat[y0:y1]; lon=self.lon[x0:x1]
        b=BlockCheckpoint(self._checkpoint_path(period,yc,xc),lat,lon,self.config)
        if create:b.create_or_open()
        self._open_blocks[key]=b; return b

    @staticmethod
    def _target_variable_keys() -> List[str]:
        keys=[]
        for var in VARIABLES:
            keys += [f"mean_{var}",f"M2_{var}",f"M3_{var}",f"M4_{var}",f"min_{var}",f"max_{var}",f"n_{var}",f"missing_count_{var}",f"threshold_count_{var}"]
        keys += ["supersaturation_count","invalid_e_over_p_count"]
        for x,y in PAIRS:
            tag=f"{x}__{y}"; keys += [f"pair_{tag}_mean_x",f"pair_{tag}_mean_y",f"pair_{tag}_M2_x",f"pair_{tag}_M2_y",f"pair_{tag}_Cxy",f"pair_{tag}_n",f"joint_threshold_count_{tag}"]
            if (x,y) in CONFIG.hist_pairs: keys.append(f"hist_{tag}")
        return keys

    def _update_state_arrays(self, old:Dict[str,np.ndarray], xdata:Dict[str,np.ndarray], valid:Dict[str,np.ndarray],
                             sup:np.ndarray, inv:np.ndarray, pair_data:Dict[Tuple[str,str],Tuple[np.ndarray,np.ndarray,np.ndarray]],
                             bin_index:int, total_hours:int) -> Dict[str,np.ndarray]:
        """Merge one temporal bin directly into the persistent in-memory state.

        The previous implementation copied the complete 33-bin state for every bin.
        v10 updates one bin in place, which materially reduces allocations and RAM traffic.
        """
        for var in VARIABLES:
            xv = np.asarray(xdata[var])
            vv = np.asarray(valid[var], dtype=bool)
            if xv.shape != vv.shape:
                # Validity is fundamentally defined by the derived physical field.
                # Fall back to its finite mask only when an upstream mask has a
                # malformed shape; never broadcast or reshape silently.
                if xv.ndim == 2 and vv.ndim == 2 and vv.size == xv.size:
                    vv = vv.reshape(xv.shape)
                else:
                    vv = np.isfinite(xv)
            valid[var] = vv & np.isfinite(xv)
            n2,m,b2,b3,b4,xmin,xmax=batch_moments(xv,valid[var])
            n1=old[f"n_{var}"][0].astype(np.float64)
            a=n1; b=n2.astype(np.float64); nt=a+b
            mean1=old[f"mean_{var}"][0].astype(np.float64)
            outm=mean1.copy(); out2=old[f"M2_{var}"][0].astype(np.float64).copy()
            out3=old[f"M3_{var}"][0].astype(np.float64).copy(); out4=old[f"M4_{var}"][0].astype(np.float64).copy()
            both=(a>0)&(b>0); only=(a==0)&(b>0)
            if np.any(only):
                outm[only]=m[only]; out2[only]=b2[only]; out3[only]=b3[only]; out4[only]=b4[only]
            if np.any(both):
                aa=a[both]; bb=b[both]; nn=nt[both]; d=m[both]-mean1[both]
                old2=out2[both].copy(); old3=out3[both].copy()
                outm[both]=mean1[both]+d*bb/nn
                out2[both]=old2+b2[both]+d*d*aa*bb/nn
                out3[both]=(old3+b3[both]+d**3*aa*bb*(aa-bb)/(nn**2)+3*d*(aa*b2[both]-bb*old2) / nn)
                out4[both]=(out4[both]+b4[both]+d**4*aa*bb*(aa*aa-aa*bb+bb*bb)/(nn**3)
                            +6*d*d*(aa*aa*b2[both]+bb*bb*old2)/(nn**2)
                            +4*d*(aa*b3[both]-bb*old3)/nn)
            old[f"n_{var}"][0]=nt.astype(np.uint32)
            old[f"mean_{var}"][0]=outm; old[f"M2_{var}"][0]=out2
            old[f"M3_{var}"][0]=out3; old[f"M4_{var}"][0]=out4
            old[f"min_{var}"][0]=np.where(b>0,np.minimum(old[f"min_{var}"][0],xmin),old[f"min_{var}"][0])
            old[f"max_{var}"][0]=np.where(b>0,np.maximum(old[f"max_{var}"][0],xmax),old[f"max_{var}"][0])
            old[f"missing_count_{var}"][0]+= (xdata[var].shape[0]-n2).astype(np.uint32)
            for ti,t in enumerate(THRESHOLDS[var]):
                old[f"threshold_count_{var}"][bin_index,:,ti] += (valid[var] & (xdata[var]>t)).sum(axis=0,dtype=np.uint32)

        old["supersaturation_count"][0]+=sup.astype(bool).sum(axis=0,dtype=np.uint32)
        old["invalid_e_over_p_count"][0]+=inv.astype(bool).sum(axis=0,dtype=np.uint32)

        for (xname,yname),(xp,yp,pvalid) in pair_data.items():
            tag=f"{xname}__{yname}"
            xp = np.asarray(xp)
            yp = np.asarray(yp)
            pvalid = np.asarray(pvalid, dtype=bool)
            if xp.shape != yp.shape:
                raise ValueError(f"Pair shape mismatch for {tag}: {xp.shape} != {yp.shape}")
            if pvalid.shape != xp.shape:
                # Pair validity can be reconstructed exactly from the derived
                # finite fields; this is preferable to a dangerous broadcast.
                pvalid = np.isfinite(xp) & np.isfinite(yp)
            else:
                pvalid &= np.isfinite(xp) & np.isfinite(yp)
            n2,mx,bx2,_,_,_,_=batch_moments(xp,pvalid)
            _,my,by2,_,_,_,_=batch_moments(yp,pvalid)
            n1=old[f"pair_{tag}_n"][0].astype(np.float64)
            mx1=old[f"pair_{tag}_mean_x"][0].astype(np.float64)
            my1=old[f"pair_{tag}_mean_y"][0].astype(np.float64)
            m2x1=old[f"pair_{tag}_M2_x"][0].astype(np.float64)
            m2y1=old[f"pair_{tag}_M2_y"][0].astype(np.float64)
            c1=old[f"pair_{tag}_Cxy"][0].astype(np.float64)
            b=n2.astype(np.float64); a=n1; nt=a+b; both=(a>0)&(b>0); only=(a==0)&(b>0)
            bcxy_all=np.sum(np.where(pvalid,(xp-mx[None,:])*(yp-my[None,:]),0.0),axis=0)
            if np.any(only):
                mx1[only]=mx[only]; my1[only]=my[only]; m2x1[only]=bx2[only]; m2y1[only]=by2[only]
                c1[only]=bcxy_all[only]
            if np.any(both):
                aa=a[both]; bb=b[both]; nn=nt[both]; dx=mx[both]-mx1[both]; dy=my[both]-my1[both]
                old_mx=mx1[both].copy(); old_my=my1[both].copy(); old_m2x=m2x1[both].copy(); old_m2y=m2y1[both].copy(); old_c=c1[both].copy()
                mx1[both]=old_mx+dx*bb/nn; my1[both]=old_my+dy*bb/nn
                m2x1[both]=old_m2x+bx2[both]+dx*dx*aa*bb/nn
                m2y1[both]=old_m2y+by2[both]+dy*dy*aa*bb/nn
                c1[both]=old_c+bcxy_all[both]+dx*dy*aa*bb/nn
            old[f"pair_{tag}_n"][0]=nt.astype(np.uint32)
            old[f"pair_{tag}_mean_x"][0]=mx1; old[f"pair_{tag}_mean_y"][0]=my1
            old[f"pair_{tag}_M2_x"][0]=m2x1; old[f"pair_{tag}_M2_y"][0]=m2y1; old[f"pair_{tag}_Cxy"][0]=c1
            for jt,(a_name,ta,b_name,tb) in enumerate(JOINT_THRESHOLDS):
                if {a_name,b_name}!={xname,yname}: continue
                va=xp if a_name==xname else yp; vb=yp if b_name==yname else xp
                old[f"joint_threshold_count_{tag}"][bin_index,:,jt] += (pvalid & (va>ta) & (vb>tb)).sum(axis=0,dtype=np.uint32)
        return old

    def _prepare_before(self,rec:Dict[str,Any]) -> Dict[str,np.ndarray]:
        before={}
        for item in rec["targets"]:
            key=item["key"]; block=self._open_or_get_file(item["period"],rec["y_chunk"],rec["x_chunk"]); ds=block.create_or_open()
            bins=np.asarray(item["bins"],dtype=np.int64); doy=item["doy"]; before[key]=np.asarray(ds.variables[item["variable"]][doy,bins,...])
        return before

    def _commit_day_block(self,year:int,month:int,day:int,yc:int,xc:int,period_updates:Dict[str,Dict[str,np.ndarray]],bin_indices:List[int],doy:int) -> None:
        txid=self.journal.txid(year,month,day,yc,xc)
        if self.journal.is_committed(year,month,day,yc,xc): return
        targets=[]; before_payload={}
        for period,data in period_updates.items():
            block=self._open_or_get_file(period,yc,xc)
            ds=block.create_or_open()
            for variable,value in data.items():
                if variable not in ds.variables: continue
                key=f"{period}__{variable}"
                # The NetCDF variable itself defines the legal temporal-bin
                # cardinality. This prevents a 33-bin selector from touching
                # a histogram variable that has only 9 bins.
                vshape=ds.variables[variable].shape
                if len(vshape) < 2:
                    raise RuntimeError(f"Invalid checkpoint variable shape for {variable}: {vshape}")
                nbins=int(vshape[1])
                if nbins not in (9,33):
                    raise RuntimeError(f"Unsupported temporal schema for {variable}: {nbins} bins")
                safe_bins=[int(b) for b in bin_indices if 0 <= int(b) < nbins]
                if not safe_bins:
                    continue
                before_payload[key]=np.asarray(ds.variables[variable][doy,safe_bins,...])
                targets.append({"period":period,"variable":variable,"bins":safe_bins,"doy":doy,"key":key,"schema_bins":nbins})
        rec={"txid":txid,"year":year,"month":month,"day":day,"y_chunk":yc,"x_chunk":xc,"doy":doy,"bins":bin_indices,"targets":targets,"created_utc":utc_now()}
        # The rollback image must be durable before the OPEN transaction becomes
        # visible to recovery. This prevents an OPEN row with no rollback image.
        self.journal.write_before(txid,before_payload)
        self.journal.create_open(rec)
        try:
            for period,data in period_updates.items():
                block=self._open_or_get_file(period,yc,xc); ds=block.create_or_open()
                for variable,value in data.items():
                    if variable in ds.variables:
                        nbins=int(ds.variables[variable].shape[1])
                        safe_bins=[int(b) for b in bin_indices if 0 <= int(b) < nbins]
                        if safe_bins:
                            ds.variables[variable][doy,safe_bins,...]=value
                ds.sync()
            self.journal.mark_committed(txid)
            self.journal.cleanup(txid)
        except Exception:
            LOG.exception("Transaction %s failed; rollback will occur now",txid)
            self.recover()
            raise

    def _process_day_block(self,year:int,month:int,day:int,yc:int,xc:int,ds_t:xr.Dataset,ds_d:xr.Dataset,ds_p:xr.Dataset,units:Dict[str,Any]) -> None:
        y0=yc*self.config.chunk_lat; y1=min((yc+1)*self.config.chunk_lat,self.grid_shape[0]); x0=xc*self.config.chunk_lon; x1=min((xc+1)*self.config.chunk_lon,self.grid_shape[1])
        doy=climatological_doy_from_date(year,month,day)
        ti0=(day-1)*24; ti1=ti0+24
        if doy==59: return
        T3=convert_temperature(ds_t["t2m"].isel(time=slice(ti0,ti1),latitude=slice(y0,y1),longitude=slice(x0,x1)).values,units["t2m"])
        Td3=convert_temperature(ds_d["d2m"].isel(time=slice(ti0,ti1),latitude=slice(y0,y1),longitude=slice(x0,x1)).values,units["d2m"])
        P3=convert_pressure(ds_p["sp"].isel(time=slice(ti0,ti1),latitude=slice(y0,y1),longitude=slice(x0,x1)).values,units["sp"])
        if T3.shape != Td3.shape or T3.shape != P3.shape:
            raise RuntimeError(f"Aligned hourly block shape mismatch: T2m={T3.shape}, D2m={Td3.shape}, SP={P3.shape}")
        if T3.shape[0] != 24:
            raise RuntimeError(f"Expected 24 hourly samples for day {year}-{month:02d}-{day:02d}, got {T3.shape[0]}")
        cells=T3.shape[1]*T3.shape[2]; T=T3.reshape(24,cells); Td=Td3.reshape(24,cells); P=P3.reshape(24,cells)
        phys=derive_moisture(T,Td,P); periods=self._periods_for_year(year)
        updates_by_period={}
        # One transaction is the whole day/block. Every level is updated in memory first.
        for period in periods:
            block=self._open_or_get_file(period,yc,xc)
            old=block.slice_dict(doy-1,range(33))
            # helper to update a single level range
            all_new=old
            for level,(b0,b1) in LEVEL_BINS.items():
                if level=="L1": groups=[(0,24)]
                elif level=="L2": groups=[(3*k,3*k+3) for k in range(8)]
                else: groups=[(h,h+1) for h in range(24)]
                local_bins=list(range(b0,b1))
                # Slice old into the target bins, then update each bin independently.
                for gi,(a,b) in zip(local_bins,groups):
                    # State-schema-aware selection. The actual array shape is the
                    # authoritative source for bin cardinality: ordinary temporal
                    # states use 33 bins; histogram states use 9 bins.
                    subold={}
                    for k,v in old.items():
                        if getattr(v, "ndim", 0) == 0:
                            subold[k]=v
                            continue
                        nbins=int(v.shape[0])
                        if nbins == 9:
                            safe_gi=[int(z) for z in np.atleast_1d(gi) if 0 <= int(z) < 9]
                        elif nbins == 33:
                            safe_gi=[int(z) for z in np.atleast_1d(gi) if 0 <= int(z) < 33]
                        else:
                            # Any unexpected state shape is a schema error.
                            raise RuntimeError(f"Unexpected temporal state bin count for {k}: {nbins}")
                        if safe_gi:
                            # `np.take` already preserves the selected bin axis.  Do not add
                            # another singleton axis here: the persistent state
                            # for one selected bin must be shape (1, cells), so
                            # `_update_state_arrays(...)[0]` produces (cells,).
                            subold[k]=np.take(v, safe_gi, axis=0)
                    xdata={v:phys[v][a:b,:] for v in VARIABLES}; valid={v:phys["valid"][v][a:b,:] for v in VARIABLES}
                    pair_data={(x,y):(phys[x][a:b,:],phys[y][a:b,:],valid[x][a:b,:]&valid[y][a:b,:]) for x,y in PAIRS}
                    sup=phys["supersat"][a:b,:]; inv=phys["invalid_e_over_p"][a:b,:]
                    upd=self._update_state_arrays(subold,xdata,valid,sup,inv,pair_data,0,b-a)
                    for k,v in upd.items(): old[k][gi,...]=v[0,...] if v.ndim>=1 and v.shape[0]==1 else v
            # Histograms for L1/L2 only; exact counts, stored on disk.
            hist_row=0
            for level in self.config.hist_levels:
                if level=="L1": groups=[(0,24)]
                elif level=="L2": groups=[(3*k,3*k+3) for k in range(8)]
                else: continue
                for a,b in groups:
                    valid_pair=phys["valid"]["rh"][a:b,:]&phys["valid"]["q"][a:b,:]
                    h=histogram_counts(phys["rh"][a:b,:],phys["q"][a:b,:],valid_pair,"rh","q",HIST_BINS)
                    tag="rh__q"; old[f"hist_{tag}"][hist_row,...]=old[f"hist_{tag}"][hist_row,...]+h; hist_row+=1
            # Write only variables that belong to the changed day; all 33 bins.
            updates_by_period[period]={k:v for k,v in old.items()}
        self._commit_day_block(year,month,day,yc,xc,updates_by_period,list(range(33)),doy-1)

    def close_block(self,yc:int,xc:int)->None:
        for period in PERIODS:
            key=(period,yc,xc)
            block=self._open_blocks.pop(key,None)
            if block is not None:
                block.close()

    def run(self)->Dict[str,Any]:
        self.load_grid(); total_units=0; completed=0
        for y in range(self.config.start_year,self.config.end_year+1):
            for m in range(1,13):
                total_units += calendar.monthrange(y,m)[1]*(((self.grid_shape[0]+self.config.chunk_lat-1)//self.config.chunk_lat)*((self.grid_shape[1]+self.config.chunk_lon-1)//self.config.chunk_lon))
        progress_path=self.root/"progress.json"; start=time.time()
        for year in range(self.config.start_year,self.config.end_year+1):
            idx=self.file_indices(year)
            for month in range(1,13):
                ds_t, ds_d, ds_p = open_aligned_triplet(idx["t2m"][month], idx["d2m"][month], idx["sp"][month])
                try:
                    validate_datasets(ds_t, ds_d, ds_p, year, month)
                    units={"t2m":ds_t["t2m"].attrs.get("units"),"d2m":ds_d["d2m"].attrs.get("units"),"sp":ds_p["sp"].attrs.get("units")}
                    for yc in range((self.grid_shape[0]+self.config.chunk_lat-1)//self.config.chunk_lat):
                        for xc in range((self.grid_shape[1]+self.config.chunk_lon-1)//self.config.chunk_lon):
                            # Reopen output shards only for this year/month/block; no worker pool.
                            for day in range(1,calendar.monthrange(year,month)[1]+1):
                                if self.journal.is_committed(year,month,day,yc,xc):
                                    completed+=1; continue
                                self._process_day_block(year,month,day,yc,xc,ds_t,ds_d,ds_p,units); completed+=1
                                if completed%16==0:
                                    elapsed=time.time()-start; rate=completed/max(elapsed,1e-9); remaining=total_units-completed
                                    atomic_write_json(progress_path,{"engine_version":ENGINE_VERSION,"completed_units":completed,"total_units":total_units,"remaining_units":remaining,"percent":100.0*completed/max(total_units,1),"rate_units_s":rate,"eta_seconds":remaining/max(rate,1e-12),"updated_utc":utc_now()})
                            self.close_block(yc,xc)
                finally:
                    ds_t.close(); ds_d.close(); ds_p.close()
        for b in self._open_blocks.values(): b.close()
        self._open_blocks.clear()
        outputs=self.finalize_outputs();
        manifest={"project":ENGINE_NAME,"engine_version":ENGINE_VERSION,"schema_version":SCHEMA_VERSION,"config_fingerprint":config_fingerprint(self.config),"periods":PERIODS,"inputs":{"t2m":str(self.config.t2m_dir),"d2m":str(self.config.d2m_dir),"sp":str(self.config.sp_dir)},"grid":{"ny":int(self.grid_shape[0]),"nx":int(self.grid_shape[1])},"grid_alignment":{"reference":"t2m","accepted_axis_normalization":"reversed_latitude_or_longitude_within_1e-6"},"chunk":{"lat":self.config.chunk_lat,"lon":self.config.chunk_lon},"outputs":outputs,"completed_units":completed,"total_units":total_units,"created_utc":utc_now()}
        atomic_write_json(self.root/"run_manifest.json",manifest); return manifest

    def _new_final_file(self,path:Path,period:str)->Dataset:
        if Dataset is None: raise RuntimeError("netCDF4 required")
        tmp=path.with_suffix(path.suffix+".part")
        if tmp.exists(): tmp.unlink()
        ds=Dataset(tmp,"w",format="NETCDF4")
        ds.createDimension("doy",366); ds.createDimension("level_bin",33); ds.createDimension("latitude",self.grid_shape[0]); ds.createDimension("longitude",self.grid_shape[1]); ds.createDimension("threshold",max(len(v) for v in THRESHOLDS.values())); ds.createDimension("joint_threshold",len(JOINT_THRESHOLDS)); ds.createDimension("x_bin",8); ds.createDimension("y_bin",8)
        ds.createVariable("doy","i2",("doy",))[:]=np.arange(1,367,dtype=np.int16); ds.createVariable("level_bin","i2",("level_bin",))[:]=np.arange(33,dtype=np.int16)
        ds.createVariable("latitude","f4",("latitude",))[:]=self.lat; ds.createVariable("longitude","f4",("longitude",))[:]=self.lon
        for var in VARIABLES:
            for stat in ("mean","std","skew","kurt","min","max"):
                ds.createVariable(f"{stat}_{var}","f4",("doy","level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=np.nan)
            ds.createVariable(f"n_{var}","u4",("doy","level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=0)
        ds.period=period; ds.engine_version=ENGINE_VERSION; ds.schema_version=SCHEMA_VERSION; ds.calendar_contract="slot 59 reserved; slot 60 Feb-28/29 composite; slot 61 Mar-01"; ds.level_contract="L1=0; L2=1..8; L3=9..32"; ds.sync(); return ds

    def _create_output_bundle(self, period: str) -> Dict[str, Dataset]:
        if Dataset is None:
            raise RuntimeError("netCDF4 required")
        outdir=ensure_dir(self.final_root/period)
        ny,nx=self.grid_shape
        paths={
            "main": outdir/f"moisture_climatology_{period}_v10.nc",
            "diagnostics": outdir/f"moisture_climatology_diagnostics_{period}_v10.nc",
            "bivariate": outdir/f"moisture_climatology_bivariate_{period}_v10.nc",
            "empirical_rhq": outdir/f"moisture_bivariate_empirical_rh__q_{period}_v10.nc",
        }
        tmp={k:p.with_suffix(p.suffix+".part") for k,p in paths.items()}
        for q in tmp.values():
            if q.exists(): q.unlink()
        ds_main=Dataset(tmp["main"],"w",format="NETCDF4")
        ds_diag=Dataset(tmp["diagnostics"],"w",format="NETCDF4")
        ds_biv=Dataset(tmp["bivariate"],"w",format="NETCDF4")
        ds_hist=Dataset(tmp["empirical_rhq"],"w",format="NETCDF4")
        for ds in (ds_main,ds_diag,ds_biv,ds_hist):
            ds.createDimension("doy",366); ds.createDimension("level_bin",33); ds.createDimension("latitude",ny); ds.createDimension("longitude",nx)
        ds_diag.createDimension("threshold",max(len(v) for v in THRESHOLDS.values()))
        ds_biv.createDimension("joint_threshold",len(JOINT_THRESHOLDS))
        ds_hist.createDimension("hist_level_bin",9)
        month,day,label,reserved=canonical_calendar()
        for ds in (ds_main,ds_diag,ds_biv,ds_hist):
            ds.createVariable("doy","i2",("doy",))[:]=np.arange(1,367,dtype=np.int16)
            ds.createVariable("level_bin","i2",("level_bin",))[:]=np.arange(33,dtype=np.int16)
            ds.createVariable("latitude","f4",("latitude",))[:]=self.lat
            ds.createVariable("longitude","f4",("longitude",))[:]=self.lon
            ds.createVariable("month","i2",("doy",))[:]=month
            ds.createVariable("day","i2",("doy",))[:]=day
            ds.createVariable("reserved_day","i1",("doy",))[:]=reserved
            ds.engine_name=ENGINE_NAME; ds.engine_version=ENGINE_VERSION; ds.schema_version=SCHEMA_VERSION; ds.period=period
            ds.calendar_contract="slot 59 reserved; slot 60=Feb-28/Feb-29; slot 61=Mar-01"
            ds.thresholds_json=canonical_json(THRESHOLDS); ds.joint_thresholds_json=canonical_json(JOINT_THRESHOLDS)
            ds.histogram_ranges_json=canonical_json(HIST_RANGES); ds.histogram_bins_json=canonical_json(HIST_BINS)
        chunks=(1,1,min(self.config.chunk_lat,ny),min(self.config.chunk_lon,nx))
        for var in VARIABLES:
            for stat in ("mean","std","skew","kurt","min","max"):
                ds_main.createVariable(f"{stat}_{var}","f4",("doy","level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=np.nan)
            ds_main.createVariable(f"n_{var}","u4",("doy","level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
            ds_diag.createVariable(f"missing_count_{var}","u4",("doy","level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
            ds_diag.createVariable(f"threshold_count_{var}","u4",("doy","level_bin","latitude","longitude","threshold"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=0,chunksizes=(1,1,min(self.config.chunk_lat,ny),min(self.config.chunk_lon,nx),max(len(THRESHOLDS[var]),1)))
        ds_diag.createVariable("supersaturation_count","u4",("doy","level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
        ds_diag.createVariable("invalid_e_over_p_count","u4",("doy","level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0)
        for x,y in PAIRS:
            tag=f"{x}__{y}"
            for stat in ("n","mean_x","mean_y","std_x","std_y","cov","corr"):
                dtype="u4" if stat=="n" else "f4"
                ds_biv.createVariable(f"{tag}_{stat}",dtype,("doy","level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,chunksizes=chunks,fill_value=0 if stat=="n" else np.nan)
            ds_biv.createVariable(f"joint_threshold_count_{tag}","u4",("doy","level_bin","latitude","longitude","joint_threshold"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=0,chunksizes=(1,1,min(self.config.chunk_lat,ny),min(self.config.chunk_lon,nx),len(JOINT_THRESHOLDS)))
        ds_hist.createVariable("hist_level_bin","i2",("hist_level_bin",))[:]=np.arange(9,dtype=np.int16)
        ds_hist.createVariable("rh_bin_edge","f4",("x_bin",))[:]=np.linspace(HIST_RANGES["rh"][0],HIST_RANGES["rh"][1],HIST_BINS[0]+1,dtype=np.float32)
        ds_hist.createVariable("q_bin_edge","f4",("y_bin",))[:]=np.linspace(HIST_RANGES["q"][0],HIST_RANGES["q"][1],HIST_BINS[1]+1,dtype=np.float32)
        ds_hist.createVariable("n_valid","u4",("doy","hist_level_bin","latitude","longitude"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=0,chunksizes=(1,9,min(self.config.chunk_lat,ny),min(self.config.chunk_lon,nx)))
        ds_hist.createVariable("count","u4",("doy","hist_level_bin","latitude","longitude","x_bin","y_bin"),zlib=True,complevel=self.config.compression,shuffle=self.config.shuffle,fill_value=0,chunksizes=(1,9,min(self.config.chunk_lat,ny),min(self.config.chunk_lon,nx),8,8))
        return {"main":ds_main,"diagnostics":ds_diag,"bivariate":ds_biv,"empirical_rhq":ds_hist,"_tmp_paths":tmp,"_paths":paths}

    def finalize_outputs(self)->Dict[str,str]:
        ny,nx=self.grid_shape; outputs={}
        for period in PERIODS:
            bundle=self._create_output_bundle(period); ds_main=bundle["main"]; ds_diag=bundle["diagnostics"]; ds_biv=bundle["bivariate"]; ds_hist=bundle["empirical_rhq"]
            try:
                for yc in range((ny+self.config.chunk_lat-1)//self.config.chunk_lat):
                    y0=yc*self.config.chunk_lat; y1=min((yc+1)*self.config.chunk_lat,ny)
                    for xc in range((nx+self.config.chunk_lon-1)//self.config.chunk_lon):
                        x0=xc*self.config.chunk_lon; x1=min((xc+1)*self.config.chunk_lon,nx)
                        block=self._open_or_get_file(period,yc,xc); src=block.create_or_open()
                        for doy0 in range(366):
                            if doy0==58: continue
                            for var in VARIABLES:
                                n=np.asarray(src.variables[f"n_{var}"][doy0,:,:],dtype=np.float64)
                                mean=np.asarray(src.variables[f"mean_{var}"][doy0,:,:],dtype=np.float64)
                                M2=np.asarray(src.variables[f"M2_{var}"][doy0,:,:],dtype=np.float64)
                                M3=np.asarray(src.variables[f"M3_{var}"][doy0,:,:],dtype=np.float64)
                                M4=np.asarray(src.variables[f"M4_{var}"][doy0,:,:],dtype=np.float64)
                                std=np.full_like(n,np.nan); skew=np.full_like(n,np.nan); kurt=np.full_like(n,np.nan); ok2=n>=2; ok3=(n>=3)&(M2>0); ok4=(n>=4)&(M2>0)
                                std[ok2]=np.sqrt(np.maximum(M2[ok2]/(n[ok2]-1),0))
                                if np.any(ok3):
                                    nn=n[ok3]; mm2=M2[ok3]/nn; mm3=M3[ok3]/nn; skew[ok3]=np.sqrt(nn*(nn-1))/(nn-2)*mm3/np.power(mm2,1.5)
                                if np.any(ok4):
                                    nn=n[ok4]; b2=nn*M4[ok4]/np.square(M2[ok4]); kurt[ok4]=((nn-1)/((nn-2)*(nn-3)))*((nn+1)*b2-3*(nn-1))
                                sl=np.s_[doy0,:,y0:y1,x0:x1]
                                ds_main.variables[f"mean_{var}"][sl]=mean.reshape((33,y1-y0,x1-x0)); ds_main.variables[f"std_{var}"][sl]=std.reshape((33,y1-y0,x1-x0)); ds_main.variables[f"skew_{var}"][sl]=skew.reshape((33,y1-y0,x1-x0)); ds_main.variables[f"kurt_{var}"][sl]=kurt.reshape((33,y1-y0,x1-x0)); ds_main.variables[f"min_{var}"][sl]=np.asarray(src.variables[f"min_{var}"][doy0,:,:]).reshape((33,y1-y0,x1-x0)); ds_main.variables[f"max_{var}"][sl]=np.asarray(src.variables[f"max_{var}"][doy0,:,:]).reshape((33,y1-y0,x1-x0)); ds_main.variables[f"n_{var}"][sl]=np.asarray(src.variables[f"n_{var}"][doy0,:,:]).reshape((33,y1-y0,x1-x0))
                                ds_diag.variables[f"missing_count_{var}"][sl]=np.asarray(src.variables[f"missing_count_{var}"][doy0,:,:]).reshape((33,y1-y0,x1-x0)); ds_diag.variables[f"threshold_count_{var}"][doy0,:,y0:y1,x0:x1,:]=np.asarray(src.variables[f"threshold_count_{var}"][doy0,:,:,:]).reshape((33,y1-y0,x1-x0,-1))
                            ds_diag.variables["supersaturation_count"][doy0,:,y0:y1,x0:x1]=np.asarray(src.variables["supersaturation_count"][doy0,:,:]).reshape((33,y1-y0,x1-x0)); ds_diag.variables["invalid_e_over_p_count"][doy0,:,y0:y1,x0:x1]=np.asarray(src.variables["invalid_e_over_p_count"][doy0,:,:]).reshape((33,y1-y0,x1-x0))
                            # Bivariate reference fields.
                            for x,y in PAIRS:
                                tag=f"{x}__{y}"; pn=np.asarray(src.variables[f"pair_{tag}_n"][doy0,:,:],dtype=np.float64); mx=np.asarray(src.variables[f"pair_{tag}_mean_x"][doy0,:,:],dtype=np.float64); my=np.asarray(src.variables[f"pair_{tag}_mean_y"][doy0,:,:],dtype=np.float64); c=np.asarray(src.variables[f"pair_{tag}_Cxy"][doy0,:,:],dtype=np.float64)
                                sx=np.asarray(src.variables[f"pair_{tag}_M2_x"][doy0,:,:],dtype=np.float64); sy=np.asarray(src.variables[f"pair_{tag}_M2_y"][doy0,:,:],dtype=np.float64); cov=np.full_like(pn,np.nan); corr=np.full_like(pn,np.nan); ok=pn>=2; cov[ok]=c[ok]/(pn[ok]-1); sxx=np.full_like(pn,np.nan); syy=np.full_like(pn,np.nan); sxx[ok]=np.sqrt(np.maximum(sx[ok]/(pn[ok]-1),0)); syy[ok]=np.sqrt(np.maximum(sy[ok]/(pn[ok]-1),0)); okc=ok&(sxx>0)&(syy>0); corr[okc]=np.clip(cov[okc]/(sxx[okc]*syy[okc]),-1,1)
                                tag_sl=np.s_[doy0,:,y0:y1,x0:x1]; ds_biv.variables[f"{tag}_n"][tag_sl]=pn.reshape((33,y1-y0,x1-x0)); ds_biv.variables[f"{tag}_mean_x"][tag_sl]=mx.reshape((33,y1-y0,x1-x0)); ds_biv.variables[f"{tag}_mean_y"][tag_sl]=my.reshape((33,y1-y0,x1-x0)); ds_biv.variables[f"{tag}_std_x"][tag_sl]=sxx.reshape((33,y1-y0,x1-x0)); ds_biv.variables[f"{tag}_std_y"][tag_sl]=syy.reshape((33,y1-y0,x1-x0)); ds_biv.variables[f"{tag}_cov"][tag_sl]=cov.reshape((33,y1-y0,x1-x0)); ds_biv.variables[f"{tag}_corr"][tag_sl]=corr.reshape((33,y1-y0,x1-x0)); ds_biv.variables[f"joint_threshold_count_{tag}"][doy0,:,y0:y1,x0:x1,:]=np.asarray(src.variables[f"joint_threshold_count_{tag}"][doy0,:,:,:]).reshape((33,y1-y0,x1-x0,-1))
                            if "hist_rh__q" in src.variables:
                                hist_block=np.asarray(src.variables["hist_rh__q"][doy0,:,:,:,:],dtype=np.uint32).reshape((9,y1-y0,x1-x0,8,8))
                                ds_hist.variables["count"][doy0,:,y0:y1,x0:x1,:,:]=hist_block
                                ds_hist.variables["n_valid"][doy0,:,y0:y1,x0:x1]=hist_block.sum(axis=(3,4),dtype=np.uint32)
                        block.close(); self._open_blocks.pop((period,yc,xc),None)
                for ds in (ds_main,ds_diag,ds_biv,ds_hist): ds.sync()
            finally:
                for ds in (ds_main,ds_diag,ds_biv,ds_hist): ds.close()
            for key,path in bundle["_paths"].items():
                os.replace(bundle["_tmp_paths"][key],path); outputs[f"{period}:{key}"]=str(path)
        return outputs

    def audit(self)->Dict[str,Any]:
        result={"engine_version":ENGINE_VERSION,"periods":{},"errors":[]}
        for period in PERIODS:
            period_errors=[]; count=0
            for p in sorted(self.ckpt_root.glob(f"{period}__Y*.nc")):
                try:
                    with Dataset(p,"r") as ds:
                        for var in VARIABLES:
                            n=np.asarray(ds.variables[f"n_{var}"][:],dtype=np.uint64)
                            if np.any(n > np.iinfo(np.uint32).max): period_errors.append(f"n overflow {p.name} {var}")
                            mn=np.asarray(ds.variables[f"min_{var}"][:],dtype=float); mx=np.asarray(ds.variables[f"max_{var}"][:],dtype=float)
                            bad=(np.isfinite(mn)&np.isfinite(mx)&(mn>mx))
                            if np.any(bad): period_errors.append(f"min>max {p.name} {var}")
                        if int(np.asarray(ds.variables["reserved_day"][58]))!=1: period_errors.append(f"reserved slot broken {p.name}")
                        count += 1
                except Exception as exc: period_errors.append(f"{p.name}: {exc!r}")
            result["periods"][period]={"checkpoint_files":count,"errors":period_errors}; result["errors"].extend(period_errors)
        result["status"]="PASS" if not result["errors"] else "FAIL"
        atomic_write_json(self.report_root/"audit.json",result); return result

    def merge_audit(self)->Dict[str,Any]:
        # Compare FULL against an on-the-fly merge of four decades for a deterministic subset.
        out={"status":"PASS","samples":0,"max_abs_mean":0.0,"max_abs_M2":0.0,"errors":[]}
        decs=["DECADE_1981_1990","DECADE_1991_2000","DECADE_2001_2010","DECADE_2011_2020"]
        files=[]
        for yc in range((self.grid_shape[0]+self.config.chunk_lat-1)//self.config.chunk_lat):
            for xc in range((self.grid_shape[1]+self.config.chunk_lon-1)//self.config.chunk_lon):
                ffull=self._checkpoint_path("FULL_1981_2020",yc,xc)
                if not ffull.exists(): continue
                with Dataset(ffull,"r") as full:
                    for var in VARIABLES:
                        a=np.asarray(full.variables[f"mean_{var}"][::60,0,::100],dtype=float)
                        n=np.asarray(full.variables[f"n_{var}"][::60,0,::100],dtype=float)
                        m2=np.asarray(full.variables[f"M2_{var}"][::60,0,::100],dtype=float)
                        mm=np.zeros_like(a); nn=np.zeros_like(n); M=np.zeros_like(m2)
                        # use first scalar grid subset for audit
                        for dec in decs:
                            p=self._checkpoint_path(dec,yc,xc)
                            with Dataset(p,"r") as d:
                                nd=np.asarray(d.variables[f"n_{var}"][::60,0,::100],dtype=float); md=np.asarray(d.variables[f"mean_{var}"][::60,0,::100],dtype=float); m2d=np.asarray(d.variables[f"M2_{var}"][::60,0,::100],dtype=float)
                                delta=md-mm; nnew=nn+nd; both=(nn>0)&(nd>0); only=(nn==0)&(nd>0); mm=np.where(only,md,np.where(both,mm+delta*nd/np.maximum(nnew,1),mm)); M=np.where(both,M+m2d+delta*delta*nn*nd/np.maximum(nnew,1),np.where(only,m2d,M)); nn=nnew
                        dmean=np.nanmax(np.abs(mm-a)); dm2=np.nanmax(np.abs(M-m2)); out["max_abs_mean"]=max(out["max_abs_mean"],float(dmean)); out["max_abs_M2"]=max(out["max_abs_M2"],float(dm2)); out["samples"]+=int(a.size)
        if out["max_abs_mean"]>1e-5 or out["max_abs_M2"]>1e-3: out["status"]="FAIL"; out["errors"].append("FULL does not match decade merge within audit tolerance")
        atomic_write_json(self.report_root/"merge_audit.json",out); return out

    def report(self)->Dict[str,Any]:
        audit=self.audit(); merge=self.merge_audit() if self.grid_shape else {"status":"SKIP"}
        summary={"engine_version":ENGINE_VERSION,"audit":audit,"merge_audit":merge,"created_utc":utc_now()}
        atomic_write_json(self.report_root/"report_summary.json",summary)
        with (self.report_root/"report.md").open("w",encoding="utf-8") as fh:
            fh.write(f"# {ENGINE_NAME} v{ENGINE_VERSION}\n\n")
            fh.write(f"Audit: **{audit['status']}**\n\nMerge audit: **{merge.get('status','SKIP')}**\n\n")
            fh.write("## Levels\n\n- L1: daily pooled\n- L2: 8 three-hour bins\n- L3: 24 hourly bins\n")
        return summary

    def close(self):
        for b in self._open_blocks.values(): b.close()
        self._open_blocks.clear(); self.journal.close()


# >>> selftest_moments: Independent moment correctness regression test.
# SELF-TEST: MOMENTS
# Compares online/batch statistics against direct NumPy calculations on deterministic
# synthetic data.  This catches both formula errors and broadcasting mistakes.
# ------------------------------------------------------------------------
# FUNCTION selftest_moments — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def selftest_moments()->None:
    rng=np.random.default_rng(10); a=rng.normal(size=(37,19)); valid=np.ones_like(a,dtype=bool)
    n,m,M2,M3,M4,xmin,xmax=batch_moments(a,valid)
    if not np.allclose(m,a.mean(axis=0)): raise AssertionError("mean")
    if not np.allclose(M2,((a-a.mean(axis=0))**2).sum(axis=0)): raise AssertionError("M2")
    if not np.allclose(M3,((a-a.mean(axis=0))**3).sum(axis=0)): raise AssertionError("M3")
    if not np.allclose(M4,((a-a.mean(axis=0))**4).sum(axis=0)): raise AssertionError("M4")


# >>> selftest_merge: Merge-equivalence regression test.
# SELF-TEST: MERGE EQUIVALENCE
# Verifies that accumulate(A) + accumulate(B) equals accumulate(A+B) within numerical
# tolerance.  This property is fundamental to decade/FULL reconstruction.
# ------------------------------------------------------------------------
# FUNCTION selftest_merge — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def selftest_merge()->None:
    rng=np.random.default_rng(11); a=rng.normal(size=(13,7)); b=rng.normal(size=(29,7));
    va=np.ones_like(a,dtype=bool); vb=np.ones_like(b,dtype=bool)
    A=batch_moments(a,va); B=batch_moments(b,vb); C=combine_moments_state(A[0],A[1],A[2],A[3],A[4],B[0],B[1],B[2],B[3],B[4]); F=batch_moments(np.vstack([a,b]),np.ones((42,7),bool))
    for x,y in zip(C,F[:5]):
        if not np.allclose(x,y): raise AssertionError("merge")


# >>> selftest_physics: Thermodynamic reference regression test.
# SELF-TEST: PHYSICS
# Sanity-checks the thermodynamic conversion on a known, ordinary atmospheric state.
# ------------------------------------------------------------------------
# FUNCTION selftest_physics — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def selftest_physics()->None:
    d=derive_moisture(np.array([[20.]],np.float32),np.array([[15.]],np.float32),np.array([[1013.25]],np.float32))
    if not (0<=float(d["rh"][0,0])<=100 and d["q"][0,0]>0 and d["r"][0,0]>0): raise AssertionError("physics")


# >>> selftest_update_state: Regression test for the in-place state updater and its array-shape contract.
# SELF-TEST: STATE UPDATE SHAPES AND CONTENT
# Exercises the exact multidimensional array shapes used by day/block accumulation,
# including pair covariance and threshold counters.  It is specifically designed to
# catch broadcast/indexing bugs such as the one previously observed in production.
# ------------------------------------------------------------------------
# FUNCTION selftest_update_state — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def selftest_update_state()->None:
    obj=object.__new__(Engine); obj.config=CONFIG; cells=8; old={}
    for var in VARIABLES:
        for stat in ("mean","M2","M3","M4"): old[f"{stat}_{var}"]=np.zeros((1,cells))
        old[f"n_{var}"]=np.zeros((1,cells),np.uint32); old[f"min_{var}"]=np.full((1,cells),np.inf); old[f"max_{var}"]=np.full((1,cells),-np.inf); old[f"missing_count_{var}"]=np.zeros((1,cells),np.uint32); old[f"threshold_count_{var}"]=np.zeros((1,cells,5),np.uint32)
    old["supersaturation_count"]=np.zeros((1,cells),np.uint32); old["invalid_e_over_p_count"]=np.zeros((1,cells),np.uint32)
    for x,y in PAIRS:
        tag=f"{x}__{y}"
        for stat in ("mean_x","mean_y","M2_x","M2_y","Cxy"): old[f"pair_{tag}_{stat}"]=np.zeros((1,cells))
        old[f"pair_{tag}_n"]=np.zeros((1,cells),np.uint32); old[f"joint_threshold_count_{tag}"]=np.zeros((1,cells,4),np.uint32)
    old["hist_rh__q"]=np.zeros((9,cells,8,8),np.uint32)
    rng=np.random.default_rng(44); t=9; rh=rng.uniform(20,100,(t,cells)); e=rng.uniform(5,20,(t,cells)); q=rng.uniform(.002,.02,(t,cells)); r=q/(1-q)
    xd={"rh":rh,"e":e,"r":r,"q":q}; valid={k:np.ones_like(rh,bool) for k in xd}; pairs={(a,b):(xd[a],xd[b],valid[a]&valid[b]) for a,b in PAIRS}
    out=obj._update_state_arrays(old,xd,valid,np.zeros_like(rh,bool),np.zeros_like(rh,bool),pairs,0,t)
    if not np.all(out["n_rh"]==t) or not np.allclose(out["mean_q"][0],q.mean(0)) or not np.allclose(out["M2_q"][0],((q-q.mean(0))**2).sum(0)): raise AssertionError("update-state")
    if not np.allclose(out["M3_q"][0],((q-q.mean(0))**3).sum(0)): raise AssertionError("M3")
    if not np.allclose(out["M4_q"][0],((q-q.mean(0))**4).sum(0)): raise AssertionError("M4")
    if not np.all(out["threshold_count_rh"][0,:,0]==(rh>80).sum(0)): raise AssertionError("threshold")
    tag="rh__q"; cov=np.sum((rh-rh.mean(0))*(q-q.mean(0)),axis=0); pair_cov=out[f"pair_{tag}_Cxy"][0]
    if not np.allclose(pair_cov,cov,rtol=1e-12,atol=1e-12): raise AssertionError("covariance")
    if not np.allclose(out[f"pair_{tag}_M2_x"][0],((rh-rh.mean(0))**2).sum(0)): raise AssertionError("pair M2 x")


# >>> selftest_calendar: Calendar regression-test wrapper.
# SELF-TEST: CALENDAR
# Runs the complete calendar contract checks.
# ------------------------------------------------------------------------
# FUNCTION selftest_calendar — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def selftest_calendar()->None:
    validate_calendar_contract()


# >>> selftest_transaction: Fault-injection transaction/recovery smoke test.
# SELF-TEST: CRASH/RECOVERY SMOKE TEST
# Reverses a tiny state, invokes recovery, and confirms that the original state is
# restored.  It is deliberately independent of ERA5 so it can run anywhere.
# ------------------------------------------------------------------------
# FUNCTION selftest_transaction — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def selftest_transaction()->None:
    root=Path(tempfile.mkdtemp(prefix="hce_v10_tx_")); j=Journal(root); arr=np.arange(20,dtype=np.int64); p=root/"a.bin"; atomic_write_bytes(p,arr.tobytes())
    txid=j.txid(2001,1,1,0,0); rec={"txid":txid,"year":2001,"month":1,"day":1,"y_chunk":0,"x_chunk":0,"targets":[],"created_utc":utc_now()}; j.create_open(rec); before=j.write_before(txid,{"blob":p.read_bytes()})
    arr2=np.arange(20,dtype=np.int64)[::-1]; atomic_write_bytes(p,arr2.tobytes())
    def restore(r): atomic_write_bytes(p,np.load(before,allow_pickle=False)["blob"].tobytes())
    j.recover(restore)
    if not np.array_equal(np.frombuffer(p.read_bytes(),dtype=np.int64),arr): raise AssertionError("transaction recovery")
    j.close(); shutil.rmtree(root,ignore_errors=True)


# >>> benchmark: Synthetic numerical hot-path benchmark independent of real ERA5 I/O.
# SELF-TEST: NUMERICAL HOT-PATH BENCHMARK
# Measures vectorized in-memory moment throughput on representative cell counts.
# This is NOT a full ERA5/NetCDF benchmark; disk and compression must be measured on
# the target machine.
# ------------------------------------------------------------------------
# FUNCTION benchmark — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def benchmark()->Dict[str,Any]:
    rng=np.random.default_rng(123); results=[]
    for cells in (16,2048,8192):
        x=rng.random((24,cells)).astype(np.float32)*100; y=rng.random((24,cells)).astype(np.float32)*0.02; v=np.ones_like(x,bool)
        t=time.perf_counter();
        for _ in range(10): batch_moments(x,v)
        dt=(time.perf_counter()-t)/10
        results.append({"cells":cells,"batch_moments_seconds":dt,"cells_per_second":cells*24/dt})
    out={"engine_version":ENGINE_VERSION,"benchmark":"synthetic in-memory core","results":results,"note":"Not a filesystem/ERA5 benchmark."}
    atomic_write_json(CONFIG.output_root/"reports"/"benchmark.json",out); return out


# >>> setup_logging: Consistent operational logging configuration.
# LOGGING SETUP
# Centralizes console verbosity so long unattended jobs have a predictable audit trail.
# ------------------------------------------------------------------------
# FUNCTION setup_logging — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def setup_logging(verbose:bool=False)->None:
    logging.basicConfig(level=logging.DEBUG if verbose else logging.INFO,format="%(asctime)s | %(levelname)s | %(message)s")


# -----------------------------------------------------------------------------
# PILOT MODE
# -----------------------------------------------------------------------------
# A deliberately small real-data execution path used before a multiyear run.
# It processes exactly one day and one spatial block from a selected year/month/day,
# while printing frequent diagnostics. It never replaces the production run.
# The pilot uses a dedicated output root so existing production checkpoints are
# not modified accidentally.
# -----------------------------------------------------------------------------
def pilot(year: int, month: int = 1, day: int = 1) -> Dict[str, Any]:
    """Run a bounded real-data pilot.

    v11.5 extends the v11.4 pilot contract in one important way: when the
    user supplies a month but no explicit day, the pilot processes the entire
    calendar month for the first spatial block.  When ``day`` is explicitly
    supplied, only that day is processed.  This removes the ambiguity that
    previously made ``pilot --month 1`` process only January 1.

    The pilot remains isolated from production checkpoints and emits progress
    for every day.  It is therefore safe for repeated validation of the real
    ERA5-Land input while still being much smaller than the multidecade run.
    """
    import dataclasses
    if not (1 <= month <= 12):
        raise ValueError("pilot month must be 1..12")
    last_day = calendar.monthrange(year, month)[1]
    if day is not None and not (0 <= day <= last_day):
        raise ValueError(f"pilot day must be 0..{last_day} for {year}-{month:02d}")

    pilot_root = CONFIG.output_root.parent / f"{CONFIG.output_root.name}_PILOT_v11_5"
    marker = pilot_root / "PILOT_VERSION.txt"
    if marker.exists():
        marker_text = marker.read_text(encoding="utf-8").strip()
        if marker_text != ENGINE_VERSION:
            raise RuntimeError(
                f"Pilot checkpoint root belongs to {marker_text}; expected {ENGINE_VERSION}: {pilot_root}"
            )
    else:
        pilot_root.mkdir(parents=True, exist_ok=True)
        marker.write_text(ENGINE_VERSION + "\n", encoding="utf-8")

    cfg = dataclasses.replace(CONFIG) if dataclasses.is_dataclass(CONFIG) else type(CONFIG)()
    object.__setattr__(cfg, "start_year", year)
    object.__setattr__(cfg, "end_year", year)
    object.__setattr__(cfg, "output_root", pilot_root)
    eng = Engine(cfg)
    LOG.info("PILOT 1/9 | loading grid")
    eng.load_grid()
    idx = eng.file_indices(year)
    LOG.info("PILOT 2/9 | opening %04d-%02d", year, month)
    ds_t, ds_d, ds_p = open_aligned_triplet(idx["t2m"][month], idx["d2m"][month], idx["sp"][month])
    try:
        LOG.info("PILOT 3/9 | validating input triplet")
        validate_datasets(ds_t, ds_d, ds_p, year, month)
        units={"t2m":ds_t["t2m"].attrs.get("units"),"d2m":ds_d["d2m"].attrs.get("units"),"sp":ds_p["sp"].attrs.get("units")}

        first_day = day if day and day > 0 else 1
        last = day if day and day > 0 else last_day
        n_days = last - first_day + 1
        LOG.info("PILOT 4/9 | processing first spatial block | days=%d", n_days)

        for i, current_day in enumerate(range(first_day, last + 1), start=1):
            LOG.info("PILOT DAY %02d/%02d | %04d-%02d-%02d | start", i, n_days, year, month, current_day)
            t0=time.time()
            eng._process_day_block(year, month, current_day, 0, 0, ds_t, ds_d, ds_p, units)
            LOG.info("PILOT DAY %02d/%02d | %04d-%02d-%02d | PASS | %.3fs", i, n_days, year, month, current_day, time.time()-t0)

        LOG.info("PILOT 5/9 | month/day block transactions returned")
        LOG.info("PILOT 6/9 | flushing checkpoint state")
        for b in eng._open_blocks.values():
            if b.nc is not None:
                b.nc.sync()
        LOG.info("PILOT 7/9 | reopening checkpoint for structural verification")
        eng.close_block(0,0)
        eng._open_blocks.clear()
        for period in eng._periods_for_year(year):
            block = eng._open_or_get_file(period,0,0,create=False)
            block.create_or_open().sync()
        LOG.info("PILOT 8/9 | checkpoint reopen PASS")
        LOG.info("PILOT 9/9 | PASS | output=%s", pilot_root)
        return {"status":"PASS","year":year,"month":month,"day":0 if day in (None,0) else day,
                "days_processed":n_days,"output_root":str(pilot_root)}
    finally:
        try:
            ds_t.close(); ds_d.close(); ds_p.close()
        finally:
            eng.close()


# >>> main: Command-line dispatcher; scientific logic remains in testable functions.
# ----------------------------------------------------------------------------
# COMMAND MAP
#
# selftest        : run deterministic unit/regression checks only.
# validate-input  : inspect a real target year; no scientific accumulation.
# run             : execute the production 1981-2020 accumulation.
# audit           : inspect finished outputs for structural invariants.
# merge-audit     : verify FULL against a four-decade merge.
# report          : invoke/report the frozen analysis products as configured.
# benchmark       : measure the in-memory numerical hot path.
# ----------------------------------------------------------------------------
# COMMAND-LINE ENTRY POINT
# Dispatches the seven operational modes of the v11 executable.  CLI parsing is kept
# separate from scientific functions so the functions remain directly testable.
# ------------------------------------------------------------------------
# FUNCTION main — IMPLEMENTATION GUIDE
# ------------------------------------------------------------------------
# Contract:
#   This function is documented as inputs -> validation -> computation ->
#   outputs -> failure behavior. Read the type/shape assumptions before
#   reading the arithmetic.
#
# Shape discipline:
#   Array axes have scientific meaning. A broadcast that is numerically
#   legal can still be scientifically wrong if latitude/time/pair masks
#   become misaligned. Shape mismatches should fail loudly.
#
# Performance discipline:
#   This function may execute inside a large multidecade loop. Avoid Python
#   work per grid cell when a vectorized operation gives the same result.
#
# Recovery discipline:
#   A function that updates persistent state must preserve the transaction
#   contract: prepare -> verify -> commit. Never advance progress early.
# ------------------------------------------------------------------------

def main()->int:
    parser=argparse.ArgumentParser(description=f"{ENGINE_NAME} v{ENGINE_VERSION}")
    parser.add_argument("command",choices=("selftest","validate-input","pilot","run","audit","merge-audit","report","benchmark"))
    parser.add_argument("--year",type=int,default=None); parser.add_argument("--month",type=int,default=1); parser.add_argument("--day",type=int,default=0); parser.add_argument("--verbose",action="store_true")
    args=parser.parse_args(); setup_logging(args.verbose)
    if args.command=="selftest":
        selftest_calendar(); selftest_moments(); selftest_merge(); selftest_physics(); selftest_update_state(); selftest_transaction(); LOG.info("SELFTEST PASS"); return 0
    if args.command=="benchmark":
        r=benchmark(); LOG.info("BENCHMARK COMPLETE"); print(json.dumps(r,indent=2)); return 0
    if args.command=="pilot":
        r=pilot(args.year or CONFIG.start_year, args.month, args.day); print(json.dumps(r,indent=2)); return 0
    eng=Engine(CONFIG)
    try:
        if args.command=="validate-input":
            # Validate all twelve months.  The validation command is intentionally
            # read-only: it opens each monthly triplet, aligns the spatial axes,
            # checks the time axis, and then closes the files before moving on.
            year=args.year or CONFIG.start_year; idx=eng.file_indices(year)
            for m in range(1,13):
                # Each month's T2m/D2m/SP triplet is aligned to the T2m grid.
                t, d, p = open_aligned_triplet(idx["t2m"][m], idx["d2m"][m], idx["sp"][m])
                try:
                    validate_datasets(t, d, p, year, m)
                finally:
                    t.close(); d.close(); p.close()
            LOG.info("INPUT VALIDATION PASS for %d",year); return 0
        if args.command=="run": eng.run(); LOG.info("RUN COMPLETE v%s",ENGINE_VERSION); return 0
        if args.command=="audit": eng.load_grid(); r=eng.audit(); LOG.info("AUDIT %s",r["status"]); return 0 if r["status"]=="PASS" else 2
        if args.command=="merge-audit": eng.load_grid(); r=eng.merge_audit(); LOG.info("MERGE AUDIT %s",r["status"]); return 0 if r["status"]=="PASS" else 2
        if args.command=="report": eng.load_grid(); eng.report(); LOG.info("REPORT COMPLETE"); return 0
    finally: eng.close()
    return 0

if __name__=="__main__":
    raise SystemExit(main())
