#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read-only preflight for moisture_climatology_v8_1_FAST.

This script does NOT process ERA5 data. It checks configuration compatibility,
existing annual checkpoints, and reports resumable progress.
"""
from pathlib import Path
import json
import sys
import numpy as np

import moisture_climatology_v8_1_FAST as m


def safe_flag_sum(path: Path) -> tuple[int, int]:
    from netCDF4 import Dataset
    with Dataset(path, "r") as ds:
        v = ds.variables["completed_chunk"][:]
        done = int(np.ma.filled(v, 0).sum())
        total = int(np.prod(v.shape))
    return done, total


def main() -> int:
    print("=" * 78)
    print("v8.1 FAST PREFLIGHT (READ ONLY)")
    print("=" * 78)
    print(f"OUTPUT_DIR       : {m.OUTPUT_DIR}")
    print(f"CHECKPOINT_DIR   : {m.CHECKPOINT_DIR}")
    print(f"YEAR_DIR         : {m.YEAR_DIR}")
    print(f"CONFIG_HASH      : {m.CONFIG_HASH}")
    print(f"CHECKPOINT_VER   : {m.CHECKPOINT_VERSION}")
    print(f"CHUNK            : {m.CHUNK_LAT} x {m.CHUNK_LON}")
    print(f"WORKERS          : {m.MAX_WORKERS}")
    print(f"EMPIRICAL_BIVAR  : {m.BUILD_EMPIRICAL_BIVARIATE}")
    print()

    if not m.YEAR_DIR.exists():
        print("NO_CHECKPOINT_DIR")
        return 0

    ok = 0
    running = 0
    missing = 0
    for year in range(m.START_YEAR, m.END_YEAR + 1):
        final_path, json_path, part_path = m.year_paths(year)
        if m.is_year_complete(year):
            print(f"{year}: COMPLETE")
            ok += 1
            continue
        if part_path.exists():
            try:
                done, total = safe_flag_sum(part_path)
                meta = m._read_progress(json_path)
                print(f"{year}: RESUME {done}/{total} units ({100.0*done/max(total,1):.2f}%) | part={part_path.name}")
                if meta:
                    print(f"      JSON: completed={meta.get('completed_units')} total={meta.get('total_units')} hash={meta.get('config_hash')}")
                running += 1
            except Exception as exc:
                print(f"{year}: CHECKPOINT_READ_ERROR: {exc}")
        else:
            print(f"{year}: NOT_STARTED")
            missing += 1

    print()
    print(f"SUMMARY: complete={ok} resume={running} not_started={missing}")
    print("NO ERA5 PROCESSING WAS PERFORMED.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
