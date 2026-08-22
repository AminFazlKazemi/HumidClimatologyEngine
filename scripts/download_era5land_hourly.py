#!/usr/bin/env python3
"""Minimal hourly ERA5-Land fallback downloader."""
from __future__ import annotations
import argparse
from pathlib import Path
import cdsapi

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", required=True)
    ap.add_argument("--month", required=True)
    ap.add_argument("--variables", nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--area", nargs=4, type=float)
    args = ap.parse_args()

    req = {
        "variable": args.variables,
        "year": args.year,
        "month": args.month,
        "day": [f"{d:02d}" for d in range(1,32)],
        "time": [f"{h:02d}:00" for h in range(24)],
    }
    if args.area:
        req["area"] = args.area
    args.out.parent.mkdir(parents=True, exist_ok=True)
    cdsapi.Client().retrieve("reanalysis-era5-land", req, str(args.out))

if __name__ == "__main__":
    main()
