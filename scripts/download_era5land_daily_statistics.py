#!/usr/bin/env python3
"""Download official ERA5-Land post-processed daily statistics from CDS."""
from __future__ import annotations
import argparse, calendar, logging
from pathlib import Path
import cdsapi

DATASET = "derived-era5-land-daily-statistics"
DEFAULT_VARIABLES = ["2m_temperature","2m_dewpoint_temperature","surface_pressure"]
DEFAULT_STATISTICS = ["daily_mean","daily_minimum","daily_maximum"]

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger("ERA5LandDailyDownloader")

def parse_years(spec):
    if "-" in spec:
        a,b = map(int, spec.split("-",1))
        return list(range(a,b+1))
    return [int(x) for x in spec.split(",") if x.strip()]

def month_days(year, month):
    return [f"{d:02d}" for d in range(1, calendar.monthrange(year, month)[1]+1)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", required=True)
    ap.add_argument("--variables", nargs="+", default=DEFAULT_VARIABLES)
    ap.add_argument("--statistics", nargs="+", default=DEFAULT_STATISTICS,
                    choices=["daily_mean","daily_minimum","daily_maximum"])
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--time-zone", default="utc+00:00")
    ap.add_argument("--frequency", default="1_hourly", choices=["1_hourly","3_hourly","6_hourly"])
    ap.add_argument("--area", nargs=4, type=float, metavar=("N","W","S","E"))
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    client = cdsapi.Client()
    years = parse_years(args.years)
    for year in years:
        for month in range(1,13):
            for stat in args.statistics:
                var_tag = "_".join(v.replace(" ","-") for v in args.variables)
                target = args.out_dir / var_tag / stat / f"{year:04d}" / f"era5land_{var_tag}_{stat}_{year:04d}{month:02d}.zip"
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists() and not args.overwrite:
                    log.info("SKIP existing %s", target)
                    continue
                req = {
                    "variable": args.variables,
                    "year": f"{year:04d}",
                    "month": f"{month:02d}",
                    "day": month_days(year, month),
                    "daily_statistic": stat,
                    "time_zone": args.time_zone,
                    "frequency": args.frequency,
                }
                if args.area:
                    req["area"] = args.area
                log.info("DOWNLOAD | %04d-%02d | %s", year, month, stat)
                client.retrieve(DATASET, req).download(str(target))
                log.info("DONE | %s", target)

if __name__ == "__main__":
    main()
