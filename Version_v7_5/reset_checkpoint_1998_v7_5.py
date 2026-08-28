#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reset ONLY the 1998 checkpoint artifacts for HumidClimatologyEngine v7.5.

SAFE BY DEFAULT:
    python reset_checkpoint_1998_v7_5.py

ACTUAL RESET:
    python reset_checkpoint_1998_v7_5.py --yes

This script does NOT delete:
- any other year's checkpoint
- ERA5-Land input files
- final production outputs
- diagnostic files
- source code / README
"""

from pathlib import Path
import argparse
import shutil
import sys

BASE = Path(r"C:\c")
CHECKPOINT_DIR = BASE / "checkpoints_moisture_v7_5"
YEAR_DIR = CHECKPOINT_DIR / "years"

YEAR = 1998


def human_size(n: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    value = float(n)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{n} B"


def file_size(path: Path) -> int:
    try:
        return path.stat().st_size if path.is_file() else 0
    except OSError:
        return 0


def collect_targets():
    """
    Collect ONLY checkpoint files belonging to 1998.

    Expected examples:
        year_1998_083c99776caa07176970.part.nc
        year_1998_*.part.nc

    Nothing outside the 1998 year namespace is touched.
    """
    if not YEAR_DIR.exists():
        return []

    patterns = [
        f"year_{YEAR}_*.part.nc",
        f"year_{YEAR}_*.nc",
    ]

    found = {}
    for pattern in patterns:
        for p in YEAR_DIR.glob(pattern):
            if p.is_file():
                found[p.resolve()] = p

    return sorted(found.values(), key=lambda p: p.name.lower())


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reset ONLY the 1998 HumidClimatologyEngine v7.5 "
            "checkpoint files."
        )
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Actually delete the 1998 checkpoint files.",
    )
    args = parser.parse_args()

    print("=" * 100)
    print("HumidClimatologyEngine v7.5 — 1998 CHECKPOINT RESET")
    print("=" * 100)
    print(f"Base       : {BASE}")
    print(f"Checkpoint : {CHECKPOINT_DIR}")
    print(f"Year dir   : {YEAR_DIR}")
    print(f"Target year: {YEAR}")
    print()

    if not CHECKPOINT_DIR.exists():
        print("Checkpoint directory does not exist.")
        print("Nothing to reset.")
        return 0

    targets = collect_targets()

    print("SAFETY SCOPE")
    print("-" * 100)
    print(f"ONLY files matching year_{YEAR}_*.part.nc / year_{YEAR}_*.nc")
    print("will be considered.")
    print("Other years will NOT be deleted.")
    print()

    if not targets:
        print(f"No {YEAR} checkpoint files found.")
        print("Nothing to reset.")
        return 0

    reclaimable = sum(file_size(p) for p in targets)

    print(f"Found {len(targets)} checkpoint file(s):")
    for p in targets:
        print(f"  [DELETE] {p} ({human_size(file_size(p))})")

    print()
    print(f"Total size to remove: {human_size(reclaimable)}")

    if not args.yes:
        print()
        print("=" * 100)
        print("DRY RUN — NO FILES WERE DELETED")
        print("=" * 100)
        print()
        print("If the list above is correct, run:")
        print()
        print("  python reset_checkpoint_1998_v7_5.py --yes")
        print()
        return 0

    print()
    print("=" * 100)
    print("DELETING 1998 CHECKPOINTS")
    print("=" * 100)

    deleted = 0
    failed = 0

    for p in targets:
        try:
            p.unlink()
            deleted += 1
            print(f"[DELETED] {p}")
        except Exception as exc:
            failed += 1
            print(f"[FAILED ] {p}")
            print(f"          {type(exc).__name__}: {exc}")

    print()
    print("=" * 100)

    if failed:
        print("RESET FINISHED WITH ERRORS")
        print(f"Deleted : {deleted}")
        print(f"Failed  : {failed}")
        return 1

    print("RESET COMPLETE — ONLY 1998 WAS RESET")
    print(f"Deleted : {deleted} file(s)")
    print(f"Freed  : {human_size(reclaimable)}")
    print()
    print("All other year checkpoints remain untouched.")
    print("ERA5-Land input files remain untouched.")
    print("Final outputs remain untouched.")
    print("=" * 100)

    return 0


if __name__ == "__main__":
    sys.exit(main())
