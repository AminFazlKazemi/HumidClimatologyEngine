from __future__ import annotations

import argparse
import shutil
from pathlib import Path

DEFAULT_BASE = Path(r"C:\c")

# v7.5 artifacts only. No v6 artifacts are targeted.
TARGETS = [
    "checkpoints_moisture_v7_5",
    "moisture_climatology_1981_2020_v7_5.nc",
    "moisture_climatology_diagnostics_1981_2020_v7_5.nc",
    "moisture_climatology_bivariate_1981_2020_v7_5.nc",
    "moisture_climatology_run_manifest_v7_5.json",
    "moisture_bivariate_empirical_rh__q_1981_2020_v7_5.nc",
]


def resolve_targets(base: Path) -> list[Path]:
    return [base / name for name in TARGETS]


def print_plan(base: Path, targets: list[Path]) -> tuple[int, int]:
    present = 0
    total_bytes = 0
    print("\nHumidClimatologyEngine v7.5 RESET — DRY RUN")
    print(f"Base: {base}")
    print("Targets:")
    for p in targets:
        if p.is_dir():
            present += 1
            size = sum(f.stat().st_size for f in p.rglob("*") if f.is_file())
            total_bytes += size
            print(f"  [DIR ] {p}  ({size / 1024**2:.2f} MiB)")
        elif p.exists():
            present += 1
            size = p.stat().st_size
            total_bytes += size
            print(f"  [FILE] {p}  ({size / 1024**2:.2f} MiB)")
        else:
            print(f"  [MISS] {p}")
    print(f"Present targets: {present}/{len(targets)}")
    print(f"Approximate reclaimable size: {total_bytes / 1024**2:.2f} MiB")
    print("No changes made.")
    return present, total_bytes


def delete_target(p: Path) -> None:
    if p.is_symlink():
        p.unlink()
    elif p.is_dir():
        shutil.rmtree(p)
    elif p.exists():
        p.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Safely reset HumidClimatologyEngine v7.5 generated artifacts."
    )
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE,
                        help=r"Working directory (default: C:\c)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show exactly what would be deleted without deleting anything.")
    parser.add_argument("--yes", action="store_true",
                        help="Delete without interactive confirmation.")
    args = parser.parse_args()

    base = args.base.expanduser().resolve()
    targets = resolve_targets(base)
    present, _ = print_plan(base, targets)

    if args.dry_run or present == 0:
        return 0

    if not args.yes:
        answer = input("Type RESET V7.5 to continue: ").strip()
        if answer != "RESET V7.5":
            print("Reset cancelled.")
            return 2

    deleted = 0
    for p in targets:
        if p.exists() or p.is_symlink():
            delete_target(p)
            deleted += 1
            print(f"Deleted: {p}")

    print(f"\nSUCCESS: v7.5 reset complete. Deleted {deleted} target(s).")
    print("All v7.5 checkpoints, progress state, outputs, and run manifest were targeted.")
    print("v6 artifacts were not targeted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
