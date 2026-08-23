from __future__ import annotations

import argparse
import shutil
from pathlib import Path

DEFAULT_BASE = Path(r"C:\c")
RESET_TOKEN = "RESET V7.5"

# Only v7.5 generated artifacts are allowed to be removed.
# v6 artifacts are intentionally excluded.
TARGETS = (
    "checkpoints_moisture_v7_5",
    "moisture_climatology_1981_2020_v7_5.nc",
    "moisture_climatology_diagnostics_1981_2020_v7_5.nc",
    "moisture_climatology_bivariate_1981_2020_v7_5.nc",
    "moisture_climatology_run_manifest_v7_5.json",
    "moisture_bivariate_empirical_rh__q_1981_2020_v7_5.nc",
)


def _safe_size(path: Path) -> int:
    """Return an approximate file size without following directory symlinks."""
    try:
        if path.is_symlink():
            return path.stat().st_size
        if path.is_file():
            return path.stat().st_size
        if path.is_dir():
            total = 0
            for p in path.rglob("*"):
                try:
                    if p.is_symlink():
                        continue
                    if p.is_file():
                        total += p.stat().st_size
                except OSError:
                    continue
            return total
    except OSError:
        return 0
    return 0


def resolve_targets(base: Path) -> list[Path]:
    return [base / name for name in TARGETS]


def show_plan(base: Path, targets: list[Path]) -> int:
    present = 0
    total_bytes = 0

    print("\nHumidClimatologyEngine v7.5 RESET — DRY RUN")
    print(f"Base: {base}")
    print("Targets:")

    for p in targets:
        try:
            exists = p.exists() or p.is_symlink()
            if not exists:
                print(f"  [MISS] {p}")
                continue

            present += 1
            size = _safe_size(p)
            total_bytes += size
            kind = "SYMLINK" if p.is_symlink() else "DIR" if p.is_dir() else "FILE"
            print(f"  [{kind:7}] {p}  ({size / 1024**2:.2f} MiB)")
        except OSError as exc:
            print(f"  [ERROR ] {p}  ({exc})")

    print(f"Present targets: {present}/{len(targets)}")
    print(f"Approximate reclaimable size: {total_bytes / 1024**2:.2f} MiB")
    print("No changes made.")
    print("\nTo delete explicitly, run with: --yes")
    return present


def delete_target(path: Path) -> None:
    if path.is_symlink():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Safely reset HumidClimatologyEngine v7.5 generated artifacts."
    )
    parser.add_argument(
        "--base",
        type=Path,
        default=DEFAULT_BASE,
        help=r"Working directory (default: C:\c)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show exactly what would be deleted; this is also the default mode.",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Perform deletion without an interactive prompt.",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Ask for the exact confirmation token before deletion.",
    )
    args = parser.parse_args()

    base = args.base.expanduser().resolve()
    targets = resolve_targets(base)

    # Safety invariant: unless --yes or --interactive is explicitly supplied,
    # this program is always a non-destructive dry-run.
    if not args.yes and not args.interactive:
        show_plan(base, targets)
        return 0

    present = show_plan(base, targets)
    if present == 0:
        print("Nothing to delete.")
        return 0

    if args.interactive:
        try:
            answer = input(f'Type "{RESET_TOKEN}" to continue: ').strip()
        except (EOFError, KeyboardInterrupt):
            print("Reset cancelled: no interactive input available.")
            return 2
        if answer != RESET_TOKEN:
            print("Reset cancelled.")
            return 2

    deleted = 0
    failures = []
    for path in targets:
        if not (path.exists() or path.is_symlink()):
            continue
        try:
            delete_target(path)
            deleted += 1
            print(f"Deleted: {path}")
        except Exception as exc:
            failures.append((path, exc))
            print(f"ERROR deleting {path}: {exc}")

    print(f"\nDeleted: {deleted}/{present} target(s).")
    if failures:
        print("Reset finished with errors:")
        for path, exc in failures:
            print(f"  - {path}: {exc}")
        return 1

    print("SUCCESS: HumidClimatologyEngine v7.5 reset complete.")
    print("v6 artifacts were not targeted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
