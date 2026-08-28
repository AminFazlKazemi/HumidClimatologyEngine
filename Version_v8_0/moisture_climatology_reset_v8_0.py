#!/usr/bin/env python
"""
MOISTURE CLIMATOLOGY v8.0 FINAL
Checkpoint reset utility.

Removes v8 production checkpoints safely before a clean full rebuild.
"""

from pathlib import Path
import shutil
import argparse

def reset_checkpoints(path):
    p = Path(path)
    removed = 0
    for item in p.rglob("*checkpoint*"):
        if item.is_file():
            item.unlink()
            removed += 1
        elif item.is_dir():
            shutil.rmtree(item)
            removed += 1
    print(f"Removed checkpoint objects: {removed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Project/output directory")
    args = parser.parse_args()
    reset_checkpoints(args.path)
