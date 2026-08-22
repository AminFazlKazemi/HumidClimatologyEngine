from pathlib import Path
import shutil

CHECKPOINT_DIR = Path(r"C:\c\checkpoints_moisture_v6")

if CHECKPOINT_DIR.exists():
    shutil.rmtree(CHECKPOINT_DIR)
    print(f"✅ Checkpoint deleted: {CHECKPOINT_DIR}")
else:
    print(f"ℹ️ Checkpoint directory not found: {CHECKPOINT_DIR}")
from pathlib import Path
import shutil

BASE = Path(r"C:\c")

for name in [
    "checkpoints_moisture_v6",
    "moisture_climatology_1981_2020.nc",
    "moisture_climatology_diagnostics_1981_2020.nc",
]:
    p = BASE / name
    if p.is_dir():
        shutil.rmtree(p)
        print(f"✅ Deleted directory: {p}")
    elif p.exists():
        p.unlink()
        print(f"✅ Deleted file: {p}")
    else:
        print(f"ℹ️ Not found: {p}")