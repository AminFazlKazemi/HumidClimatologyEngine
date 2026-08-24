# diagnose_feb1998_hourly.py

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

# ============================================================
# FILES
# ============================================================

FILES = {
    "T2M": Path(
        r"F:\Kazemi\era5\land\T2m\T2m199802.nc"
    ),
    "D2M": Path(
        r"F:\Kazemi\era5\land\Dew_Point_Temperature\Dew_Point_Temperature199802.nc"
    ),
    "SP": Path(
        r"F:\Kazemi\era5\land\Surface_Pressure\Surface_Pressure199802.nc"
    ),
}

OUT = Path("diagnose_1998_02_hourly.csv")

# ============================================================
# SETTINGS
# ============================================================

# اگر تعداد finite یک ساعت نسبت به baseline بیشتر از این مقدار
# افت کند، آن ساعت مشکوک اعلام می‌شود.
DROP_THRESHOLD_PERCENT = 0.05

# ============================================================
# HEADER
# ============================================================

print("=" * 110)
print("FEBRUARY 1998 — HOURLY DATA FORENSICS")
print("=" * 110)
print("هدف: یافتن دقیق روز و ساعت خرابی/ناقص بودن داده")
print()


# ============================================================
# ANALYZE ONE FILE
# ============================================================

def analyze_file(name, path):

    print()
    print("=" * 110)
    print(f"{name}: {path}")
    print("=" * 110)

    if not path.exists():
        print("!!! FILE NOT FOUND !!!")
        return None

    print(f"Size: {path.stat().st_size / 1024**2:.2f} MiB")

    ds = xr.open_dataset(path)

    # --------------------------------------------------------
    # variable
    # --------------------------------------------------------

    candidates = {
        "T2M": "t2m",
        "D2M": "d2m",
        "SP": "sp",
    }

    var = candidates[name]

    if var not in ds:
        print(f"!!! VARIABLE {var} NOT FOUND !!!")
        ds.close()
        return None

    da = ds[var]

    print(f"Variable : {var}")
    print(f"Dims     : {da.dims}")
    print(f"Shape    : {da.shape}")
    print(f"Units    : {da.attrs.get('units', '')}")

    # --------------------------------------------------------
    # time
    # --------------------------------------------------------

    time = pd.DatetimeIndex(ds.time.values)

    print(f"First    : {time[0]}")
    print(f"Last     : {time[-1]}")
    print(f"Hours    : {len(time)}")

    # Expected February 1998
    expected = pd.date_range(
        "1998-02-01 00:00:00",
        "1998-02-28 23:00:00",
        freq="h",
    )

    missing_times = expected.difference(time)
    extra_times = time.difference(expected)

    if len(missing_times):
        print("!!! MISSING TIMESTAMPS !!!")
        for t in missing_times:
            print("   ", t)

    if len(extra_times):
        print("!!! EXTRA TIMESTAMPS !!!")
        for t in extra_times:
            print("   ", t)

    # --------------------------------------------------------
    # HOURLY SCAN
    # --------------------------------------------------------

    records = []

    print()
    print("Scanning every hour ...")

    for i, t in enumerate(time):

        # Load exactly one hourly 301x301 field
        field = da.isel(time=i).values

        total = field.size

        finite = np.isfinite(field).sum()
        nonfinite = total - finite

        finite_pct = 100.0 * finite / total

        if finite:
            finite_values = field[np.isfinite(field)]
            vmin = float(np.min(finite_values))
            vmax = float(np.max(finite_values))
        else:
            vmin = np.nan
            vmax = np.nan

        records.append({
            "file": name,
            "time": t,
            "date": t.date(),
            "hour": t.hour,
            "total": total,
            "finite": int(finite),
            "nonfinite": int(nonfinite),
            "finite_pct": finite_pct,
            "min": vmin,
            "max": vmax,
        })

    ds.close()

    df = pd.DataFrame(records)

    # --------------------------------------------------------
    # BASELINE
    # --------------------------------------------------------
    #
    # Important:
    # The original files have ~10.9458% non-finite values.
    # This is probably a persistent spatial mask.
    #
    # Therefore we do NOT call every NaN an error.
    #
    # We calculate the normal hourly finite-count pattern and
    # search for sudden deviations.
    # --------------------------------------------------------

    median_finite_pct = df["finite_pct"].median()

    df["drop_from_median_pct"] = (
        median_finite_pct - df["finite_pct"]
    )

    df["relative_drop_pct"] = (
        100.0
        * df["drop_from_median_pct"]
        / median_finite_pct
    )

    df["suspicious"] = (
        df["relative_drop_pct"]
        > DROP_THRESHOLD_PERCENT
    )

    # --------------------------------------------------------
    # REPORT
    # --------------------------------------------------------

    print()
    print(f"Normal median finite : {median_finite_pct:.8f}%")
    print(
        f"Minimum finite       : "
        f"{df['finite_pct'].min():.8f}%"
    )
    print(
        f"Maximum finite       : "
        f"{df['finite_pct'].max():.8f}%"
    )

    suspicious = df[df["suspicious"]].copy()

    print()
    print("-" * 110)

    if suspicious.empty:
        print("NO SUSPICIOUS HOURLY DROP FOUND.")
    else:

        print(
            f"!!! {len(suspicious)} SUSPICIOUS HOURS FOUND !!!"
        )

        print()

        for _, r in suspicious.iterrows():

            print(
                f"{r['time']} | "
                f"finite={r['finite']:,}/{r['total']:,} | "
                f"finite={r['finite_pct']:.6f}% | "
                f"drop={r['relative_drop_pct']:.6f}% | "
                f"min={r['min']:.6f} | "
                f"max={r['max']:.6f}"
            )

    return df


# ============================================================
# RUN ALL THREE
# ============================================================

all_results = {}

for name, path in FILES.items():
    result = analyze_file(name, path)

    if result is not None:
        all_results[name] = result


# ============================================================
# COMBINE
# ============================================================

if not all_results:
    raise RuntimeError("No files could be analyzed.")

combined = pd.concat(
    all_results.values(),
    ignore_index=True
)

combined.to_csv(
    OUT,
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# COMMON HOURLY TIMELINE
# ============================================================

print()
print()
print("=" * 110)
print("CROSS-FILE COMPARISON")
print("=" * 110)

pivot = combined.pivot_table(
    index="time",
    columns="file",
    values="finite_pct",
)

# ------------------------------------------------------------
# Difference from each file's median
# ------------------------------------------------------------

for name in FILES:
    if name in pivot.columns:
        med = pivot[name].median()

        pivot[f"{name}_drop_pct"] = (
            100.0
            * (med - pivot[name])
            / med
        )


# ============================================================
# COMMON ANOMALIES
# ============================================================

drop_columns = [
    c for c in pivot.columns
    if c.endswith("_drop_pct")
]

if drop_columns:

    # هر فایل باید افت قابل توجه داشته باشد
    threshold = DROP_THRESHOLD_PERCENT

    for name in FILES:
        col = f"{name}_drop_pct"

        if col in pivot:
            pivot[f"{name}_BAD"] = (
                pivot[col] > threshold
            )

    bad_cols = [
        f"{name}_BAD"
        for name in FILES
        if f"{name}_BAD" in pivot
    ]

    pivot["NUMBER_OF_BAD_FILES"] = (
        pivot[bad_cols].sum(axis=1)
    )

    # --------------------------------------------------------
    # Print hours where at least one file is abnormal
    # --------------------------------------------------------

    suspicious_all = pivot[
        pivot["NUMBER_OF_BAD_FILES"] > 0
    ]

    print()

    if suspicious_all.empty:

        print("NO CROSS-FILE HOURLY ANOMALY FOUND.")

    else:

        print(
            f"FOUND {len(suspicious_all)} "
            f"SUSPICIOUS TIMESTAMPS."
        )

        print()

        for t, row in suspicious_all.iterrows():

            bad = []

            for name in FILES:

                flag = f"{name}_BAD"

                if flag in row and bool(row[flag]):
                    bad.append(name)

            print(
                f"{t}  -->  "
                f"BAD: {', '.join(bad)}"
            )


# ============================================================
# DAILY SUMMARY
# ============================================================

print()
print()
print("=" * 110)
print("DAILY SUMMARY")
print("=" * 110)

for name, df in all_results.items():

    print()
    print(f"--- {name} ---")

    daily = df.groupby("date").agg(
        hours=("time", "count"),
        min_finite_pct=("finite_pct", "min"),
        mean_finite_pct=("finite_pct", "mean"),
        max_finite_pct=("finite_pct", "max"),
        suspicious_hours=("suspicious", "sum"),
    )

    for date, r in daily.iterrows():

        if r["suspicious_hours"] > 0:

            print(
                f"{date} | "
                f"hours={int(r['hours'])} | "
                f"min={r['min_finite_pct']:.6f}% | "
                f"mean={r['mean_finite_pct']:.6f}% | "
                f"suspicious_hours="
                f"{int(r['suspicious_hours'])}"
            )


# ============================================================
# SAVE CROSS-FILE RESULT
# ============================================================

cross_out = Path(
    "diagnose_1998_02_cross_file_hourly.csv"
)

pivot.to_csv(
    cross_out,
    encoding="utf-8-sig"
)


# ============================================================
# FINAL
# ============================================================

print()
print("=" * 110)
print("FORENSIC CHECK FINISHED")
print("=" * 110)

print(f"Hourly report : {OUT}")
print(f"Cross-file    : {cross_out}")

print()
print("IMPORTANT:")
print(
    "If the same timestamp is abnormal in T2M, D2M and SP, "
    "the problem is very likely in the source/download."
)
print(
    "If only one variable is abnormal, that specific file "
    "should be re-downloaded."
)
print("=" * 110)