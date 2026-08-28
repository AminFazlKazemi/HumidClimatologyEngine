# diagnose_1998_2_content.py

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

# ============================================================
# SETTINGS
# ============================================================

# درصدی از میدان که اگر یک مقدار واحد داشته باشد، مشکوک است
CONSTANT_FIELD_THRESHOLD = 99.99

# اگر بیش از این درصد سلول‌ها نسبت به ساعت قبل دقیقاً یکسان باشند
# ساعت مشکوک می‌شود.
IDENTICAL_THRESHOLD = 99.99

# جهش‌های بسیار بزرگ متوالی
# اینها فقط FLAG هستند و به تنهایی به معنی خرابی نیستند.
MAX_REASONABLE_JUMP = {
    "T2M": 15.0,       # K/hour
    "D2M": 15.0,       # K/hour
    "SP": 30000.0,     # Pa/hour
}

# ============================================================
# HELPERS
# ============================================================

def robust_stats(x):
    """Statistics on finite values only."""

    x = x[np.isfinite(x)]

    if x.size == 0:
        return {
            "min": np.nan,
            "max": np.nan,
            "mean": np.nan,
            "std": np.nan,
            "median": np.nan,
            "unique": 0,
        }

    return {
        "min": float(np.min(x)),
        "max": float(np.max(x)),
        "mean": float(np.mean(x)),
        "std": float(np.std(x)),
        "median": float(np.median(x)),
        "unique": int(np.unique(x).size),
    }


def analyze_file(name, path):

    print()
    print("=" * 120)
    print(f"{name} CONTENT FORENSICS")
    print("=" * 120)

    if not path.exists():
        print("!!! FILE NOT FOUND !!!")
        return None

    ds = xr.open_dataset(path)

    var = {
        "T2M": "t2m",
        "D2M": "d2m",
        "SP": "sp",
    }[name]

    if var not in ds:
        print(f"!!! VARIABLE {var} NOT FOUND !!!")
        ds.close()
        return None

    da = ds[var]

    times = pd.DatetimeIndex(ds.time.values)

    print(f"Variable : {var}")
    print(f"Shape    : {da.shape}")
    print(f"First    : {times[0]}")
    print(f"Last     : {times[-1]}")
    print(f"Hours    : {len(times)}")

    records = []

    previous = None

    for i, t in enumerate(times):

        field = np.asarray(
            da.isel(time=i).values,
            dtype=np.float64
        )

        finite_mask = np.isfinite(field)
        finite = field[finite_mask]

        total = field.size
        nfinite = finite.size
        nonfinite = total - nfinite

        if nfinite == 0:

            stats = {
                "min": np.nan,
                "max": np.nan,
                "mean": np.nan,
                "std": np.nan,
                "median": np.nan,
                "unique": 0,
            }

            constant_pct = 100.0
            identical_pct = 0.0
            mean_abs_change = np.nan
            max_abs_change = np.nan

        else:

            stats = robust_stats(finite)

            # ------------------------------------------------
            # CONSTANT FIELD
            # ------------------------------------------------

            counts = np.unique(
                finite,
                return_counts=True
            )[1]

            constant_pct = (
                100.0 * counts.max() / nfinite
            )

            # ------------------------------------------------
            # COMPARE WITH PREVIOUS HOUR
            # ------------------------------------------------

            if previous is None:

                identical_pct = np.nan
                mean_abs_change = np.nan
                max_abs_change = np.nan

            else:

                valid_pair = (
                    np.isfinite(field)
                    & np.isfinite(previous)
                )

                n_pair = int(valid_pair.sum())

                if n_pair == 0:

                    identical_pct = 0.0
                    mean_abs_change = np.nan
                    max_abs_change = np.nan

                else:

                    a = field[valid_pair]
                    b = previous[valid_pair]

                    diff = np.abs(a - b)

                    identical_pct = (
                        100.0
                        * np.count_nonzero(a == b)
                        / n_pair
                    )

                    mean_abs_change = float(
                        np.mean(diff)
                    )

                    max_abs_change = float(
                        np.max(diff)
                    )

        # ----------------------------------------------------
        # FLAGS
        # ----------------------------------------------------

        flags = []

        if nfinite == 0:
            flags.append("ALL_NONFINITE")

        if constant_pct >= CONSTANT_FIELD_THRESHOLD:
            flags.append("CONSTANT_FIELD")

        if (
            not np.isnan(identical_pct)
            and identical_pct >= IDENTICAL_THRESHOLD
        ):
            flags.append("IDENTICAL_TO_PREVIOUS")

        if (
            not np.isnan(max_abs_change)
            and max_abs_change > MAX_REASONABLE_JUMP[name]
        ):
            flags.append("LARGE_JUMP")

        records.append({
            "file": name,
            "time": t,
            "finite": nfinite,
            "nonfinite": nonfinite,
            "finite_pct": 100.0 * nfinite / total,
            "min": stats["min"],
            "max": stats["max"],
            "mean": stats["mean"],
            "std": stats["std"],
            "median": stats["median"],
            "unique": stats["unique"],
            "constant_pct": constant_pct,
            "identical_previous_pct": identical_pct,
            "mean_abs_change": mean_abs_change,
            "max_abs_change": max_abs_change,
            "flags": "|".join(flags),
        })

        previous = field

        # ----------------------------------------------------
        # PROGRESS
        # ----------------------------------------------------

        if (i + 1) % 100 == 0:
            print(f"  scanned {i + 1}/{len(times)} hours")

    ds.close()

    df = pd.DataFrame(records)

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print("-" * 120)
    print("CONTENT SUMMARY")
    print("-" * 120)

    print(
        f"mean of hourly means : "
        f"{df['mean'].mean():.8f}"
    )

    print(
        f"minimum hourly min   : "
        f"{df['min'].min():.8f}"
    )

    print(
        f"maximum hourly max   : "
        f"{df['max'].max():.8f}"
    )

    print(
        f"maximum constant %   : "
        f"{df['constant_pct'].max():.8f}"
    )

    print(
        f"maximum identical %  : "
        f"{df['identical_previous_pct'].max():.8f}"
    )

    print(
        f"maximum hour-to-hour "
        f"change              : "
        f"{df['max_abs_change'].max():.8f}"
    )

    # ========================================================
    # FLAGGED HOURS
    # ========================================================

    flagged = df[df["flags"] != ""].copy()

    print()
    print("-" * 120)
    print("FLAGGED HOURS")
    print("-" * 120)

    if flagged.empty:

        print("NO CONTENT ANOMALY FLAGGED.")

    else:

        print(
            f"{len(flagged)} flagged hours found."
        )

        for _, r in flagged.iterrows():

            print(
                f"{r['time']} | "
                f"mean={r['mean']:.6f} | "
                f"min={r['min']:.6f} | "
                f"max={r['max']:.6f} | "
                f"unique={r['unique']:,} | "
                f"constant={r['constant_pct']:.4f}% | "
                f"identical_prev="
                f"{r['identical_previous_pct']:.4f}% | "
                f"jump={r['max_abs_change']:.6f} | "
                f"FLAGS={r['flags']}"
            )

    return df


# ============================================================
# RUN
# ============================================================

results = {}

for name, path in FILES.items():

    result = analyze_file(name, path)

    if result is not None:
        results[name] = result


# ============================================================
# SAVE INDIVIDUAL REPORTS
# ============================================================

for name, df in results.items():

    filename = f"forensic_1998_02_{name}.csv"

    df.to_csv(
        filename,
        index=False,
        encoding="utf-8-sig"
    )

    print(
        f"\nSaved: {filename}"
    )


# ============================================================
# CROSS-FILE HOURLY COMPARISON
# ============================================================

print()
print()
print("=" * 120)
print("CROSS-FILE CONTENT COMPARISON")
print("=" * 120)

if results:

    combined = pd.concat(
        results.values(),
        ignore_index=True
    )

    # --------------------------------------------------------
    # Means
    # --------------------------------------------------------

    means = combined.pivot_table(
        index="time",
        columns="file",
        values="mean"
    )

    print()
    print("Hourly means:")
    print(means.to_string())

    # --------------------------------------------------------
    # Flag timestamps appearing in any file
    # --------------------------------------------------------

    flagged_times = {}

    for name, df in results.items():

        bad = df[df["flags"] != ""]

        for _, r in bad.iterrows():

            t = r["time"]

            if t not in flagged_times:
                flagged_times[t] = []

            flagged_times[t].append(
                f"{name}:{r['flags']}"
            )

    print()
    print("-" * 120)
    print("TIMESTAMPS WITH CONTENT FLAGS")
    print("-" * 120)

    if not flagged_times:

        print("NONE")

    else:

        for t in sorted(flagged_times):

            print(
                f"{t} --> "
                + " ; ".join(flagged_times[t])
            )


# ============================================================
# FINAL DIAGNOSTIC
# ============================================================

print()
print()
print("=" * 120)
print("FINAL DIAGNOSTIC")
print("=" * 120)

print("""
Interpretation:

1. ALL_NONFINITE
   Entire hourly field is invalid.

2. CONSTANT_FIELD
   Almost the entire spatial field contains one identical value.

3. IDENTICAL_TO_PREVIOUS
   The complete field is essentially repeated from the previous hour.

4. LARGE_JUMP
   A very large hour-to-hour change occurred.
   This is a warning, not automatically a corruption.

5. If the SAME timestamp is flagged in T2M + D2M + SP:
   investigate the original ERA5 download / source file.

6. If ONLY ONE variable is flagged:
   that specific input file becomes the primary suspect.

7. If NOTHING is flagged:
   the problem is probably NOT a simple corruption of the
   February input values, and we should investigate the
   processing engine / station matching / date handling next.
""")

print("=" * 120)
print("FORENSIC CONTENT CHECK FINISHED")
print("=" * 120)