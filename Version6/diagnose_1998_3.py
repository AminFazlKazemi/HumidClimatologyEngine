from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd

FILE = Path(
    r"F:\Kazemi\era5\land\Dew_Point_Temperature"
    r"\Dew_Point_Temperature199802.nc"
)

# ساعت‌های مشکوک از تست قبلی
TARGET_TIMES = pd.to_datetime([
    "1998-02-03 06:00",
    "1998-02-17 12:00",
    "1998-02-17 13:00",
    "1998-02-18 07:00",
    "1998-02-19 07:00",
])

# آستانه‌های تغییر مکانی
# K/hour
THRESHOLDS = [2.0, 5.0, 8.0, 10.0, 15.0]

print("=" * 120)
print("D2M SPATIAL EVENT FORENSICS — FEBRUARY 1998")
print("=" * 120)

ds = xr.open_dataset(FILE)
da = ds["d2m"]

times = pd.DatetimeIndex(ds["time"].values)

print("File:", FILE)
print("Shape:", da.shape)
print("Time:", times[0], "->", times[-1])

# ------------------------------------------------------------
# helper
# ------------------------------------------------------------

def get_field(t):
    idx = np.where(times == t)[0]
    if len(idx) != 1:
        raise RuntimeError(f"Timestamp not found uniquely: {t}")
    return np.asarray(
        da.isel(time=int(idx[0])).values,
        dtype=np.float64
    )


def summarize_difference(a, b, label):

    valid = np.isfinite(a) & np.isfinite(b)

    d = np.abs(a - b)

    print()
    print(f"--- {label} ---")

    if not np.any(valid):
        print("NO VALID OVERLAP")
        return

    dv = d[valid]

    print(f"Valid cells           : {valid.sum():,}")
    print(f"Mean abs difference   : {np.mean(dv):.6f} K")
    print(f"Median abs difference : {np.median(dv):.6f} K")
    print(f"95th percentile       : {np.percentile(dv,95):.6f} K")
    print(f"99th percentile       : {np.percentile(dv,99):.6f} K")
    print(f"Maximum difference    : {np.max(dv):.6f} K")

    for threshold in THRESHOLDS:

        pct = 100.0 * np.count_nonzero(
            valid & (d > threshold)
        ) / valid.sum()

        print(
            f"Cells with |Δ| > {threshold:5.1f} K : "
            f"{pct:10.6f}%"
        )

    # locations of largest changes
    masked = np.where(valid, d, -np.inf)

    flat_indices = np.argpartition(
        masked.ravel(),
        -10
    )[-10:]

    flat_indices = flat_indices[
        np.argsort(masked.ravel()[flat_indices])[::-1]
    ]

    lat = ds["latitude"].values
    lon = ds["longitude"].values

    print("\nTop 10 largest spatial changes:")

    for fi in flat_indices:

        iy, ix = np.unravel_index(
            fi,
            masked.shape
        )

        print(
            f"  lat={lat[iy]:8.3f}, "
            f"lon={lon[ix]:8.3f}, "
            f"|Δ|={masked[iy,ix]:10.6f} K"
        )


# ------------------------------------------------------------
# main analysis
# ------------------------------------------------------------

for target in TARGET_TIMES:

    print()
    print("=" * 120)
    print("TARGET:", target)
    print("=" * 120)

    pos = np.where(times == target)[0]

    if len(pos) != 1:
        print("TARGET NOT FOUND")
        continue

    i = int(pos[0])

    if i == 0 or i == len(times) - 1:
        print("Cannot compare both neighbors.")
        continue

    previous_time = times[i - 1]
    next_time = times[i + 1]

    previous = get_field(previous_time)
    current = get_field(target)
    following = get_field(next_time)

    print("Previous:", previous_time)
    print("Target  :", target)
    print("Next    :", next_time)

    # --------------------------------------------------------
    # target vs previous
    # --------------------------------------------------------

    summarize_difference(
        previous,
        current,
        "TARGET vs PREVIOUS HOUR"
    )

    # --------------------------------------------------------
    # next vs target
    # --------------------------------------------------------

    summarize_difference(
        current,
        following,
        "NEXT HOUR vs TARGET"
    )

    # --------------------------------------------------------
    # two-sided anomaly
    # --------------------------------------------------------

    valid = (
        np.isfinite(previous)
        & np.isfinite(current)
        & np.isfinite(following)
    )

    d1 = np.abs(current - previous)
    d2 = np.abs(following - current)

    # conservative two-sided anomaly:
    # both sides must show a large transition
    anomaly = valid & (d1 > 5.0) & (d2 > 5.0)

    pct = (
        100.0 * anomaly.sum() / valid.sum()
        if valid.sum() else np.nan
    )

    print()
    print("--- TWO-SIDED EVENT TEST ---")
    print(
        f"Cells with >5 K change on BOTH sides: "
        f"{anomaly.sum():,}"
    )
    print(
        f"Fraction of valid grid: {pct:.6f}%"
    )

    # --------------------------------------------------------
    # compare current value to average of neighbors
    # --------------------------------------------------------

    neighbor_mean = 0.5 * (
        previous + following
    )

    valid_mid = valid

    residual = np.abs(
        current - neighbor_mean
    )

    rv = residual[valid_mid]

    print()
    print("--- INTERPOLATED-MIDPOINT TEST ---")

    if rv.size:

        print(
            f"Mean |current - midpoint|   : "
            f"{np.mean(rv):.6f} K"
        )

        print(
            f"Median                       : "
            f"{np.median(rv):.6f} K"
        )

        print(
            f"95th percentile              : "
            f"{np.percentile(rv,95):.6f} K"
        )

        print(
            f"99th percentile              : "
            f"{np.percentile(rv,99):.6f} K"
        )

        for threshold in [2, 5, 8, 10, 15]:

            p = (
                100.0 *
                np.count_nonzero(
                    valid_mid & (residual > threshold)
                )
                / valid_mid.sum()
            )

            print(
                f"|residual| > {threshold:2d} K : "
                f"{p:.6f}%"
            )

    # --------------------------------------------------------
    # spatial extent classification
    # --------------------------------------------------------

    p15 = (
        100.0 *
        np.count_nonzero(valid & ((d1 > 15) | (d2 > 15)))
        / valid.sum()
    )

    p8 = (
        100.0 *
        np.count_nonzero(valid & ((d1 > 8) | (d2 > 8)))
        / valid.sum()
    )

    print()
    print("--- PRELIMINARY CLASSIFICATION ---")

    if p15 > 20:
        print("!!! EXTENSIVE DOMAIN-WIDE ANOMALY !!!")
    elif p8 > 10:
        print("!!! LARGE-SCALE ANOMALY !!!")
    elif p8 > 1:
        print("WARNING: REGIONAL / MODERATE ANOMALY")
    else:
        print("Likely localized / physically plausible event.")

# ------------------------------------------------------------
# close
# ------------------------------------------------------------

ds.close()

print()
print("=" * 120)
print("D2M SPATIAL FORENSICS FINISHED")
print("=" * 120)

print("""
Interpretation:

- If a very large fraction of the 301x301 grid changes abruptly,
  especially on both sides of the flagged hour, suspect the source
  file / retrieval.

- If the change is concentrated in a geographically coherent region,
  the event may be meteorologically real.

- If only a small fraction of cells show large changes, it is much
  less consistent with whole-file corruption.

- The strongest evidence of corruption is a broad, abrupt,
  spatially incoherent jump followed by an abrupt reversal.
""")