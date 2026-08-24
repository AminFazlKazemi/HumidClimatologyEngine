from pathlib import Path
import os
import traceback
import numpy as np

FILE = Path(
    r"C:\C\checkpoints_moisture_v7_5\years"
    r"\year_1998_083c99776caa07176970.part.nc"
)

print("=" * 100)
print("NETCDF FORENSIC CHECK")
print("=" * 100)

# ---------------------------------------------------------------------
# 1. FILE LEVEL
# ---------------------------------------------------------------------

print("\n[1] FILE")
print("-" * 100)

if not FILE.exists():
    print("ERROR: FILE DOES NOT EXIST")
    raise SystemExit(1)

size = FILE.stat().st_size

print("Path :", FILE)
print("Size :", f"{size:,} bytes")
print("Size :", f"{size / 1024**2:.2f} MiB")
print("Size :", f"{size / 1024**3:.3f} GiB")

if size < 100 * 1024 * 1024:
    print("WARNING: file is unusually small")


# ---------------------------------------------------------------------
# 2. netCDF4 LOW-LEVEL TEST
# ---------------------------------------------------------------------

print("\n[2] NETCDF4 LOW-LEVEL TEST")
print("-" * 100)

try:
    import netCDF4

    ds = netCDF4.Dataset(FILE, "r")

    print("OPEN: OK")
    print("Format:", ds.file_format)

    print("\nDimensions:")
    for name, dim in ds.dimensions.items():
        print(
            f"  {name:25s} "
            f"size={len(dim):12,d} "
            f"unlimited={dim.isunlimited()}"
        )

    print("\nVariables:")
    for name, var in ds.variables.items():
        print(
            f"  {name:35s} "
            f"dtype={str(var.dtype):10s} "
            f"shape={str(var.shape)}"
        )

    print("\nGlobal attributes:")
    for a in ds.ncattrs():
        try:
            print(f"  {a}: {getattr(ds, a)}")
        except Exception:
            print(f"  {a}: <UNREADABLE>")

except Exception as e:
    print("\n!!! NETCDF4 OPEN/HEADER ERROR !!!")
    print(type(e).__name__, ":", e)
    traceback.print_exc()

    print("\nThis strongly suggests file corruption/truncation.")
    raise SystemExit(2)


# ---------------------------------------------------------------------
# 3. VARIABLE FORENSICS
# ---------------------------------------------------------------------

print("\n[3] VARIABLE FORENSICS")
print("-" * 100)

variables = list(ds.variables.keys())

for name in variables:

    var = ds.variables[name]

    print("\n" + name)
    print("  dtype :", var.dtype)
    print("  shape :", var.shape)

    try:
        units = getattr(var, "units", "")
        print("  units :", units)
    except Exception:
        pass

    try:
        long_name = getattr(var, "long_name", "")
        if long_name:
            print("  name  :", long_name)
    except Exception:
        pass

    # Skip huge multidimensional variables here;
    # coordinates and small variables are tested fully.
    try:
        if var.size <= 10_000_000:
            x = var[:]

            if np.ma.isMaskedArray(x):
                data = x.compressed()
                masked = np.count_nonzero(x.mask)
            else:
                data = np.asarray(x)
                masked = 0

            if np.issubdtype(data.dtype, np.number):

                finite = np.isfinite(data)

                n_total = data.size
                n_finite = np.count_nonzero(finite)
                n_nan = np.count_nonzero(np.isnan(data))
                n_inf = np.count_nonzero(np.isinf(data))

                print("  total :", f"{n_total:,}")
                print("  finite:", f"{n_finite:,}")
                print("  masked:", f"{masked:,}")
                print("  NaN   :", f"{n_nan:,}")
                print("  Inf   :", f"{n_inf:,}")

                if n_finite:
                    d = data[finite]
                    print("  min   :", np.min(d))
                    print("  max   :", np.max(d))
                    print("  mean  :", np.mean(d))

        else:
            print(
                f"  FULL ARRAY TEST SKIPPED "
                f"({var.size:,} values)"
            )

    except Exception as e:
        print("  !!! READ ERROR !!!")
        print(" ", type(e).__name__, ":", e)


# ---------------------------------------------------------------------
# 4. IDENTIFY IMPORTANT VARIABLES
# ---------------------------------------------------------------------

print("\n[4] IMPORTANT VARIABLES")
print("-" * 100)

keywords = {
    "pressure": [
        "pressure", "pres", "sp", "msl", "surface_pressure"
    ],
    "temperature": [
        "temperature", "temp", "t2m", "2m_temperature"
    ],
    "dewpoint": [
        "dew", "dewpoint", "dew_point", "d2m", "2m_dewpoint"
    ],
}

for group, keys in keywords.items():

    print(f"\n{group.upper()}")

    found = []

    for name in variables:
        low = name.lower()

        if any(k in low for k in keys):
            found.append(name)

    if found:
        for name in found:
            print("  ", name)
    else:
        print("   NOT FOUND")


# ---------------------------------------------------------------------
# 5. FULL READ TEST FOR IMPORTANT VARIABLES
# ---------------------------------------------------------------------

print("\n[5] FULL READ TEST")
print("-" * 100)

important_names = []

for name in variables:

    low = name.lower()

    if any(
        k in low
        for k in (
            "pressure",
            "pres",
            "temperature",
            "temp",
            "t2m",
            "dew",
            "d2m",
        )
    ):
        important_names.append(name)

for name in sorted(set(important_names)):

    var = ds.variables[name]

    print(f"\nTesting: {name}")
    print("  shape:", var.shape)

    try:

        # Read first element
        if var.ndim > 0:
            first_index = (0,) * var.ndim
            first = var[first_index]
            print("  first element: OK")

        # Read last element
        if var.ndim > 0:
            last_index = tuple(s - 1 for s in var.shape)
            last = var[last_index]
            print("  last element : OK")

        # Read first time slice if possible
        if var.ndim >= 1 and var.shape[0] > 0:

            idx = [slice(None)] * var.ndim
            idx[0] = 0

            _ = var[tuple(idx)]

            print("  first slice  : OK")

        # Read last time slice
        if var.ndim >= 1 and var.shape[0] > 0:

            idx = [slice(None)] * var.ndim
            idx[0] = var.shape[0] - 1

            _ = var[tuple(idx)]

            print("  last slice   : OK")

    except Exception as e:

        print("  !!! DATA READ FAILURE !!!")
        print(" ", type(e).__name__, ":", e)


# ---------------------------------------------------------------------
# 6. CHUNKED READ
# ---------------------------------------------------------------------

print("\n[6] CHUNKED READ TEST")
print("-" * 100)

for name in sorted(set(important_names)):

    var = ds.variables[name]

    if var.ndim == 0:
        continue

    print(f"\n{name}")

    n = var.shape[0]

    # 10 positions through the complete time dimension
    positions = np.linspace(
        0,
        max(0, n - 1),
        min(10, n),
        dtype=int
    )

    failed = False

    for i in positions:

        try:

            idx = [slice(None)] * var.ndim
            idx[0] = int(i)

            _ = var[tuple(idx)]

            print(f"  slice {i:8,d}: OK")

        except Exception as e:

            failed = True

            print(
                f"  slice {i:8,d}: FAILED -> "
                f"{type(e).__name__}: {e}"
            )

    if not failed:
        print("  RESULT: ALL TESTED SLICES OK")


# ---------------------------------------------------------------------
# 7. xarray TEST
# ---------------------------------------------------------------------

print("\n[7] XARRAY TEST")
print("-" * 100)

try:

    import xarray as xr

    print("Opening with xarray...")

    xds = xr.open_dataset(
        FILE,
        engine="netcdf4",
        decode_times=True,
        cache=False,
    )

    print("XARRAY OPEN: OK")
    print(xds)

    print("\nDimensions:")
    print(dict(xds.sizes))

    print("\nCoordinates:")
    for name in xds.coords:
        print(
            f"  {name}: "
            f"shape={xds[name].shape}, "
            f"dtype={xds[name].dtype}"
        )

    print("\nTime test:")

    time_candidates = [
        x for x in xds.coords
        if x.lower() in ("time", "valid_time", "date")
        or "time" in x.lower()
    ]

    for name in time_candidates:

        try:
            print(" ", name)
            print("   first:", xds[name].values[0])
            print("   last :", xds[name].values[-1])
        except Exception as e:
            print("   ERROR:", e)

except Exception as e:

    print("\n!!! XARRAY ERROR !!!")
    print(type(e).__name__, ":", e)
    traceback.print_exc()


# ---------------------------------------------------------------------
# 8. CLOSE
# ---------------------------------------------------------------------

try:
    ds.close()
except Exception:
    pass

print("\n" + "=" * 100)
print("FORENSIC CHECK FINISHED")
print("=" * 100)