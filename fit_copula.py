import os
import sys
import json
import time
import numpy as np
import xarray as xr
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor, as_completed
from netCDF4 import Dataset

# افزودن مسیر فایل اصلی به sys.path
sys.path.append(str(Path(r"K:\kazemi\papers\temperature_interpolation\HumidClimatologyEngine")))

# ایمپورت توابع مورد نیاز
from moisture_climatology_v8_0_final_single_pass import (
    extract_centered_five_day_moisture_window,
    fit_window_models,
    DOY_COUNT,
    START_YEAR,
    END_YEAR,
    YEAR_DIR,
    OUTPUT_DIR,
    CONFIG_HASH
)

# ---------------------------
# تنظیمات
# ---------------------------
OUTPUT_NC = OUTPUT_DIR / "moisture_copula_parameters_1981_2020_v8_0.nc"
MAX_WORKERS = 4          # تعداد هسته‌های موازی
CHUNK_DAYS = 5           # تعداد روزهایی که هر کارگر پردازش می‌کند
LAT_RANGE = None         # مثلاً (0, 50) برای محدود کردن عرض جغرافیایی
LON_RANGE = None         # مثلاً (0, 100)
STRIDE_LAT = 1           # گام برداشتن از سلول‌ها (1 = همه)
STRIDE_LON = 1

# ---------------------------
# توابع کمکی
# ---------------------------
def get_grid_info():
    """خواندن مختصات شبکه از یک فایل نمونه."""
    sample_year = START_YEAR
    sample_file = next(Path(r"F:\Kazemi\era5\land\T2m").glob(f"*{sample_year}01*.nc"))
    with xr.open_dataset(sample_file) as ds:
        lat = ds.latitude.values
        lon = ds.longitude.values
    return lat, lon

def process_day(doy_idx, lat, lon):
    """پردازش یک روز اقلیمی برای همه‌ی سلول‌های انتخاب‌شده."""
    # نتایج این روز
    day_results = {
        'doy': doy_idx + 1,
        'lat_indices': [],
        'lon_indices': [],
        'rh_dist_name': [],
        'rh_params': [],
        'q_dist_name': [],
        'q_params': [],
        'rho': []
    }

    # تاریخ مرکزی این DOY (با استفاده از سال غیر کبیسه 2001 به عنوان مرجع)
    import datetime as dt
    base = dt.date(2001, 1, 1)
    target = base + dt.timedelta(days=doy_idx)
    target_str = target.strftime("%Y-%m-%d")  # فقط برای ساخت پنجره، سال مهم نیست چون فایل‌ها بر اساس ماه/سال واقعی انتخاب می‌شوند؟

    # برای هر سلول انتخابی
    for j, lat_idx in enumerate(range(0, len(lat), STRIDE_LAT)):
        if LAT_RANGE is not None and not (LAT_RANGE[0] <= lat[lat_idx] <= LAT_RANGE[1]):
            continue
        for i, lon_idx in enumerate(range(0, len(lon), STRIDE_LON)):
            if LON_RANGE is not None and not (LON_RANGE[0] <= lon[lon_idx] <= LON_RANGE[1]):
                continue

            try:
                # استخراج پنجره ۵ روزه
                window_vals, meta = extract_centered_five_day_moisture_window(
                    target_date=target_str,
                    lat_index=lat_idx,
                    lon_index=lon_idx
                )
                # برازش مدل‌ها و کاپولا
                model_result = fit_window_models(window_vals)
                rh_best = model_result["marginals"]["rh"]["best"]
                q_best = model_result["marginals"]["q"]["best"]
                rho = model_result["copula_rh_q"]["rho"]

                # ذخیره‌ی نتایج
                day_results['lat_indices'].append(lat_idx)
                day_results['lon_indices'].append(lon_idx)
                day_results['rh_dist_name'].append(rh_best["name"])
                day_results['rh_params'].append(rh_best["params"])
                day_results['q_dist_name'].append(q_best["name"])
                day_results['q_params'].append(q_best["params"])
                day_results['rho'].append(rho)

            except Exception as e:
                # در صورت خطا (مثلاً داده ناکافی) مقادیر NaN قرار دهید
                day_results['lat_indices'].append(lat_idx)
                day_results['lon_indices'].append(lon_idx)
                day_results['rh_dist_name'].append("")
                day_results['rh_params'].append({})
                day_results['q_dist_name'].append("")
                day_results['q_params'].append({})
                day_results['rho'].append(np.nan)

    return day_results

def save_results_to_netcdf(all_day_results, lat, lon):
    """ذخیره‌ی نتایج در فایل NetCDF."""
    # ابعاد
    nlat = len(range(0, len(lat), STRIDE_LAT))
    nlon = len(range(0, len(lon), STRIDE_LON))

    with Dataset(OUTPUT_NC, 'w', format='NETCDF4') as ds:
        ds.createDimension('doy', DOY_COUNT)
        ds.createDimension('latitude', nlat)
        ds.createDimension('longitude', nlon)
        ds.createDimension('param', 5)  # حداکثر ۵ پارامتر

        # متغیرهای مختصات
        ds.createVariable('doy', 'i2', ('doy',))[:] = np.arange(1, DOY_COUNT+1, dtype=np.int16)
        # مختصات انتخاب‌شده
        selected_lat = lat[::STRIDE_LAT] if LAT_RANGE is None else lat[(lat>=LAT_RANGE[0]) & (lat<=LAT_RANGE[1])][::STRIDE_LAT]
        selected_lon = lon[::STRIDE_LON] if LON_RANGE is None else lon[(lon>=LON_RANGE[0]) & (lon<=LON_RANGE[1])][::STRIDE_LON]
        ds.createVariable('latitude', 'f4', ('latitude',))[:] = selected_lat
        ds.createVariable('longitude', 'f4', ('longitude',))[:] = selected_lon

        # متغیرهای ذخیره‌سازی
        rh_dist = ds.createVariable('rh_dist_name', str, ('doy','latitude','longitude'))
        q_dist = ds.createVariable('q_dist_name', str, ('doy','latitude','longitude'))
        rh_params = ds.createVariable('rh_params', 'f4', ('doy','latitude','longitude','param'), fill_value=np.nan)
        q_params = ds.createVariable('q_params', 'f4', ('doy','latitude','longitude','param'), fill_value=np.nan)
        rho = ds.createVariable('gaussian_copula_rho', 'f4', ('doy','latitude','longitude'), fill_value=np.nan)

        # حلقه روی روزها و پر کردن
        for doy_idx, day_res in enumerate(all_day_results):
            for idx in range(len(day_res['lat_indices'])):
                lat_pos = day_res['lat_indices'][idx] // STRIDE_LAT
                lon_pos = day_res['lon_indices'][idx] // STRIDE_LON
                rh_dist[doy_idx, lat_pos, lon_pos] = day_res['rh_dist_name'][idx]
                q_dist[doy_idx, lat_pos, lon_pos] = day_res['q_dist_name'][idx]
                rho[doy_idx, lat_pos, lon_pos] = day_res['rho'][idx]

                # تبدیل پارامترها به آرایه‌ی ۵تایی
                rh_p = day_res['rh_params'][idx]
                q_p = day_res['q_params'][idx]
                # ترتیب پارامترها را بر اساس نوع توزیع می‌توانید استاندارد کنید؛ اینجا فقط ۵ مقدار اول را می‌گذاریم
                rh_vals = []
                for key in ['loc','scale','shape','a','b']:
                    rh_vals.append(rh_p.get(key, np.nan))
                q_vals = []
                for key in ['loc','scale','shape','a','b']:
                    q_vals.append(q_p.get(key, np.nan))
                rh_params[doy_idx, lat_pos, lon_pos, :] = rh_vals
                q_params[doy_idx, lat_pos, lon_pos, :] = q_vals

    print(f"نتایج در {OUTPUT_NC} ذخیره شد.")

def main():
    lat, lon = get_grid_info()
    print(f"ابعاد شبکه: {len(lat)}x{len(lon)}")
    print(f"تعداد سلول‌های انتخابی: {len(range(0, len(lat), STRIDE_LAT)) * len(range(0, len(lon), STRIDE_LON))}")

    # فهرست روزها برای پردازش
    doy_indices = list(range(DOY_COUNT))
    # تقسیم به chunk برای موازی‌سازی
    chunks = [doy_indices[i:i+CHUNK_DAYS] for i in range(0, len(doy_indices), CHUNK_DAYS)]

    all_day_results = [None] * DOY_COUNT

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {}
        for chunk in chunks:
            # هر chunk شامل چند روز است، اما برای سادگی هر روز جداگانه ارسال می‌شود
            for d in chunk:
                futures[executor.submit(process_day, d, lat, lon)] = d

        for future in as_completed(futures):
            d = futures[future]
            try:
                res = future.result()
                all_day_results[d] = res
                print(f"DOY {d+1} انجام شد.")
            except Exception as e:
                print(f"خطا در DOY {d+1}: {e}")

    save_results_to_netcdf(all_day_results, lat, lon)

if __name__ == "__main__":
    main()