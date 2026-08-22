# Download tools

The preferred acquisition route is the official CDS dataset
`derived-era5-land-daily-statistics`.

Example:

```bash
python scripts/download_era5land_daily_statistics.py \
  --years 2021 \
  --variables 2m_temperature 2m_dewpoint_temperature surface_pressure \
  --statistics daily_mean daily_minimum daily_maximum \
  --out-dir ./data/raw
```

Use the hourly fallback only when a custom aggregation or a variable not
represented by the daily-statistics product is required.
