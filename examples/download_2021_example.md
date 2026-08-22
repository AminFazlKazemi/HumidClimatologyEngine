# 2021 download and legacy-data example

## Official CDS route

```bash
python scripts/download_era5land_daily_statistics.py \
  --years 2021 \
  --variables 2m_temperature 2m_dewpoint_temperature surface_pressure \
  --statistics daily_mean daily_minimum daily_maximum \
  --out-dir ./data/raw
```

## Legacy manual data

The project previously had manually calculated 2021 daily statistics. Keep that
archive with explicit metadata. Before substituting the official CDS-derived
product, compare daily definitions, time-zone convention, sample frequency,
missing values and leap-day handling.

A simple provenance record:

```yaml
year: 2021
source_dataset: ERA5-Land
processing: historical_manual_daily_aggregation
time_zone: UTC
notes: "Legacy 2021 statistics retained for continuity."
```
