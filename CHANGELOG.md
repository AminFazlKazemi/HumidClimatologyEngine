# Changelog

# v11.5 — Public Release

## [11.5.0] - Public Release

### Added
- Direct hourly processing of ERA5-Land T2m, D2m and Surface Pressure.
- Three temporal statistical levels: L1 daily pooled, L2 eight 3-hour bins, and L3 twenty-four hourly bins (33 bins total).
- Empirical moisture products for RH, vapor pressure (e), mixing ratio (r), and specific humidity (q).
- Online statistics through fourth order with counts, means, M2, M3, M4, extrema, missing counts and threshold counts.
- Pair-specific dependence statistics and empirical RH × q joint histogram support.
- Four decadal products plus FULL 1981–2020.
- Transactional checkpointing, durable COMMITTED-state tracking, restart/recovery, audit and merge-audit workflows.
- Detailed runtime progress reporting and real-data pilot execution.

### Validation
- D2m latitude-axis normalization to the T2m reference grid validated on real ERA5-Land input.
- January 2011 full-month pilot: 31/31 days PASS.
- Checkpoint flush PASS and checkpoint reopen PASS.

### Documentation
- Scientific and engineering reference, user/production runbook, v8 comparison, applications guide, analytical toolkit specification, and public README synchronized to v11.5.

### Notes
- v11.5 is the public release identity of HumidClimatologyEngine.
- Historical v8/v10/v7.5 materials remain historical references and are not the current release identity.

## [8.0.0] - Historical Release
- Single-pass empirical production architecture.
- Improved restart workflow.

## Historical packaging history

## 0.6.1
- Production bug-fix release for the v6 engine.
- Removed the duplicate `log_progress()` definition that caused `TypeError: multiple values for argument stage`.
- Standardized stage logging on `phase=` metadata.
- Fixed the merged-statistics shape validation f-string.
- Preserved the v6 checkpoint/config contract so existing 1981–2020 annual checkpoints can be reused.
- Added explicit source-level syntax/consistency validation to the patch workflow.

## 0.6.0
- Initial HumidClimatologyEngine repository packaging.
- Joint (T, Td, logP) climatology engine.
- Welford/Pébay accumulation and merging.
- Batch Monte Carlo and physical moisture derivation.
- Daily/annual checkpoints and SHA-256 verification.
- Detailed documentation and teaching notebook.
- ERA5-Land CDS daily-statistics downloader.
