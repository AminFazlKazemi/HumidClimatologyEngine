# Changelog

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
