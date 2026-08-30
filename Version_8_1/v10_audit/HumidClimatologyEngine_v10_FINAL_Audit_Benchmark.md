# HumidClimatologyEngine v10.0 FINAL — Audit & Benchmark

## Scope

This report audits the replacement v10 implementation intended to replace the historical v8 runner in the same repository filename:

`moisture_climatology_v8_0_FINAL_SINGLE_PASS_.py`

The v10 input contract remains the same as the historical v8 workflow: hourly ERA5-Land T2m, 2 m dew-point temperature, and surface pressure. The v10 software version is 10.0.0.

## P0 correctness/performance issues fixed

1. **Python-object quantile hot path removed.** The prior per-cell centroid/sketch object graph was removed from the production accumulation path. Quantiles can be produced by the analysis/second-pass layer instead of forcing millions of Python objects into RAM.
2. **State duplication across pairs removed.** Univariate moments are stored once per variable; pair states hold only pair-specific covariance/dependence information.
3. **Per-state pickle explosion removed.** The new checkpoint is one chunked NetCDF state shard per period/spatial block instead of separate pickle files for every period × level × pair × bin.
4. **Transaction truth separated from progress.** SQLite WAL + `synchronous=FULL` is the source of committed-work truth. Progress JSON is telemetry only.
5. **Rollback image made durable before OPEN transaction is recorded.** This avoids an OPEN record that has no rollback image after a crash.
6. **Day × spatial-block is the durable work unit.** A full day for one spatial block updates L1, L2, L3 and all configured periods/pairs in memory, then performs one durable transaction.
7. **NetCDF block descriptors are not left open for all blocks.** Active block files are closed at the end of the month/block work unit, avoiding excessive open file handles.
8. **Calendar metadata corrected.** Slot 59 is reserved, slot 60 is the Feb-28/Feb-29 composite, and slot 61 is Mar-01.
9. **Histogram behavior corrected.** Out-of-range observations are no longer silently clamped into edge bins; the histogram range is an explicit statistical support choice.
10. **Final outputs expanded.** Main scalar statistics, diagnostics, bivariate reference parameters, and RH×q empirical products are all represented by dedicated output files.
11. **Version identity unified.** Engine metadata is `10.0.0`; checkpoint/schema metadata is v10-specific.
12. **Validation tightened.** Time axes must be strictly hourly and grid coordinates must match and be monotonic.

## Crash-safety contract

```text
before-image durable
        ↓
OPEN transaction durable
        ↓
NetCDF scientific state update
        ↓
NetCDF sync()
        ↓
SQLite COMMIT durable
        ↓
rollback image cleanup
```

If power is lost before SQLite COMMIT, startup recovery restores the before-image and the work unit is replayed. A work unit is not considered complete because a progress percentage advanced.

## Tests performed

- Python bytecode compilation: **PASS**
- `compileall`: **PASS**
- AST parse: **PASS**
- Calendar contract: **PASS**
- Fourth-order moment batch calculations versus direct NumPy reference: **PASS**
- Pébay moment merge versus direct accumulation: **PASS**
- Physics smoke test for RH/r/q: **PASS**
- Transaction rollback smoke test: **PASS**
- Production state-update smoke test: **PASS**
- Histogram count conservation smoke test: **PASS**
- No executable `ProcessPoolExecutor` import/use: **PASS**
- No executable `ThreadPoolExecutor` import/use: **PASS**

## Hot-path benchmark

Synthetic in-memory benchmark of the vectorized moment core:

| Cells | 24-hour block | Approx. cells/s | Notes |
|---:|---:|---:|---|
| 16 | 10 repeated runs | 2.6–4.2 million | cache dominated |
| 2,048 | 10 repeated runs | 6.8–8.5 million | stable vectorized regime |
| 8,192 | 10 repeated runs | 4.4–7.8 million | large-block regime |

A full logical day/block hot-path simulation for **8,192 cells × 24 hours × 33 L1/L2/L3 bins** was measured at approximately **0.79 s** for the in-memory statistical update layer in one run, equivalent to roughly **8.18 million cell-hour-bin operations/s** under that synthetic workload.

## Important benchmark limitation

The execution environment used for this audit does not provide the `netCDF4` Python package and has no network access to install it. Therefore a real ERA5-Land NetCDF I/O benchmark and a complete 1981–2020 end-to-end benchmark could not be performed here.

Accordingly, this report does **not** claim a measured full-production wall-clock time on the user's machine.

## Release status

The source is code-complete for the v10 architecture and passes all environment-independent tests executed here. Final production acceptance still requires one controlled real-data pilot on the target ERA5-Land archive, followed by the complete 1981–2020 run and output audit.
