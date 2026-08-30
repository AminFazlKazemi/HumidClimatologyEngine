# HumidClimatologyEngine v10.0 — Final Audit / Benchmark v2

## Source
`moisture_climatology_v8_0_FINAL_SINGLE_PASS_.py` (software identity: v10.0.0)

## Correctness fixes applied

1. Corrected the canonical 366-slot calendar metadata so slot 59 is reserved and slot 60 is the Feb-28/Feb-29 composite.
2. Added pair-specific M2 states (`M2_x`, `M2_y`) so pairwise correlation uses the variance of the same paired-valid population; marginal M2 can differ when masks differ.
3. Corrected pair covariance initialization for the first batch: Cxy now receives the batch centered cross-product rather than zero.
4. Corrected in-place M3/M4 merge logic to use the pre-update M2/M3 values, avoiding mutation-order errors.
5. Eliminated full 33-bin state copies for every temporal bin; the state is updated in-place bin-by-bin.
6. Replaced histogram `np.add.at` with flat-index `np.bincount` in the hot path.
7. Added stronger calendar, M3/M4, covariance, and pair-M2 self-tests.
8. Startup loads the grid before recovery so recovery can safely reopen checkpoint shards.
9. Added threshold/histogram metadata to final outputs and histogram `n_valid`.
10. Kept transaction truth in SQLite with `synchronous=FULL`; progress JSON is informational only.

## Performance architecture

- Single-process core; no process/thread executor.
- Streaming month input.
- One spatial block resident at a time.
- One day × block scientific transaction.
- In-place state updates rather than repeated complete-state copies.
- Histogram counting via `bincount`.
- NetCDF chunks aligned to the spatial block.
- L1/L2 histogram by default; L3 histogram is intentionally not persisted in the hot path.

## Automated checks in this environment

- Python compilation: PASS
- Calendar self-test: PASS
- Batch moments vs direct reference: PASS
- Pébay merge: PASS
- Physics: PASS
- State update including M3/M4/threshold/covariance: PASS
- Transaction rollback smoke test: PASS
- Synthetic in-memory benchmark: PASS

## Synthetic benchmark

| cells | mean batch-moment time (s) | cells/sec (24 time rows) |
|---:|---:|---:|
| 16 | 0.000087 | 4.40M |
| 2,048 | 0.007099 | 6.92M |
| 8,192 | 0.025955 | 7.57M |

These are in-memory NumPy measurements only; they are not ERA5 filesystem benchmarks.

## Remaining external validation

A true release acceptance run still requires the actual ERA5-Land folders on the target Windows machine. That run must verify:

- all 1981–2020 months exist and pass timestamp/grid/unit validation;
- the real 301×301 (or target) grid completes;
- interruption/restart equivalence on actual NetCDF checkpoints;
- FULL 1981–2020 equals the merge of the four decade products within declared floating-point tolerance;
- final NetCDF files reopen and pass metadata/QC audit.

These are environmental validations, not claims that can be proven from the source file alone.
