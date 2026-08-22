# Production run

1. Download/audit inputs.
2. Run `pytest`.
3. Run the built-in mandatory tests.
4. Perform N convergence tests.
5. Generate annual Welford checkpoints.
6. Merge annual states.
7. Run daily Monte Carlo.
8. Validate RH/e/r/q ranges.
9. Inspect diagnostics.
10. Archive configuration, commit, input hashes and output hashes.
