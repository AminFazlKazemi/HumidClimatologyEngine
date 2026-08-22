# Contributing

Scientific changes should identify whether they affect data acquisition,
thermodynamics, statistical modeling, numerical accumulation, or performance.

Before a pull request:

```bash
python -m compileall src
pytest
```

Do not commit raw ERA5-Land files, outputs, checkpoints, credentials, or
machine-specific paths.
