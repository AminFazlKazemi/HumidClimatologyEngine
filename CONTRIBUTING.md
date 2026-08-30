# Contributing

HumidClimatologyEngine v11.5 is the current public release baseline.


Scientific changes should identify whether they affect data acquisition,
thermodynamics, statistical modeling, numerical accumulation, or performance.

Before a pull request:

```bash
python -m py_compile humid_climatology_engine_v11.5.py
pytest
```

Do not commit raw ERA5-Land files, outputs, checkpoints, credentials, or
machine-specific paths.
