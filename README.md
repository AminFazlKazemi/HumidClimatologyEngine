# HumidClimatologyEngine

> **Production-grade multivariate moisture climatology engine for ERA5-Land**  
> Joint statistical modeling of **2 m temperature (T2m), 2 m dew-point temperature (Td2m), and surface pressure (P)** in the transformed state space **(T, Td, ln P)**, followed by physically constrained Monte Carlo propagation to **relative humidity, vapor pressure, mixing ratio, and specific humidity**.

[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-pytest-0A9EDC?logo=pytest&logoColor=white)](tests/)
[![ERA5--Land](https://img.shields.io/badge/data-ERA5--Land-1F4E79)](https://cds.climate.copernicus.eu/)
[![NetCDF](https://img.shields.io/badge/output-NetCDF4-4B8BBE)](https://www.unidata.ucar.edu/software/netcdf/)

---

## 1. Why this project exists

HumidClimatologyEngine is designed for research workflows where a moisture climatology must be more than a simple arithmetic average.

Instead of independently averaging humidity-related variables, the engine preserves the empirical dependence between:

\[
X=(T,T_d,\ln P)
\]

and propagates that joint uncertainty through nonlinear thermodynamic transformations.

The production chain is:

```text
ERA5-Land
   │
   ├── 2 m temperature (T)
   ├── 2 m dew-point temperature (Td)
   └── surface pressure (P)
          │
          ▼
paired-valid observations
          │
          ▼
Welford / covariance accumulation
          │
          ▼
annual sufficient-statistic checkpoints
          │
          ▼
parallel merge across years
          │
          ▼
DOY × grid-cell multivariate state
          │
          ▼
PSD-controlled covariance / Cholesky
          │
          ▼
streaming Monte Carlo in (T, Td, lnP)
          │
          ▼
physical moisture transformation
          │
          ├── RH
          ├── vapor pressure
          ├── mixing ratio
          └── specific humidity
          │
          ▼
Pébay higher moments + MC uncertainty
          │
          ▼
daily checkpoints
          │
          ▼
streaming NetCDF finalization
          │
          ▼
validated climatology + diagnostics
```

The project prioritizes:

- **scientific traceability**
- **numerical defensiveness**
- **bounded memory use**
- **restartability**
- **deterministic reproducibility**
- **explicit diagnostics**
- **independent validation**
- **long-running workstation safety**

---

# 2. Scientific scope

The engine currently implements a trivariate statistical representation:

\[
X =
\begin{bmatrix}
T\\
T_d\\
\ln P
\end{bmatrix}
\sim \mathcal N(\mu,\Sigma)
\]

for each climatological day and grid cell.

The sampled state is then transformed into moisture quantities.

### Relative humidity

\[
RH = 100\frac{e_s(T_d)}{e_s(T)}
\]

### Vapor pressure

\[
e=e_s(T_d)
\]

### Mixing ratio

\[
r=\frac{0.622e}{P-e}
\]

### Specific humidity

\[
q=\frac{r}{1+r}
\]

The exact saturation-vapor-pressure implementation is part of the production code and should be treated as part of the model definition when reproducing published results.

---

# 3. Core design principles

## 3.1 Never hide a numerical failure

The covariance layer explicitly checks positive-semidefiniteness and Cholesky feasibility.

The implementation does **not** silently replace an invalid covariance matrix with an identity matrix.

This matters scientifically: an identity fallback would erase the estimated dependence structure and could produce apparently valid but scientifically different results.

Diagnostics retain information about covariance validity and minimum eigenvalues.

---

## 3.2 Never build the giant Monte Carlo cube

The engine does not construct a tensor such as:

```text
N_SAMPLES × latitude × longitude
```

for an entire day.

Instead:

```text
day
  → cell chunk
      → sample batch
          → transform
              → online statistics
```

The dominant temporary workload is therefore controlled approximately by:

```text
CELL_CHUNK_SIZE × SAMPLE_BATCH_SIZE
```

rather than the entire global grid.

This is one of the central memory-safety features of the project.

---

## 3.3 Never keep the full multi-decadal raw history in memory

The accumulation stage works through the source data incrementally.

For the default 1981–2020 configuration:

```text
year
  → month
      → day
          → valid paired observations
              → sufficient statistics
```

At the end of a year, the sufficient statistics are checkpointed.

Raw 40-year histories are not retained by the statistical accumulator.

---

## 3.4 Restart is a first-class feature

The workflow is deliberately designed for jobs that may run for many hours or days.

It uses:

- annual checkpoints
- daily checkpoints
- configuration hashes
- SHA-256 verification
- atomic checkpoint writes
- deterministic random streams
- compatibility metadata

A completed valid checkpoint can be reused after interruption instead of recomputing the corresponding stage.

---

# 4. Calendar convention

The production climatology has **366 slots**.

The project uses the following explicit convention:

| Climatological DOY | Meaning |
|---:|---|
| 1–58 | Jan 1 → Feb 27 |
| **59** | **Reserved** |
| **60** | **Feb 28 + Feb 29 composite** |
| 61–366 | Mar 1 → Dec 31 |

This means:

```text
Leap year:
Feb 28 → DOY 60
Feb 29 → DOY 60
Mar 1  → DOY 61

Non-leap year:
Feb 28 → DOY 60
Mar 1  → DOY 61
```

The reserved slot is intentional and must not be silently reinterpreted downstream.

The repository contains explicit tests for this mapping.

---

# 5. Statistical engine

## 5.1 Paired-valid observations

An observation contributes to the multivariate state only when the required variables form a valid joint sample.

The core condition is:

```text
T finite
Td finite
P finite
P > 0
```

This avoids constructing a covariance matrix from mismatched populations.

---

## 5.2 Welford accumulation

The engine uses numerically stable online accumulation for:

- count
- mean
- variance-related quantities
- covariance terms

The method avoids repeatedly storing all historical observations simply to calculate final moments.

---

## 5.3 Parallel merge

Annual sufficient-statistic states can be merged using parallel/merge formulas.

Conceptually:

```text
1981 state ─┐
1982 state ─┤
1983 state ─┤
...         ├──► merged climatological state
2020 state ─┘
```

This is important for restartability and memory efficiency.

---

## 5.4 Higher moments

The transformed moisture variables retain:

- mean
- sample standard deviation
- bias-corrected skewness
- Fisher excess kurtosis

The implementation uses central-moment accumulators based on Pébay-style online/mergeable formulas.

Interpretation of Fisher excess kurtosis:

```text
≈ 0  → approximately normal tail weight
> 0  → heavier tails
< 0  → lighter tails
```

The final conventions are aligned with the intended SciPy-style definitions:

```python
scipy.stats.skew(..., bias=False)
scipy.stats.kurtosis(..., bias=False, fisher=True)
```

---

# 6. Monte Carlo engine

The default production configuration is:

```text
N_SAMPLES          = 5000
CELL_CHUNK_SIZE    = 1024
SAMPLE_BATCH_SIZE  = 256
MAX_WORKERS        = 2
RANDOM_SEED        = 20260821
```

These are configuration choices, not universal scientific constants.

A recommended sensitivity experiment is:

```text
N = 500
N = 1000
N = 2000
N = 5000
N = 10000
```

Then compare:

- mean RH
- RH standard deviation
- RH skewness
- RH kurtosis
- mean q
- mean r
- mean vapor pressure
- diagnostic failure rates

A value such as 5000 should be retained only when convergence is adequate for the intended scientific conclusions.

---

# 7. Monte Carlo uncertainty

The diagnostic product reports Monte Carlo standard-error estimates for selected means:

```text
mc_se_mean_rh
mc_se_mean_e
mc_se_mean_r
mc_se_mean_q
```

For a Monte Carlo mean:

\[
SE(\bar{x})\approx \frac{s}{\sqrt N}
\]

where \(s\) is the Monte Carlo standard deviation and \(N\) is the number of valid realizations.

This provides an explicit way to distinguish:

```text
scientific variability
        from
Monte Carlo sampling error
```

---

# 8. Physical and numerical safeguards

The engine contains multiple layers of defensive checks.

### Input checks

- required dimensions
- time coordinate compatibility
- latitude compatibility
- longitude compatibility
- month completeness
- duplicate month detection
- pressure positivity
- finite-value screening

### Statistical checks

- minimum observation counts
- covariance validity
- PSD tolerance
- eigenvalue diagnostics
- Cholesky feasibility

### Transformation checks

- pressure validity
- vapor-pressure validity
- denominator validity for mixing ratio
- finite transformed values
- physical bounds

### Output checks

- dimensions
- coordinates
- variable presence
- bounds
- diagnostic consistency
- checkpoint integrity

The philosophy is:

> **A calculation should fail loudly when its assumptions are violated, rather than silently producing plausible-looking numbers.**

---

# 9. Checkpoint architecture

Two main checkpoint layers are used.

## 9.1 Annual checkpoints

A year produces:

```text
year_YYYY_<config_hash>.npz
year_YYYY_<config_hash>.json
```

The metadata records information such as:

- year
- schema version
- configuration hash
- SHA-256
- grid shape
- timestamp

A checkpoint is reused only when its integrity and compatibility checks succeed.

---

## 9.2 Daily checkpoints

Daily Monte Carlo results use:

```text
day_001_<config_hash>.npz
day_001_<config_hash>.json
...
day_366_<config_hash>.npz
day_366_<config_hash>.json
```

This allows a partially completed climatology to resume from the last valid daily state.

---

## 9.3 Atomic writes

Checkpoint writes are designed around temporary files followed by atomic replacement.

The intended failure model is:

```text
calculate
   ↓
write temporary artifact
   ↓
flush / synchronize
   ↓
atomic replace
   ↓
write metadata / checksum
```

A process termination during writing should therefore be much less likely to leave a misleadingly complete checkpoint.

---

# 10. Reproducibility model

Reproducibility metadata includes core numerical configuration such as:

```text
start year
end year
DOY count
Monte Carlo sample count
cell chunk size
sample batch size
worker count
random seed
schema version
PSD tolerances
```

The configuration is hashed.

The implementation also records software/environment metadata such as:

- Python version
- NumPy version
- SciPy version
- xarray version
- netCDF4 version
- script SHA-256
- configuration hash
- creation time

For serious publication work, preserve the exact input data version and preprocessing recipe as well.

---

# 11. Input data: ERA5-Land

The project is designed for ERA5-Land data obtained from the **Copernicus Climate Data Store (CDS)** or an equivalent reproducible source.

Relevant variables are:

```text
2m temperature
2m dewpoint temperature
surface pressure
```

Common ERA5-Land short names are:

```text
t2m
d2m
sp
```

The current production reader expects the repository's configured variable names and file layout; if your NetCDF files use different names, normalize them during preprocessing or adapt the reader layer.

> **Important:** the engine expects daily fields for climatology construction. A monthly climatological mean is not equivalent to an individual daily field.

---

# 12. Recommended preprocessing

A transparent workflow is:

```text
ERA5-Land hourly
       │
       ▼
monthly local NetCDF
       │
       ▼
daily aggregation
       │
       ├── T2m daily
       ├── Td2m daily
       └── P daily
       │
       ▼
coordinate/unit validation
       │
       ▼
HumidClimatologyEngine
```

Recommended directory organization:

```text
era5/
└── land/
    └── daily/
        ├── T2m/
        │   ├── era5land_t2m_1981_01.nc
        │   ├── era5land_t2m_1981_02.nc
        │   └── ...
        ├── Dew_Point_Temperature/
        │   ├── era5land_d2m_1981_01.nc
        │   ├── era5land_d2m_1981_02.nc
        │   └── ...
        └── Surface_Pressure/
            ├── era5land_sp_1981_01.nc
            ├── era5land_sp_1981_02.nc
            └── ...
```

For 1981–2020:

```text
40 years × 12 months = 480 monthly files per variable
```

The file indexer is intentionally strict about monthly completeness and duplicate months.

---

# 13. Input contract

Before a production run, verify:

```python
import xarray as xr

ds = xr.open_dataset("example.nc")

print(ds)
print(ds.data_vars)
print(ds.dims)
print(ds.latitude.values[:5])
print(ds.longitude.values[:5])
print(ds.time.values[:5])
```

The T, Td, and P datasets must have compatible:

- time coordinates
- latitude coordinates
- longitude coordinates
- dimensions

The engine validates these relationships before combining the variables.

---

# 14. Installation

## Requirements

Python:

```text
>= 3.11
```

Install the project dependencies:

```bash
pip install -r requirements.txt
```

or:

```bash
pip install .
```

For development/testing:

```bash
pip install -e ".[test]"
```

For notebook support:

```bash
pip install -e ".[notebook]"
```

For all optional development/notebook dependencies:

```bash
pip install -e ".[all]"
```

Core dependencies include:

```text
numpy
scipy
xarray
netCDF4
tqdm
psutil
cdsapi
PyYAML
```

---

# 15. Repository layout

The current repository intentionally keeps the main production implementation compact:

```text
HumidClimatologyEngine/
├── README.md
├── CHANGELOG.md
├── CITATION.cff
├── CONTRIBUTING.md
├── SECURITY.md
├── LICENSE
├── pyproject.toml
├── requirements.txt
├── environment.yml
│
├── configs/
│   └── example_paths.yaml
│
├── docs/
│   └── methodology.md
│
├── examples/
│   ├── download_2021_example.md
│   ├── production_run.md
│   └── HumidClimatologyEngine_tutorial.ipynb
│
├── scripts/
│   ├── download_era5land_daily_statistics.py
│   ├── download_era5land_hourly.py
│   └── README.md
│
├── src/
│   ├── __init__.py
│   └── moisture_climatology_v6.py
│
├── moisture_climatology.py
│
└── tests/
    ├── README.md
    └── test_repository.py
```

### Current architecture

The current release is deliberately centered on:

```text
moisture_climatology.py
        │
        ▼
src/moisture_climatology_v6.py
```

The long-term architectural direction can split physics, statistics, I/O, checkpoints, Monte Carlo, calendar handling, and validation into separate modules, but the README describes the implementation that actually ships rather than pretending that a future modular architecture already exists.

---

# 16. Running the engine

The simplest production entry point is:

```bash
python moisture_climatology.py
```

The wrapper delegates to:

```text
src/moisture_climatology_v6.py
```

The main configuration currently lives in the production module.

For reproducible projects, prefer storing the final configuration in version control and recording its configuration hash with the generated product.

---

# 17. Windows and Spyder

The code is written with Windows multiprocessing considerations in mind.

A typical Spyder launch is:

```python
runfile(
    r"K:\path\to\HumidClimatologyEngine\moisture_climatology.py",
    wdir=r"K:\path\to\HumidClimatologyEngine",
)
```

For long production jobs, executing the repository as a script is preferable to repeatedly copying individual code cells into an interactive namespace.

Because each worker can consume significant memory, begin conservatively:

```python
MAX_WORKERS = 2
```

Then benchmark before increasing parallelism.

---

# 18. Memory strategy

The processing hierarchy is:

```text
40 years
  → one year
      → one month
          → one day
              → one cell chunk
                  → one sample batch
```

The engine therefore avoids:

```text
all years × all days × all cells × all samples
```

being materialized at once.

This matters especially on Windows systems where process workers may have substantial independent memory footprints.

### Practical tuning knobs

```python
CELL_CHUNK_SIZE
SAMPLE_BATCH_SIZE
MAX_WORKERS
N_SAMPLES
```

A safe tuning strategy is:

1. start with the default values;
2. observe RAM during a representative workload;
3. increase workers only if memory headroom remains large;
4. benchmark larger cell chunks;
5. benchmark larger sample batches;
6. compare runtime against peak RAM;
7. retain the configuration that is both stable and scientifically adequate.

---

# 19. Progress and observability

The runtime logger is designed to provide operational visibility during long runs.

Where available, it reports:

```text
RAM used
RAM available
RAM percentage
CPU usage
stage/phase
year
month
day
DOY
chunk
Monte Carlo progress
elapsed time
```

Telemetry is deliberately non-fatal: failure of a progress/telemetry probe must not interrupt the scientific calculation.

---

# 20. Output products

## Main climatology

Default:

```text
moisture_climatology_1981_2020.nc
```

The main product contains, for each climatological day and grid cell:

### Relative humidity

```text
mean_rh
std_rh
skew_rh
kurt_rh
```

### Vapor pressure

```text
mean_vapor_pressure
std_vapor_pressure
skew_vapor_pressure
kurt_vapor_pressure
```

### Mixing ratio

```text
mean_mixing_ratio
std_mixing_ratio
skew_mixing_ratio
kurt_mixing_ratio
```

### Specific humidity

```text
mean_specific_humidity
std_specific_humidity
skew_specific_humidity
kurt_specific_humidity
```

---

## Diagnostic product

Default:

```text
moisture_climatology_diagnostics_1981_2020.nc
```

Representative diagnostics include:

```text
supersaturation_fraction
invalid_e_over_p_fraction
invalid_covariance_fraction
min_eigenvalue
mc_se_mean_rh
mc_se_mean_e
mc_se_mean_r
mc_se_mean_q
valid_sample_count
corr_T_Td
corr_T_logP
corr_Td_logP
valid_observation_count
```

These diagnostics are essential for distinguishing:

```text
real climatological behavior
        from
numerical pathology
        from
covariance failure
        from
insufficient Monte Carlo sampling
        from
missing observations
```

---

# 21. Reading the NetCDF

```python
import xarray as xr

ds = xr.open_dataset(
    r"C:\c\moisture_climatology_1981_2020.nc"
)

print(ds)
```

Read a field:

```python
rh = ds["mean_rh"].sel(doy=200)
```

Specific humidity:

```python
q = ds["mean_specific_humidity"].sel(doy=200)
```

Mixing ratio:

```python
r = ds["mean_mixing_ratio"].sel(doy=200)
```

Composite February day:

```python
rh_feb = ds["mean_rh"].sel(doy=60)
```

Remember:

```text
DOY 59 = reserved
DOY 60 = Feb-28/Feb-29 composite
```

---

# 22. Basic visualization

```python
import matplotlib.pyplot as plt
import xarray as xr

ds = xr.open_dataset(
    r"C:\c\moisture_climatology_1981_2020.nc"
)

field = ds["mean_rh"].sel(doy=200)

field.plot(
    figsize=(10, 7),
    robust=True,
)

plt.title("Mean Relative Humidity — DOY 200")
plt.tight_layout()
plt.show()
```

For publication maps, users should add an appropriate cartographic projection, geographic boundaries, units, colorbar labeling, and figure metadata.

---

# 23. Point climatology

```python
lat0 = 35.0
lon0 = 51.0

series = (
    ds["mean_rh"]
    .sel(latitude=lat0, longitude=lon0, method="nearest")
)

series.plot(figsize=(12, 5))
```

This produces the annual climatological cycle at the nearest grid cell.

---

# 24. Example production workflow

A defensible production sequence is:

```text
1. Preserve raw ERA5-Land inputs
2. Document the exact dataset/version
3. Build daily T/Td/P fields
4. Validate units and coordinates
5. Run repository tests
6. Run a small synthetic/physics check
7. Run Monte Carlo convergence experiments
8. Run the target climatology
9. Inspect annual checkpoints
10. Inspect daily checkpoints
11. Inspect diagnostics
12. Validate the final NetCDF
13. Compare selected cells against an independent implementation
14. Archive configuration + hashes + environment
15. Publish code + metadata + outputs as appropriate
```

---

# 25. Validation strategy

Validation should occur at several independent levels.

## 25.1 Repository tests

The current test suite includes checks for:

- leap-day mapping
- moisture physical transformation
- source-level progress-function consistency

Run:

```bash
pytest
```

---

## 25.2 Syntax validation

```bash
python -m compileall src
```

---

## 25.3 Physics ground truth

Use an intentionally simple independent implementation for selected scalar test cases.

The goal is to avoid validating an optimized vectorized implementation solely against itself.

---

## 25.4 Statistical validation

Compare online/mergeable calculations with trusted reference calculations on controlled synthetic arrays.

Validate:

- mean
- variance
- covariance
- skewness
- kurtosis

---

## 25.5 Monte Carlo convergence

Compare multiple sample sizes and quantify:

- absolute differences
- relative differences
- spatial maximum error
- spatial mean error
- diagnostic failure rates

---

## 25.6 Final-product validation

At minimum check:

```text
coordinates
dimensions
finite values
RH bounds
q bounds
r validity
vapor pressure positivity
valid sample counts
covariance diagnostics
Monte Carlo standard errors
```

---

# 26. Restart semantics

The intended execution graph is:

```text
START
  │
  ├── valid annual checkpoint?
  │       ├── yes → reuse
  │       └── no  → calculate
  │
  ▼
MERGE ANNUAL STATES
  │
  ├── valid daily checkpoint?
  │       ├── yes → reuse
  │       └── no  → calculate
  │
  ▼
FINALIZE NETCDF
  │
  ▼
VALIDATE OUTPUT
  │
  ▼
DONE
```

This makes the workflow suitable for:

- workstation restarts
- power interruptions
- controlled performance experiments
- memory tuning
- long Monte Carlo jobs

---

# 27. Checkpoint safety rules

Do not mix checkpoints produced with incompatible:

- model definitions
- schema versions
- configuration hashes
- random seeds
- Monte Carlo sample counts
- PSD tolerances

Do not delete annual checkpoints until the final product has been independently validated.

For a deliberately different experiment, use a distinct output/checkpoint location or ensure that the configuration identity changes.

---

# 28. Scientific limitations

This engine is powerful, but it is not assumption-free.

## 28.1 Gaussian joint model

The model assumes an approximately multivariate-normal representation in:

\[
(T,T_d,\ln P)
\]

Real atmospheric states may deviate from this assumption.

This can matter especially for:

- extremes
- tails
- strongly skewed regimes
- mountainous environments
- cold conditions
- unusual pressure regimes
- compound events

---

## 28.2 Reanalysis uncertainty

ERA5-Land is a reanalysis product, not direct truth.

Uncertainty can arise from:

- model physics
- observations assimilated into the reanalysis
- spatial representativeness
- temporal aggregation
- terrain representation

---

## 28.3 Monte Carlo error

A finite sample count produces finite sampling error.

Increasing `N_SAMPLES` reduces Monte Carlo error but increases runtime.

---

## 28.4 Thermodynamic approximation

The exact saturation-vapor-pressure formulation and phase treatment are model choices.

These choices should be reported in publications and tested in sensitivity experiments when scientifically relevant.

---

## 28.5 Calendar convention

The Feb-28/Feb-29 pooling convention is deliberate but not universal.

Alternative climatological calendar definitions may produce slightly different results and should be treated as sensitivity dimensions when necessary.

---

# 29. Recommended sensitivity matrix

For a research publication, consider evaluating:

| Dimension | Suggested alternatives |
|---|---|
| Monte Carlo | 500 / 1000 / 2000 / 5000 / 10000 |
| Pressure state | P / ln P |
| Saturation formulation | water-only / phase-aware |
| February treatment | pooled / alternative Feb treatment |
| Joint distribution | Gaussian / mixture / copula |
| Spatial domain | full / regional |
| Temporal period | baseline / sensitivity period |

The final method should be chosen from quantitative validation, not runtime alone.

---

# 30. Why publish RH, e, r, and q together?

The four products are complementary.

### Relative humidity

Useful for:

- saturation state
- environmental diagnostics
- fog/cloud-related applications
- human-environment applications

### Specific humidity

Useful for:

- water-vapor transport
- atmospheric moisture budgets
- mass-based moisture analysis

### Mixing ratio

Useful for:

- thermodynamic calculations
- atmospheric moisture diagnostics

### Vapor pressure

Useful as:

- a physically interpretable intermediate
- an audit point in the transformation chain
- a useful variable in thermodynamic analyses

Publishing all four also reduces the need for downstream users to reconstruct nonlinear transformations from a single summarized quantity.

---

# 31. Performance philosophy

The project does **not** optimize for one benchmark number at the expense of scientific traceability.

The target is:

```text
correct
+
reproducible
+
restartable
+
memory-aware
+
diagnostic-rich
+
fast enough for production
```

The most important optimization is often not a micro-optimization in NumPy.

It is controlling the amount of data simultaneously resident in memory.

---

# 32. Operational tuning

A practical tuning sequence is:

### Step 1 — baseline

Use:

```text
N_SAMPLES=5000
CELL_CHUNK_SIZE=1024
SAMPLE_BATCH_SIZE=256
MAX_WORKERS=2
```

### Step 2 — observe

Record:

```text
peak RAM
average RAM
CPU utilization
I/O behavior
runtime per day
runtime per year
```

### Step 3 — tune one parameter at a time

For example:

```text
workers: 1 → 2 → 3
chunk:   512 → 1024 → 2048
batch:   128 → 256 → 512
```

### Step 4 — verify scientific equivalence

Never accept a performance change merely because it is faster.

Re-run representative cells/days and verify the numerical results and diagnostics.

---

# 33. Data integrity and provenance

For publication-grade output, archive:

```text
[ ] exact ERA5-Land dataset/version
[ ] CDS request or download recipe
[ ] raw-data checksums
[ ] daily aggregation method
[ ] variable names
[ ] variable units
[ ] coordinate convention
[ ] spatial resolution
[ ] spatial extent
[ ] time period
[ ] Monte Carlo sample count
[ ] random seed
[ ] chunk size
[ ] batch size
[ ] worker count
[ ] schema version
[ ] configuration hash
[ ] source/script SHA-256
[ ] Python version
[ ] NumPy version
[ ] SciPy version
[ ] xarray version
[ ] netCDF4 version
[ ] diagnostic NetCDF
[ ] final output checksum
```

This turns the output from an isolated file into an auditable research artifact.

---

# 34. Security and repository hygiene

Never commit:

```text
CDS API keys
API tokens
passwords
private keys
credentials
raw restricted data
machine-specific secrets
```

Also avoid committing:

```text
large raw ERA5-Land archives
temporary checkpoint directories
generated NetCDF products
Python cache directories
`.pytest_cache`
machine-specific output paths
```

Use the repository's `.gitignore` and keep credentials outside version control.

---

# 35. Reproducible publication recipe

A strong publication package can contain:

```text
code/
data-access/
preprocessing/
configuration/
validation/
outputs/
diagnostics/
figures/
metadata/
```

At minimum, preserve:

```text
source code
configuration
environment
data provenance
checksums
validation results
final NetCDF
diagnostic NetCDF
```

If the raw data cannot be redistributed, publish the exact acquisition recipe and checksums instead.

---

# 36. Citation

When this software contributes to a scientific result, cite:

1. **ERA5-Land** and the appropriate dataset documentation;
2. **HumidClimatologyEngine** and its exact software version;
3. the statistical/thermodynamic methodology used in the study;
4. the preprocessing workflow.

The repository includes:

```text
CITATION.cff
```

Update its repository URL and release version before public publication if they are still placeholders.

Suggested methodological wording:

> Daily moisture climatologies were generated from ERA5-Land 2-m air temperature, 2-m dew-point temperature, and surface pressure using HumidClimatologyEngine. The joint atmospheric state was represented in the transformed space \((T,T_d,\ln P)\), propagated through a deterministic Monte Carlo workflow, and converted to relative humidity, vapor pressure, mixing ratio, and specific humidity. The climatological calendar used a composite February 28/29 day.

Adapt this wording to the exact version and configuration used in the paper.

---

# 37. Versioning

The repository currently identifies the production package as:

```text
0.6.1
```

The core engine identifies itself as a v6 implementation.

See:

```text
CHANGELOG.md
```

for release history.

For scientific publication, always record:

```text
software version
git commit
configuration hash
input-data version
output checksum
```

A software version alone is not sufficient to reproduce a configured scientific run.

---

# 38. Troubleshooting

## Problem: missing monthly input

Symptom:

```text
Missing months for YYYY
```

Check:

```text
file names
year
month
directory
```

The indexer expects exactly one recognizable file for each month.

---

## Problem: duplicate monthly input

Symptom:

```text
Duplicate month YYYY-MM
```

Remove or rename ambiguous duplicates before production.

Do not silently select one of multiple files.

---

## Problem: coordinate mismatch

Check:

```python
ds_t.time
ds_d.time
ds_p.time

ds_t.latitude
ds_d.latitude
ds_p.latitude

ds_t.longitude
ds_d.longitude
ds_p.longitude
```

The engine intentionally rejects incompatible grids rather than guessing how to align them.

---

## Problem: excessive RAM

Reduce, in roughly this order:

```text
MAX_WORKERS
SAMPLE_BATCH_SIZE
CELL_CHUNK_SIZE
```

Only then consider reducing `N_SAMPLES`, because sample count affects scientific precision.

---

## Problem: slow Monte Carlo stage

Check:

```text
CPU utilization
disk I/O
worker count
chunk size
batch size
N_SAMPLES
```

Do not increase workers blindly: on Windows, more processes can increase memory pressure and reduce total throughput.

---

## Problem: invalid covariance diagnostics

Do not suppress the warning merely to finish the run.

Inspect:

```text
min_eigenvalue
invalid_covariance_fraction
correlations
valid_observation_count
```

Then determine whether the issue is:

- insufficient observations
- numerical precision
- genuinely ill-conditioned covariance
- problematic input data
- a modeling limitation

---

# 39. What this project is — and is not

### It is

- a multivariate climatology engine;
- a Monte Carlo uncertainty-propagation engine;
- a physically explicit moisture transformation pipeline;
- a restartable long-running scientific workflow;
- a diagnostic-rich NetCDF producer;
- a reproducibility-oriented research codebase.

### It is not

- a substitute for ERA5-Land documentation;
- a guarantee that the Gaussian model is physically exact;
- a black-box humidity calculator;
- a replacement for independent validation;
- a universal solution for extreme-tail modeling.

---

# 40. Future development roadmap

The current single-engine implementation can evolve toward a more modular architecture:

```text
src/humidclimatology/
├── calendar.py
├── physics.py
├── statistics.py
├── covariance.py
├── monte_carlo.py
├── checkpoints.py
├── io.py
├── validation.py
├── diagnostics.py
└── cli.py
```

Potential future capabilities include:

- formal command-line configuration;
- schema-validated YAML configuration;
- richer dataset adapters;
- Dask-aware execution;
- cloud/object-storage workflows;
- richer NetCDF encoding/compression controls;
- independent reference implementations;
- automated convergence reports;
- copula/mixture alternatives;
- extreme-tail validation;
- provenance manifests;
- reproducibility bundles;
- benchmark reports;
- CI matrices across supported Python versions.

These are future directions, not claims about the current release.

---

# 41. Project quality checklist

Before calling a run **production-ready**, verify:

- [ ] source compiles;
- [ ] repository tests pass;
- [ ] input units are documented;
- [ ] input coordinates are validated;
- [ ] all months are present;
- [ ] duplicate files are absent;
- [ ] leap-day behavior is confirmed;
- [ ] a scalar physics reference agrees with vectorized calculations;
- [ ] Welford statistics agree with an independent reference;
- [ ] Monte Carlo convergence is quantified;
- [ ] covariance diagnostics are acceptable;
- [ ] checkpoint hashes validate;
- [ ] restart behavior has been tested;
- [ ] final NetCDF dimensions are correct;
- [ ] final variables are finite where expected;
- [ ] physical bounds are checked;
- [ ] diagnostics have been inspected;
- [ ] configuration hash is archived;
- [ ] source hash/commit is archived;
- [ ] environment versions are archived;
- [ ] output checksum is archived.

---

# 42. Minimal quick start

```bash
# 1. Clone
git clone <your-repository-url>
cd HumidClimatologyEngine

# 2. Create an environment
python -m venv .venv

# Windows
.venv\Scripts\activate

# Linux/macOS
source .venv/bin/activate

# 3. Install
python -m pip install --upgrade pip
pip install -e ".[all]"

# 4. Run tests
pytest

# 5. Run syntax validation
python -m compileall src

# 6. Configure paths and production parameters

# 7. Run
python moisture_climatology.py
```

---

# 43. Output contract at a glance

```text
INPUT
  ERA5-Land daily T2m + Td2m + surface pressure

MODEL
  X = (T, Td, lnP)
  multivariate normal
  PSD / Cholesky validation

SAMPLING
  deterministic seed
  chunked cells
  batched samples

TRANSFORMATION
  RH
  vapor pressure
  mixing ratio
  specific humidity

STATISTICS
  mean
  std
  bias-corrected skewness
  Fisher excess kurtosis

UNCERTAINTY
  Monte Carlo standard errors
  valid sample counts

RESILIENCE
  annual checkpoints
  daily checkpoints
  SHA-256
  config hash
  atomic writes
  restart

OUTPUT
  main NetCDF
  diagnostic NetCDF
  provenance metadata
```

---

# 44. Final perspective

HumidClimatologyEngine is best understood as a **scientific computation pipeline**, not merely a humidity formula.

Its central strength is the complete chain:

```text
auditable input
      +
paired multivariate statistics
      +
dependence preservation
      +
numerically stable accumulation
      +
physical transformation
      +
streaming Monte Carlo
      +
higher-order distributional statistics
      +
uncertainty diagnostics
      +
atomic checkpoints
      +
restartability
      +
independent validation
      +
reproducible provenance
```

That chain is what makes a long-running climatology easier to:

- inspect,
- reproduce,
- debug,
- validate,
- publish,
- extend,
- and defend scientifically.

---

## License

MIT. See [`LICENSE`](LICENSE).

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Security

See [`SECURITY.md`](SECURITY.md).

## Methodology

See [`docs/methodology.md`](docs/methodology.md).

## Production workflow

See [`examples/production_run.md`](examples/production_run.md).

## Tutorial

See [`examples/HumidClimatologyEngine_tutorial.ipynb`](examples/HumidClimatologyEngine_tutorial.ipynb).

---

## Maintainer

**Amin Fazlkazemi**

Before public release, replace any remaining repository URL placeholders such as:

```text
https://github.com/YOUR-USERNAME/HumidClimatologyEngine
```

with the final canonical repository URL.
