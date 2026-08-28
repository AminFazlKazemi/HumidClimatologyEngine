# HumidClimatologyEngine

> **ERA5-Land hourly empirical moisture climatology, day-resolved probability modelling, and bivariate dependence analysis**

![HumidClimatologyEngine v7.5 — hourly empirical moisture climatology and bivariate probability modelling](V7_5.png)

[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![ERA5-Land](https://img.shields.io/badge/data-ERA5--Land-1F4E79)](https://cds.climate.copernicus.eu/)
[![NetCDF4](https://img.shields.io/badge/output-NetCDF4-4B8BBE)](https://www.unidata.ucar.edu/software/netcdf/)

## Contents

- [1. Scope and release status](#1-scope-and-release-status)
- [2. v6 versus v7](#2-v6-versus-v7)
- [3. Scientific data flow](#3-scientific-data-flow)
- [4. ERA5-Land input contract](#4-era5-land-input-contract)
- [5. Time indexing and the 5-day centred window](#5-time-indexing-and-the-5-day-centred-window)
- [6. Climatological calendar](#6-climatological-calendar)
- [7. Thermodynamic calculations](#7-thermodynamic-calculations)
- [8. Empirical moment engine](#8-empirical-moment-engine)
- [9. Two-dimensional probability framework](#9-two-dimensional-probability-framework)
- [10. Candidate univariate distributions](#10-candidate-univariate-distributions)
- [11. Bimodal Normal: five-parameter model](#11-bimodal-normal-five-parameter-model)
- [12. Copula and dependence layer](#12-copula-and-dependence-layer)
- [13. Model fitting and selection](#13-model-fitting-and-selection)
- [14. Spatial chunking and station/grid-cell extraction](#14-spatial-chunking-and-stationgrid-cell-extraction)
- [15. Checkpoints and power-failure recovery](#15-checkpoints-and-power-failure-recovery)
- [16. Progress reporting and ETA](#16-progress-reporting-and-eta)
- [17. Production outputs](#17-production-outputs)
- [18. Bivariate dominance report](#18-bivariate-dominance-report)
- [19. Provenance and reproducibility](#19-provenance-and-reproducibility)
- [20. Validation and QA](#20-validation-and-qa)
- [21. Performance and memory](#21-performance-and-memory)
- [22. Known limitations and explicit non-claims](#22-known-limitations-and-explicit-non-claims)
- [23. Repository / release layout](#23-repository--release-layout)
- [24. Operational runbook](#24-operational-runbook)
- [25. Configuration reference](#25-configuration-reference)
- [26. Example commands](#26-example-commands)
- [27. Scientific interpretation](#27-scientific-interpretation)
- [28. Release checklist](#28-release-checklist)
- [Appendix A — Formula sheet](#appendix-a--formula-sheet)
- [Appendix B — Output schema](#appendix-b--output-schema)
- [Appendix C — Failure classes](#appendix-c--failure-classes)

---

## 1. Scope and release status

HumidClimatologyEngine is a research-grade workflow for constructing **day-resolved moisture climatology from ERA5-Land hourly fields**.

The principal production period is:

```text
Target climatology: 1981-01-01 through 2020-12-31
Climatological slots: 366
Primary frequency: hourly
Primary variables: t2m, d2m, sp
```

The production philosophy is **empirical first**:

```text
ERA5-Land hourly T2m + D2m + SP
        |
        v
exact timestamp alignment + physical validity
        |
        v
RH / e / r / q for every valid hourly state
        |
        +----------------------+
        |                      |
        v                      v
empirical moments       empirical 2-D probability
(DOY x cell)            (DOY x cell x pair)
        |
        v
optional 5-day centred distribution-model layer
        |
        v
Normal / Skew-Normal / Pearson III / Beta /
Bimodal Normal / copula candidates
        |
        v
selection + diagnostics + visual dominance report
```

### Current implementation status

| Capability | v7.5 code status | Scientific role |
|---|---|---|
| Hourly empirical RH/e/r/q moments | Integrated | Primary production product |
| Annual disk-backed checkpoints | Integrated | Restart / power-failure protection |
| Per-DOY spatial completion bitmap | Integrated | Transaction-like progress state |
| Global progress / remaining units / ETA | Integrated | Operations |
| Empirical 2-D histogram/PDF | Integrated | Primary joint reference |
| Bivariate Gaussian evaluator | Integrated | Reference candidate only |
| 5-day centred extraction | Implemented as fitting/query layer | Distribution fitting |
| Normal / Skew-Normal / Pearson III / Beta fitting | Implemented in fitting layer | Candidate marginals |
| Bimodal Normal 5-parameter fitting | Implemented in fitting layer | Candidate multimodal model |
| Gaussian copula estimator | Implemented in fitting layer | Dependence candidate |
| Bivariate dominance report generator | Separate tool | Visual QA / reporting |
| Full automatic 40-year per-cell model selection in `main()` | **Not claimed** | Requires explicit second-pass orchestration |

The last row is deliberate. A production README must describe what the code actually executes, not what a future orchestrator is expected to execute.

---

## 2. v6 versus v7

### v6: historical / educational baseline

`moisture_climatology_v6.py` is intentionally retained. Its input contract is based on **daily statistics**, followed by a joint statistical representation of `(T, Td, ln P)` and Monte-Carlo propagation.

Conceptually:

```text
hourly source
   |
   v
DAILY AGGREGATION
   |
   v
statistics of (T, Td, ln P)
   |
   v
joint Gaussian approximation
   |
   v
Monte Carlo
   |
   v
RH / e / r / q
```

This is valuable for teaching and benchmarking:

- paired multivariate statistics;
- covariance construction;
- positive-semidefinite checks;
- Cholesky factorization;
- stochastic propagation;
- Monte-Carlo sampling error;
- restartable statistical states.

It is **not** the preferred production method for direct hourly moisture climatology because aggregation occurs before the nonlinear moisture transformations.

### v7: production direction

v7 starts from hourly paired observations and transforms the physical state before estimating the moisture distributions:

```text
hourly T2m + D2m + SP
        |
        v
thermodynamics
        |
        v
RH / e / r / q
        |
        v
empirical statistics + empirical joint probability
```

This preserves much more of the observed hourly structure and avoids making a Gaussian reconstruction the primary definition of the moisture distribution.

### Why v6 is not deleted

v6 remains part of the scientific history of the project. It should be cited as:

> **Historical / educational daily-statistical Gaussian-Monte-Carlo baseline.**

The project should not silently rewrite history by deleting it.

---

## 3. Scientific data flow

The production pipeline has two distinct layers.

### Layer A — primary climatology

```text
files
  -> timestamp validation
  -> grid validation
  -> unit normalization
  -> climatological-day mapping
  -> hourly physical transformation
  -> Welford/Pébay accumulation
  -> annual checkpoint
  -> multi-year merge
  -> NetCDF climatology
```

### Layer B — probability modelling

```text
5-day centred hourly sample
  -> valid paired sample
  -> candidate marginal fits
  -> candidate dependence fits
  -> diagnostics
  -> model ranking
  -> frozen selection result
  -> bivariate dominance report
```

The layers are intentionally separated so that the empirical climatology remains valid even if a particular parametric distribution later proves inadequate.

---

## 4. ERA5-Land input contract

### Required variables

| Variable | ERA5-Land field | Internal meaning |
|---|---|---|
| 2 m air temperature | `t2m` | T |
| 2 m dew-point temperature | `d2m` | Td |
| Surface pressure | `sp` | P |

Typical ERA5-Land source units are Kelvin for `t2m`/`d2m` and Pascal for `sp`; the engine converts them to degrees Celsius and hPa before thermodynamic calculations.

### Timestamp is authoritative

The loader must use the **actual datetime coordinate inside the NetCDF**, not the filename, not an assumed hour sequence, and not an inferred DOY.

The current input layer accepts `valid_time` first and can fall back to `time`. This matters because edge files can have a structure such as:

```text
Dimensions:
    valid_time = 48
    latitude   = 301
    longitude  = 301

Coordinates:
    valid_time = datetime64[ns]
    latitude   = 50.0 ... 20.0
    longitude  = 35.0 ... 65.0

Data variables:
    d2m(valid_time, latitude, longitude)
    t2m(valid_time, latitude, longitude)
```

### Grid contract

Before scientific processing, the three variables must be checked for:

- latitude existence and ordering;
- longitude existence and ordering;
- matching coordinate values;
- matching dimensions;
- matching timestamps;
- unit metadata;
- finite / missing-value semantics.

A mismatch must fail closed rather than silently regrid or guess.

---

## 5. Time indexing and the 5-day centred window

### Window definition

For a target date `D`, the fitting window is:

```text
D-2, D-1, D, D+1, D+2
```

Using hourly data, a complete five-day window spans:

```text
D-2 00:00  ->  D+2 23:00
```

which corresponds to **120 hourly timestamps** before missing-value filtering.

The target day remains the statistical label. The surrounding four days are padding used to stabilise the local distribution fit.

### Why the window is applied before fitting

The window must be extracted from the raw hourly series. Applying a five-day window after a DOY summary would destroy the temporal information needed by Pearson III, Skew-Normal, Bimodal Normal, Beta, and copula fitting.

Correct:

```text
raw hourly series
   -> centred 5-day extraction
   -> distribution fit
```

Incorrect:

```text
hourly
   -> daily/DOY summary
   -> 5-day averaging of summaries
   -> distribution fit
```

### Boundary padding

The target period is 1981–2020.

Additional files outside the target period are used only as temporal padding.

For the beginning of the record, the supplied edge file is:

```text
K:\kazemi\papers\temperature_interpolation\19801230-19801231.nc
```

This file supplies 1980-12-30 and 1980-12-31 for the first complete 1981 windows.

At the other end, data are available through 2021-06, which is more than sufficient to complete the final 2020 centred windows.

Padding observations are **not** counted as climatology years outside 1981–2020. They exist only to complete local fitting windows.

### Window completeness metadata

Every fitting query should retain:

```text
window_days_requested
window_days_available
window_completeness_fraction
paired_observation_count
first_timestamp
last_timestamp
```

A complete window is expected to have five calendar days and up to 120 hourly timestamps per variable before validity filtering.

---

## 6. Climatological calendar

The production calendar has 366 slots.

| Slot | Meaning |
|---:|---|
| 1–58 | January 1 through February 27 |
| 59 | Reserved |
| 60 | February 28 + February 29 composite |
| 61–366 | March 1 through December 31 |

Operational mapping:

```text
Leap year:
  Feb 28 -> 60
  Feb 29 -> 60
  Mar 01 -> 61

Non-leap year:
  Feb 28 -> 60
  Mar 01 -> 61
```

Slot 59 is reserved. The February 28/29 pooling is a formal model decision and must be preserved by downstream tools.

---

## 7. Thermodynamic calculations

The engine works internally with:

```text
T_C  = T_K  - 273.15
Td_C = Td_K - 273.15
P_hPa = P_Pa / 100
```

### Saturation vapor pressure

For `T >= 0 °C`:

```text
es(T) = 6.112 * exp(17.67*T / (T + 243.5))
```

For `T < 0 °C`:

```text
es(T) = 6.112 * exp(22.46*T / (T + 272.62))
```

### Vapor pressure

```text
e = es(Td)
```

### Relative humidity

```text
RH_raw = 100 * es(Td) / es(T)
```

The reported RH is bounded to `[0, 100]`. Supersaturation is counted separately rather than being hidden.

### Mixing ratio

```text
r = 0.622 * e / (P - e)
```

Only physically valid states with:

```text
e > 0
P > 0
e < P
```

are allowed into the `r` / `q` empirical sample.

### Specific humidity

```text
q = r / (1 + r)
```

The project keeps physical validity filtering separate from display clipping. For example, clipping RH to 100% does not make an invalid `e/P` state physically valid.

---

## 8. Empirical moment engine

For each climatological day and spatial cell, the primary production product accumulates:

```text
n
mean
M2
M3
M4
```

for:

```text
RH
e
r
q
```

The update strategy is Welford/Pébay-style online accumulation. It is numerically stable, mergeable, and does not require retaining the entire multidecadal sample history.

### Final estimators

For `n >= 2`:

```text
sample variance = M2 / (n - 1)
sample std      = sqrt(sample variance)
```

The production finalization additionally derives:

- mean;
- sample standard deviation;
- bias-corrected skewness;
- Fisher excess kurtosis.

### Mergeability

Annual states can be merged into the 1981–2020 multi-year state without replaying the original hourly archive.

This is central to both parallel processing and restartability.

---

## 9. Two-dimensional probability framework

The bivariate layer is deliberately **not locked to a bivariate normal distribution**.

### Primary joint product: empirical 2-D probability

The primary joint distribution is built directly from hourly paired observations.

For the current production configuration:

```text
pair = (RH, q)
bins = 8 x 8
```

The stored grid counts are normalised by sample count and bin area to produce a piecewise-constant density.

Therefore the primary object is:

```text
p_empirical(x, y | DOY, grid cell)
```

rather than:

```text
p_gaussian(x, y | DOY, grid cell)
```

### Reference Gaussian evaluator

A vectorized bivariate Gaussian PDF evaluator exists because it is useful as a comparison baseline:

```text
bivariate_gaussian_pdf(
    x, y,
    mean_x, std_x,
    mean_y, std_y,
    rho
)
```

This evaluator is **not** the definition of the empirical joint distribution.

### Why empirical first

The actual moisture distribution may be:

- skewed;
- heavy-tailed;
- bounded;
- multimodal;
- regime-dependent;
- asymmetric in its dependence structure.

A single Gaussian surface cannot be assumed to reproduce all of those features.

---

## 10. Candidate univariate distributions

The 5-day fitting layer evaluates candidate distributions per target date and local series.

### Normal

Useful as a transparent baseline and as a reference for improvement metrics.

### Skew-Normal

Useful when a single dominant regime has asymmetric shape.

The implementation retains two auditable variants:

1. **ClimateProcessingEngine-style moment parameterisation** using sample mean, standard deviation and sample skewness;
2. **SciPy maximum-likelihood fit**.

Both can be compared by the same likelihood information criteria.

### Pearson Type III

Pearson III is included because it can reproduce skewed continuous distributions and performed well for temperature in the companion `ClimateProcessingEngine` workflow.

It should be evaluated rather than assumed optimal for moisture variables.

### Beta

Beta is a candidate only when the variable is naturally bounded on `[0, 1]`.

Examples:

```text
RH_fraction = RH / 100
q           in [0, 1]
```

Endpoint handling requires care because ordinary Beta support is open at the boundaries. The current fitting layer uses a small epsilon transformation before fitting.

### Bimodal Normal

The two-component Gaussian mixture is included as a five-parameter model and is documented separately below.

### Selection principle

No family is declared universally correct.

The fitting layer returns the full candidate table so that the winning family and the alternatives remain auditable.

---

## 11. Bimodal Normal: five-parameter model

The project uses a **two-component Gaussian mixture**, not a two-piece normal distribution.

The five independent parameters are:

```text
w1
mu1
sigma1
mu2
sigma2
```

with:

```text
w2 = 1 - w1
```

### Fitting method

The current implementation uses an EM-style `GaussianMixture` fit with:

```text
n_components = 2
n_init       = 10
max_iter     = 1000
tol          = 1e-4
reg_covar    = 1e-6
random_state = 20260821
```

Components are sorted by mean after fitting so that:

```text
mu1 <= mu2
```

This makes stored parameters stable and interpretable.

### Retained diagnostics

The fit records:

- log-likelihood;
- AIC;
- AICc;
- BIC;
- component weights;
- component means;
- component standard deviations;
- separation metric comparable to Ashman-type separation;
- overlap coefficient;
- EM configuration.

### Scientific use

Bimodal Normal is valuable when a local distribution contains two distinct regimes. It must not be used merely because it has five parameters and a lower raw likelihood; model complexity is explicitly penalised by AICc/BIC and should also be reviewed with multimodality diagnostics.

---

## 12. Copula and dependence layer

The joint distribution is separated into two conceptual parts:

```text
marginal distributions
        +
dependence model
        =
joint distribution
```

### Pseudo-observations

The fitting layer converts each marginal to rank-based pseudo-observations:

```text
u = (rank(x) - 0.5) / n
v = (rank(y) - 0.5) / n
```

The Gaussian copula transform then uses the inverse standard normal transformation:

```text
z_x = Phi^-1(u)
z_y = Phi^-1(v)
```

and estimates dependence from the resulting transformed ranks.

### Current copula status

The current fitting layer contains a Gaussian copula estimator. It is a **candidate dependence model**, not a scientific axiom.

Future candidate families may include:

- t-copula;
- tail-dependent Archimedean copulas;
- empirical copula;
- other copulas selected by diagnostic performance.

### Important distinction

The existence of a copula fit does **not** mean the marginal distributions are Gaussian. For example:

```text
RH   -> Beta
q    -> Beta
dependence -> Gaussian copula
```

is a legitimate joint model.

Likewise:

```text
RH   -> Bimodal Normal
q    -> Pearson III
dependence -> copula
```

can be evaluated if scientifically justified.

---

## 13. Model fitting and selection

### Minimum sample size

The fitting layer uses:

```text
FIT_MIN_OBS = 30
```

for the candidate distribution fits.

This is an estimability guard, not a universal statistical truth. A future publication should report sensitivity to the minimum-sample threshold.

### Information criteria

For `k` fitted parameters and `n` observations:

```text
AIC  = 2k - 2 logL
BIC  = k log(n) - 2 logL
AICc = AIC + 2k(k+1) / (n-k-1)
```

The current selection order is:

```text
minimum AICc
        -> tie-break by BIC
```

The complete candidate table remains available for audit.

### What should be added before publication-scale model selection

A production release should augment the information criteria with:

- independent goodness-of-fit diagnostics;
- tail diagnostics;
- parameter stability checks;
- multimodality diagnostics;
- sensitivity to the 5-day window;
- sensitivity to minimum sample count.

That prevents a purely information-criterion-driven choice from being mistaken for proof of physical adequacy.

---

## 14. Spatial chunking and station/grid-cell extraction

The production accumulation engine uses spatial chunks:

```text
CHUNK_LAT = 32
CHUNK_LON = 64
```

The processing hierarchy is therefore approximately:

```text
year
  -> month
     -> climatological DOY slot
        -> latitude chunk
           -> longitude chunk
              -> hourly slices
                 -> moisture transformation
                 -> online accumulation
```

### Why this matters

Loading an entire regional 40-year hourly archive into memory is unnecessary and unsafe.

The chunk strategy bounds the number of grid cells resident in the working arrays and permits deterministic checkpointing at a much finer granularity.

### Station / grid-cell query

The 5-day fitting layer is designed to extract only the requested grid point or station-nearest cell and its required time window.

The extraction logic:

```text
requested date + coordinate
        |
        v
identify spatial block / cell
        |
        v
read only intersecting input files
        |
        v
use actual time coordinate
        |
        v
select D-2 ... D+2
        |
        v
align T2m / D2m / SP by exact timestamp
```

This avoids reopening the complete regional grid for every local model fit.

### Edge-file handling

The special 1980-12-30/31 combined T2m+D2m file is explicitly supported by the fitting/query layer because it is the supplied source for the initial 1981 centered windows.

---

## 15. Checkpoints and power-failure recovery

Power-loss recovery is a **scientific requirement**, not merely a convenience.

### Checkpoint unit

The annual checkpoint stores a completion flag for every:

```text
DOY × latitude-chunk × longitude-chunk
```

This is substantially finer than a year-level flag.

### Transaction order

For every spatial chunk:

```text
1. read hourly slices
2. calculate moisture
3. update sufficient statistics
4. write chunk into checkpoint
5. ds.sync()
6. mark completed_chunk = 1
7. ds.sync()
8. update progress JSON atomically
```

The completion flag is therefore written **after** the statistical state has been persisted.

### What happens after a power failure

If power is lost after a chunk is written but before its completion flag is committed, the chunk is recomputed.

If the completion flag is already committed, the chunk is skipped on restart.

Thus a failure does not invalidate the whole year.

### Checkpoint compatibility

A checkpoint is reusable only when its:

- schema version;
- configuration hash;
- grid shape;
- chunk shape;
- variable contract

match the current run.

### Progress state

A companion JSON record stores:

```text
completed_units
total_units
remaining_units
progress_percent
last_completed_doy
updated_utc
configuration hash
```

This JSON is written atomically.

---

## 16. Progress reporting and ETA

The engine explicitly reports progress at both the **year level** and **global run level**.

Example pattern:

```text
PROGRESS | Year 1987 | DOY 173 chunks 31/40 | 77.50% | 31/40 units | remaining 9 | rate 0.84 units/s | ETA 0.00 h
```

Global messages report:

```text
GLOBAL PROGRESS | 64.28% | 12345/19200 units | remaining 6855 | active years 2
```

At the end of accumulation:

```text
ACCUMULATION COMPLETE | 100.00% | ... | remaining 0
```

The quantities are explicit:

- percentage complete;
- completed units;
- total units;
- remaining units;
- processing rate;
- ETA where it is estimable;
- active/completed years.

This is designed for long unattended runs and for transparent postmortem analysis after an interruption.

---

## 17. Production outputs

### Main climatology

```text
moisture_climatology_1981_2020_v7_5.nc
```

Conceptual dimensions:

```text
doy x latitude x longitude
```

Primary variables:

```text
mean_rh   std_rh   skew_rh   kurt_rh
mean_e    std_e    skew_e    kurt_e
mean_r    std_r    skew_r    kurt_r
mean_q    std_q    skew_q    kurt_q
```

### Diagnostics

```text
moisture_climatology_diagnostics_1981_2020_v7_5.nc
```

Important diagnostics include:

```text
n_obs
supersaturation_fraction
invalid_e_over_p_fraction
```

### Bivariate reference parameters

```text
moisture_climatology_bivariate_1981_2020_v7_5.nc
```

The current default pair is:

```text
(RH, q)
```

The stored reference state includes per-DOY/per-cell paired count, means, covariance, and correlation.

### Empirical bivariate PDFs

For fixed-range pairs, the engine can build a separate empirical 2-D product:

```text
moisture_bivariate_empirical_RH__q_1981_2020_v7_5.nc
```

The current default configuration uses an `8 x 8` histogram grid for this empirical surface.

### Run manifest

```text
moisture_climatology_run_manifest_v7_5.json
```

The manifest records configuration, source hash, output paths, and hashes of generated products where available.

---

## 18. Bivariate dominance report

The separate tool:

```text
bivariate_dominance_report.py
```

turns a frozen model-selection product into a visual report.

### Expected selection input

```text
best_model_code(doy, latitude, longitude)
```

plus a global mapping:

```json
{"1":"Empirical-2D","2":"Gaussian-Copula","3":"t-Copula","4":"Beta-Copula","5":"Bimodal-Copula"}
```

The reporting tool is intentionally distribution-agnostic. It does not refit models and does not assume Gaussian copulas.

### Report products

It generates:

1. model legend and provenance page;
2. monthly model-dominance shares;
3. DOY-by-model spatial occurrence heatmap;
4. monthly spatial maps of the winning model;
5. PNG files for publication / QA use;
6. a combined PDF report.

### Scientific questions answered

**Monthly dominance:**

> Which bivariate model families dominate during each month?

**DOY heatmap:**

> On which climatological days does each model dominate across the domain?

**Spatial maps:**

> Where does each bivariate model family dominate during a given month?

### Separation of responsibilities

```text
fitting engine
    |
    v
best_model_code
    |
    v
bivariate_dominance_report.py
    |
    +--> PDF
    +--> PNG maps
    +--> monthly shares
    +--> DOY heatmap
```

The report generator must remain deterministic and must not alter the model-selection result.

---

## 19. Provenance and reproducibility

A publication-grade run should archive at least:

```text
software filename/version
schema version
configuration hash
source SHA-256
input inventory
input SHA-256 where practical
calendar convention
window definition
model-fitting configuration
random seed(s)
Python version
NumPy version
SciPy version
scikit-learn version where used
xarray version
netCDF4 version
checkpoint hashes
final NetCDF hashes
validation log
report-generator version/hash
```

### Run manifest principle

The README describes defaults. The **executed run manifest is authoritative**.

A publication must never copy configuration values from documentation when the exact values can be recovered from the run manifest.

### Reproducibility claim

“Same seed” alone is not sufficient to claim bit-for-bit reproducibility across arbitrary parallel execution. Reduction order, library versions, and floating-point behavior must also be controlled or the tolerance must be explicitly declared.

---

## 20. Validation and QA

The validation hierarchy is:

### Level 1 — calendar

Test leap and non-leap mappings, including the reserved slot.

### Level 2 — scalar physics

Compare the implemented saturation and moisture formulae against an independent scalar reference.

### Level 3 — vectorized physics

Verify numerical equivalence between scalar and vector implementations within an explicit tolerance.

### Level 4 — moments

Compare Welford/Pébay online estimates with trusted NumPy/SciPy batch estimates.

### Level 5 — merge

Verify that:

```text
accumulate(A) + accumulate(B)
```

and:

```text
accumulate(A + B)
```

agree to the declared numerical tolerance.

### Level 6 — bivariate empirical PDF

Verify that:

```text
sum(bin_counts) = valid_pair_count
```

and that integrated piecewise density is one within numerical tolerance.

### Level 7 — distribution fitting

Use synthetic distributions with known shape, skewness, and multimodality. Confirm that the candidate fitting layer behaves as expected.

### Level 8 — restart

Run a small case, interrupt it, restart it, and compare against uninterrupted execution.

### Level 9 — serialization

Verify dimensions, coordinates, units, fill values, metadata, and checksums.

### Level 10 — report

Verify that the dominance report uses only the frozen selection array and that all model codes are mapped correctly.

---

## 21. Performance and memory

The primary optimisation goal is:

> **minimum wall-clock time subject to numerical correctness, memory stability, and reproducibility.**

### Current controls

```text
MAX_WORKERS = 2
CHUNK_LAT   = 32
CHUNK_LON   = 64
PROGRESS_FLUSH_CHUNKS = 16
PROGRESS_LOG_EVERY_CHUNKS = 8
```

### Tuning order

When memory is constrained:

1. reduce spatial chunk size;
2. reduce worker count;
3. verify hidden dtype copies;
4. inspect I/O and serialization;
5. only then consider changing other scientific settings.

Operational settings must never silently change the scientific model.

---

## 22. Known limitations and explicit non-claims

### The empirical histogram is resolution-dependent

An `8 x 8` grid is a practical reference surface, not a universal optimal bandwidth choice. Histogram sensitivity should be tested before publication.

### The bivariate Gaussian is not the truth

It is a reference evaluator only.

### Current copula layer is not a complete universal copula search

The current fitting layer contains a Gaussian copula estimator. A publication-scale comparison should add and benchmark other dependence families where tail dependence or asymmetry is scientifically important.

### Candidate selection is not the same as physical validation

A lower AICc or BIC does not prove physical truth. It identifies a better statistical compromise under the candidate set and sampling assumptions.

### Five-day window is a modelling choice

The centred `±2 day` window is intended to stabilise local distribution fitting while preserving seasonal locality. It should be sensitivity-tested against alternative window widths for publication-critical results.

### Endpoint behaviour of Beta is special

RH can contain exact 0 or 100% values. Standard Beta fitting requires careful endpoint handling; the current layer uses an epsilon transformation. A zero/one-inflated model is a future extension where endpoint mass is scientifically important.

### Padding is not climatology

1980 and 2021 data used for window completion are not additional climatology years.

### Full model-selection orchestration is a separate pass

The current v7.5 `main()` builds the empirical climatology and empirical bivariate products. The five-day distribution/coplanula layer is exposed as fitting functions and extraction utilities; it should be orchestrated as a controlled second pass before claiming a full 40-year, every-cell, every-DOY candidate-selection product.

---

## 23. Repository / release layout

A publication-oriented repository should contain:

```text
HumidClimatologyEngine/
|
+-- moisture_climatology_v7_5.py
+-- bivariate_dominance_report.py
+-- moisture_climatology_reset.py
+-- README.md
+-- HumidClimatologyEngine_v7_5_Detailed_Scientific_Engineering_Reference.pdf
+-- HumidClimatologyEngine_v6_vs_v7_5_Detailed_Comparison.pdf
+-- HumidClimatologyEngine_v7_5_Bivariate_Dominance_Report_Documentation.pdf
+|
+ +-- tests/
+ +-- configs/
+ +-- docs/
+ +-- reports/
+ +-- examples/
+ +-- CHANGELOG.md
+ +-- CITATION.cff
+ +-- LICENSE
+```
+
+The tree above is the recommended release layout. Files that are not actually shipped must not be represented as if they already exist in a package archive.

---

## 24. Operational runbook

### Before starting

- freeze the code version;
- freeze configuration;
- verify input inventories;
- verify timestamps and units;
- verify the common grid;
- confirm the 1980 edge file exists for the first windows;
- confirm adequate 2021 padding for the final windows.

### During accumulation

Monitor:

```text
GLOBAL PROGRESS
YEAR PROGRESS
completed units
remaining units
rate
ETA
RAM / CPU
checkpoint timestamp
```

### After interruption

1. leave the partial checkpoint files intact;
2. rerun the same code and configuration;
3. the completion bitmap determines which DOY/chunk units are skipped;
4. incomplete units are recomputed;
5. the year is accepted only when all spatial-day units are committed.

### Before publication

Archive:

```text
main NetCDF
statistics diagnostics
bivariate parameters
empirical 2-D PDFs
selection product
bivariate dominance PDF/PNGs
run manifest
source hash
input manifest
validation report
README
scientific PDF
```

---

## 25. Configuration reference

| Parameter | Current value | Role |
|---|---:|---|
| `START_YEAR` | 1981 | Target climatology start |
| `END_YEAR` | 2020 | Target climatology end |
| `DOY_COUNT` | 366 | Calendar slots |
| `MAX_WORKERS` | 2 | Process concurrency |
| `CHUNK_LAT` | 32 | Spatial latitude chunk |
| `CHUNK_LON` | 64 | Spatial longitude chunk |
| `PROGRESS_FLUSH_CHUNKS` | 16 | JSON progress flush cadence |
| `PROGRESS_LOG_EVERY_CHUNKS` | 8 | Progress log cadence |
| `BIVARIATE_PAIRS` | `(('rh','q'),)` | Primary bivariate pair |
| `BIVARIATE_NX` | 8 | Empirical x bins |
| `BIVARIATE_NY` | 8 | Empirical y bins |
| `BUILD_EMPIRICAL_BIVARIATE` | True | Build empirical 2-D PDF |
| `WINDOW_HALF_WIDTH_DAYS` | 2 | Centred local fit window |
| `WINDOW_SIZE_DAYS` | 5 | Window length |
| `FIT_MIN_OBS` | 30 | Candidate-fit minimum sample |
| `BIMODAL_N_INIT` | 10 | Two-component EM starts |
| `BIMODAL_MAX_ITER` | 1000 | EM iteration limit |
| `BIMODAL_TOL` | `1e-4` | EM convergence tolerance |
| `BIMODAL_REG_COVAR` | `1e-6` | EM covariance regularisation |

Treat the executable configuration and run manifest as authoritative for an actual release.

---

## 26. Example commands

### Run the main climatology

```bash
python moisture_climatology_v7_5.py
```

### Generate the bivariate dominance report

```bash
python bivariate_dominance_report.py \
  moisture_bivariate_model_selection_1981_2020_v7_5.nc \
  --output-pdf bivariate_distribution_dominance_report.pdf \
  --output-dir bivariate_dominance_figures
```

### Inspect a station/grid-cell window

The fitting layer is designed to expose a target date, coordinate/cell, and variable through the centred 5-day extraction utilities. For large-scale station workflows, wrap those functions in a station-ID table rather than reading the entire regional grid repeatedly.

### Reset generated artifacts

```bash
python moisture_climatology_reset.py --dry-run
```

Use a non-interactive deletion only after reviewing the dry-run output.

---

## 27. Scientific interpretation

The four primary moisture fields are not interchangeable.

### RH

Intuitive, but strongly nonlinear because both temperature and dew point enter the ratio.

### Vapor pressure

Closer to actual water-vapor partial pressure and less directly temperature-normalised than RH.

### Mixing ratio

Sensitive to the `P - e` denominator and therefore to physically invalid or near-invalid states.

### Specific humidity

Bounded and useful for many transport and budget applications.

The bivariate products answer a different question:

> **How are two moisture variables jointly distributed at a given climatological day and location?**

The purpose of keeping an empirical surface alongside parametric models is to distinguish:

```text
data-driven evidence
```

from:

```text
parametric approximation
```

That distinction should remain explicit in every publication and downstream analysis.

---

## 28. Release checklist

Before tagging a production release:

### Data

- [ ] All target months are present.
- [ ] Edge padding for 1981 is present.
- [ ] 2021 padding is sufficient for the 2020 endpoint.
- [ ] `valid_time`/`time` coordinates are valid and monotonic.
- [ ] T2m, D2m and SP grids agree.
- [ ] Units are verified.

### Scientific core

- [ ] Calendar tests pass.
- [ ] Physical reference tests pass.
- [ ] Online moments match offline references.
- [ ] Merge equivalence passes.
- [ ] Empirical 2-D PDF normalisation passes.
- [ ] Bimodal synthetic recovery passes.
- [ ] Skew-Normal candidate fit passes.
- [ ] Pearson III candidate fit passes.
- [ ] Beta endpoint handling passes.
- [ ] Copula estimator passes synthetic dependence checks.

### Operations

- [ ] Checkpoint/restart test passes.
- [ ] Completion bitmap is monotonic.
- [ ] Progress percentage is correct.
- [ ] Remaining-unit count is correct.
- [ ] ETA does not crash the run when rate is unavailable.
- [ ] RAM/CPU logging is non-fatal.

### Outputs

- [ ] Main NetCDF dimensions are correct.
- [ ] Diagnostic NetCDF dimensions are correct.
- [ ] Bivariate output is internally consistent.
- [ ] `_FillValue` metadata is correct.
- [ ] Units and long names are present.
- [ ] Run manifest is complete.
- [ ] SHA-256 hashes are archived.

### Documentation

- [ ] README matches the executable code.
- [ ] Scientific PDF matches the README.
- [ ] v6 historical status is preserved.
- [ ] v6 limitations are explicitly documented.
- [ ] Bivariate dominance report documentation is included.
- [ ] No future/unimplemented capability is presented as shipped production behavior.

---

# Appendix A — Formula sheet

```text
T_C   = T_K - 273.15
Td_C  = Td_K - 273.15
P_hPa = P_Pa / 100
```

```text
es_water(T) = 6.112 * exp(17.67*T / (T + 243.5))
```

```text
es_ice(T) = 6.112 * exp(22.46*T / (T + 272.62))
```

```text
e  = es(Td)
RH = 100 * e / es(T)
r  = 0.622 * e / (P - e)
q  = r / (1 + r)
```

```text
sample variance = M2 / (n - 1)
```

```text
AIC  = 2k - 2 logL
BIC  = k log(n) - 2 logL
AICc = AIC + 2k(k+1)/(n-k-1)
```

Bimodal Normal:

```text
f(x) = w1*N(x | mu1, sigma1)
     + (1-w1)*N(x | mu2, sigma2)
```

with stored parameters:

```text
w1, mu1, sigma1, mu2, sigma2
```

---

# Appendix B — Output schema

### Main climatology

```text
doy
latitude
longitude
month
day

n_obs
mean_rh
std_rh
skew_rh
kurt_rh

mean_e
std_e
skew_e
kurt_e

mean_r
std_r
skew_r
kurt_r

mean_q
std_q
skew_q
kurt_q
```

### Diagnostics

```text
n_obs
supersaturation_fraction
invalid_e_over_p_fraction
```

### Bivariate reference parameters

For each configured pair:

```text
pair_<x>__<y>_mean_x
pair_<x>__<y>_mean_y
pair_<x>__<y>_Cxy
```

plus the derived correlation fields in the final bivariate product where applicable.

### Empirical 2-D PDF

```text
x_edges
 y_edges
counts / density
valid_pair_count
DOY
latitude
longitude
```

The exact variable names in the final file are governed by the executable schema and should be checked from the generated NetCDF rather than inferred from this document.

---

# Appendix C — Failure classes

A production log should make the failure class identifiable.

### Input / provenance failure

Examples:

```text
missing month
inconsistent coordinate
unit mismatch
unexpected time coordinate
checksum mismatch
```

### Calendar failure

Examples:

```text
wrong leap-year mapping
reserved slot populated incorrectly
post-February shift
```

### Statistical failure

Examples:

```text
insufficient paired observations
non-finite moments
merge inconsistency
```

### Physical failure

Examples:

```text
e <= 0
P <= 0
e >= P
non-finite saturation pressure
```

### Checkpoint failure

Examples:

```text
schema mismatch
configuration hash mismatch
incomplete spatial-day bitmap
corrupt JSON progress record
```

### Model-fitting failure

Examples:

```text
insufficient sample
non-positive scale
failed likelihood evaluation
invalid Beta support
EM non-convergence
```

### Reporting failure

Examples:

```text
missing best_model_code
unknown model code
missing latitude/longitude
invalid model_names mapping
```

---

## Final scientific position

HumidClimatologyEngine deliberately maintains two methodological generations:

```text
v6
Historical / educational
Daily statistics -> Joint Gaussian -> Monte Carlo

v7
Production direction
Hourly data -> Direct thermodynamic transformation -> Empirical climatology
                                      |
                                      +--> empirical 2-D probability
                                      |
                                      +--> 5-day local distribution fitting
                                      |
                                      +--> candidate marginals + copula dependence
                                      |
                                      +--> Bimodal / skewed / bounded alternatives
```

The project should prefer **evidence over a fixed distributional assumption**. The empirical 2-D product is the non-parametric reference; Normal, Skew-Normal, Pearson III, Beta, Bimodal Normal, and copula families are modelling candidates whose adequacy must be demonstrated for the relevant day, location, and scientific question.

The README is intentionally explicit about what is implemented, what is a reference candidate, and what remains a second-pass orchestration task. That distinction is part of the scientific reproducibility contract.


## v8.0 FINAL — Single-Pass Empirical Production Engine

Version 8.0 introduces an engineering upgrade of the hourly empirical workflow.

### Main improvements

- Single-pass generation of the physical moisture variables and empirical products.
- The RH–q bivariate empirical probability mass function is designed to be accumulated during the same hourly processing stream.
- Avoids a second full ERA5-Land scan only for bivariate statistics.
- Checkpoint architecture is extended for long production runs and restart safety.
- Keeps the non-parametric empirical approach: no Gaussian assumption is imposed on joint moisture behaviour.

### Production philosophy

The engine separates:
- physical conversion (T, Td, pressure → moisture variables),
- online statistical accumulation,
- empirical probability products,
- final NetCDF publishing.

The objective is a reproducible climate-processing workflow suitable for multi-decadal hourly ERA5-Land datasets.

### Release note

v8.0 is primarily an internal architecture and performance release. It improves scalability and computational efficiency rather than changing the scientific definition of the climatology.
