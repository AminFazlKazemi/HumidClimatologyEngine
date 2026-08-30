# HumidClimatologyEngine

> **ERA5-Land hourly empirical moisture climatology and day-resolved bivariate probability framework**

[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB)](https://www.python.org/)
[![Data](https://img.shields.io/badge/data-ERA5--Land-1F4E79)](https://cds.climate.copernicus.eu/)
[![Output](https://img.shields.io/badge/output-NetCDF4-4B8BBE)](https://www.unidata.ucar.edu/software/netcdf/)

## 1. What this project is

HumidClimatologyEngine produces climatological distributions of atmospheric moisture variables from ERA5-Land 2-m air temperature, 2-m dew-point temperature, and surface pressure.

The current production branch is **v7.5 Hourly Empirical Production**. It is intentionally different from the earlier v6 method.

### v7.5 production method

```text
ERA5-Land hourly T2m + D2m + SP
        |
        v
paired hourly validation
        |
        v
thermodynamic transformation
        |
        +--> RH
        +--> vapor pressure e
        +--> mixing ratio r
        +--> specific humidity q
        |
        v
online empirical moments
(n, mean, M2, M3, M4)
        |
        +--> 366-slot daily climatology
        |
        +--> day/grid-cell bivariate probability parameters
```

The primary v7.5 moisture products are **empirical**. They do not use the v6 Gaussian/Monte-Carlo reconstruction for their primary statistics.

---

## Historical v6 status

`moisture_climatology_v6.py` is intentionally retained as the **historical/educational baseline**. Its workflow starts from **daily statistics**, models the daily state `(T, Td, ln P)` with a joint Gaussian approximation, and propagates that model by Monte Carlo. This is useful pedagogically for covariance, PSD, Cholesky, Monte Carlo, and uncertainty propagation, but it is **not the preferred production method for direct hourly moisture climatology** because the early daily aggregation removes within-day information needed by nonlinear humidity transformations.

The v7.5 production path instead reads hourly T2m, D2m, and surface pressure and transforms every valid hourly state directly.

## 2. v6 is retained intentionally

`moisture_climatology_v6.py` is **not deleted**.

It is retained as a **historical and educational baseline** because it demonstrates a useful methodological idea: represent the atmospheric state jointly as `(T, Td, ln P)`, estimate covariance, propagate a stochastic model, and then transform the sampled state into nonlinear moisture diagnostics.

However, v6 is not the preferred production method for direct hourly moisture climatology.

### The main v6 methodological limitation

v6 operates on **daily statistics**, not the original hourly paired atmospheric states. That means information about the within-day joint timing of T2m, D2m, and surface pressure has already been compressed before the nonlinear moisture transformation.

Its architecture is approximately:

```text
hourly ERA5-Land
      |
      v
 daily aggregation
      |
      v
 statistics of (T, Td, ln P)
      |
      v
 joint Gaussian model
      |
      v
 Monte Carlo
      |
      v
 RH / e / r / q
```

This is a valid **teaching and methodological baseline**, but it should not be mistaken for a direct empirical distribution of hourly moisture states.

### Other v6 limitations

* Gaussianity is a structural approximation and may be weak for multimodal regimes, tails, and regime transitions.
* Monte-Carlo sampling adds sampling error, especially for higher moments.
* Covariance quality controls the stochastic propagation because the dependence structure is estimated before the nonlinear transformation.
* A bivariate Gaussian derived from v6 covariance is a model-based PDF, not an empirical two-dimensional distribution.

The purpose of keeping v6 is therefore historical transparency, methodological comparison, and education.

---

## 3. Why v7.5 is the production branch

The central scientific improvement is the order of operations:

```text
hourly paired atmospheric state
          |
          v
nonlinear moisture calculation
          |
          v
empirical distributional accumulation
```

instead of calculating moisture after first compressing the atmospheric state into daily summaries.

For nonlinear functions such as RH, mixing ratio, and specific humidity, this preserves the actual hourly relationship between temperature, dew point, and pressure as much as the available hourly reanalysis permits.

---

## 4. Thermodynamic formulation

Let temperature and dew point be in degrees Celsius and pressure in hPa.

### Saturation vapor pressure

For T >= 0 C:

`es(T) = 6.112 * exp(17.67*T / (T + 243.5))`

For T < 0 C:

`es(T) = 6.112 * exp(22.46*T / (T + 272.62))`

The code uses water phase for T >= 0 C and ice phase for T < 0 C.

### Vapor pressure

`e = es(Td)`

### Relative humidity

`RH_raw = 100 * es(Td) / es(T)`

RH is clipped to the reporting interval 0-100 percent, while supersaturation is retained as a diagnostic.

### Mixing ratio

`r = 0.622 * e / (P - e)`

The state is accepted only when `e > 0`, `P > 0`, and `e < P`.

### Specific humidity

`q = r / (1 + r)`

---

## 5. Climatological calendar

The project uses a formal 366-slot calendar.

| Slot | Meaning |
|---:|---|
| 1-58 | Jan 1 through Feb 27 |
| 59 | Reserved |
| 60 | Feb 28 + Feb 29 composite |
| 61-366 | Mar 1 through Dec 31 |

Operationally:

```text
Leap year:
  Feb 28 -> 60
  Feb 29 -> 60
  Mar 01 -> 61

Non-leap year:
  Feb 28 -> 60
  Mar 01 -> 61
```

This is a scientific model contract, not merely an array-indexing convenience.

---

## 6. Empirical statistics

For each **DOY x grid cell** and each moisture variable, v7.5 maintains:

* `n`
* `mean`
* `M2`
* `M3`
* `M4`

using numerically stable online Welford/Pébay-style updates.

The final products are:

* mean
* sample standard deviation
* bias-corrected skewness
* Fisher excess kurtosis

Annual states are mergeable, so the 1981-2020 climatology is built from annual sufficient statistics rather than replaying raw multi-decadal history during the final merge.

---

## 7. Hourly memory strategy

v7.5 does not materialize a complete monthly `time x latitude x longitude` cube for all three variables.

The processing order is:

```text
time slice
   -> spatial block
      -> T / Td / P conversion
         -> RH / e / r / q
            -> online statistics
```

The default spatial block is controlled by the internal chunk dimensions. This bounds temporary arrays while keeping the annual checkpoint compact.

The annual accumulator is still proportional to `DOY_COUNT x grid cells`; for very large grids, the next scaling step is spatial sharding with the same merge algebra.

---

## 8. Two-dimensional probability, by day and grid cell

This is a first-class project requirement.

### What v7.5 stores

The default configured pairs are:

* `(RH, q)`
* `(RH, r)`

For every climatological day and every grid cell, v7.5 stores:

* paired-valid count `n`
* `mean_x`
* `mean_y`
* `std_x`
* `std_y`
* covariance `cov`
* Pearson correlation `corr`

The parameters are stored in a separate NetCDF file:

```text
moisture_climatology_bivariate_1981_2020_v7_5.nc
```

### Reference bivariate PDF

The code provides a vectorized evaluator:

`bivariate_gaussian_pdf(x, y, mean_x, std_x, mean_y, std_y, rho)`

It evaluates the standard bivariate Gaussian density using the empirical per-day/per-cell parameters.

This is a **reference probability model**, not a claim that the empirical joint distribution is Gaussian.

### Why the distinction matters

A smooth Gaussian reference PDF is useful for:

* joint-event calculations,
* probability contours,
* conditional calculations,
* threshold integration,
* comparisons between days and locations.

But the primary empirical hourly observations remain the scientific source of the marginal moments.

---

## 9. Beta distribution: candidate, not imposed assumption

Beta distributions are potentially useful for bounded variables such as:

* `RH / 100`
* `q`

but v7.5 does **not** silently force a Beta distribution onto them.

This is deliberate. RH can contain boundary behavior near 0 or 100 percent, and q can be strongly concentrated near zero. A plain Beta density on `(0,1)` is therefore not automatically adequate.

The intended future comparison framework is:

```text
Empirical 2-D reference
        |
        +-- KDE / adaptive histogram
        +-- copula-based model
        +-- Beta marginal candidates for bounded variables
        +-- Gaussian / t / alternative copulas
```

The chosen distribution should be justified by per-DOY/per-cell diagnostics rather than globally imposed.

---

## 10. Input contract and validation

Required variables:

| Quantity | Variable | Expected source unit |
|---|---|---|
| 2-m temperature | `t2m` | K or deg C |
| 2-m dew point | `d2m` | K or deg C |
| surface pressure | `sp` | Pa or hPa |

The implementation checks:

* required dimensions and variables,
* non-empty time axes,
* strictly increasing time,
* hourly 60-minute spacing,
* exact time-coordinate agreement between T2m, D2m, and SP,
* exact latitude-coordinate agreement,
* exact longitude-coordinate agreement,
* requested year/month consistency,
* supported source units,
* complete monthly file inventory.

The rule is **fail closed**: ambiguous inputs are rejected instead of guessed.

---

## 11. Diagnostics

The primary diagnostic product contains:

* `valid_observation_count`
* `supersaturation_fraction`
* `invalid_e_over_p_fraction`

These distinguish missing coverage from physically invalid states.

The bivariate product additionally carries the per-pair sample count and dependence parameters.

---

## 12. Checkpoints and restart

Annual checkpoint state is stored as NetCDF with a JSON sidecar.

The accepted checkpoint contract includes:

* year
* schema version
* configuration hash
* grid shape
* completed native dates
* UTC timestamp
* SHA-256 checksum

A checkpoint is reused only when its metadata and checksum are compatible with the active configuration.

The hourly state for a completed climatological slot is rebuilt deterministically if interruption occurs before the slot is marked complete.

---

## 13. Provenance

Every final run produces a manifest containing:

* implementation name
* schema version
* configuration hash
* script SHA-256
* period
* calendar policy
* input directories
* configured bivariate pairs
* SHA-256 hashes of final main, diagnostic, and bivariate outputs
* creation UTC timestamp

For a publication release, also archive the ERA5-Land acquisition identifiers, the environment specification, and the exact input inventory/checksums.

---

## 14. Output files

For the default 1981-2020 run:

```text
moisture_climatology_1981_2020_v7_5.nc
moisture_climatology_diagnostics_1981_2020_v7_5.nc
moisture_climatology_bivariate_1981_2020_v7_5.nc
moisture_climatology_run_manifest_v7_5.json
```

### Main product

For each of RH, e, r, and q:

```text
mean_<var>
std_<var>
skew_<var>
kurt_<var>
```

Dimensions:

```text
doy x latitude x longitude
```

### Bivariate product

For each configured pair, for example `rh__q`:

```text
rh__q_n
rh__q_mean_x
rh__q_mean_y
rh__q_std_x
rh__q_std_y
rh__q_cov
rh__q_corr
```

---

## 15. QA and tests

v7.5 includes tests for:

* leap-day calendar mapping,
* fourth-order online moment accumulation,
* covariance accumulation,
* Pébay merge equivalence,
* thermodynamic physical bounds,
* bivariate Gaussian PDF evaluation.

The current pure scientific/unit tests pass without requiring the NetCDF runtime library. The full production run additionally requires `netCDF4`.

### Required publication QA

Before a final scientific release, also run:

1. online-vs-offline statistics on a trusted synthetic dataset;
2. online-vs-offline statistics on representative real ERA5-Land subsets;
3. interrupted-vs-uninterrupted restart equivalence;
4. independent scalar thermodynamic reference checks;
5. numerical stress tests for near-zero variance and extreme tails;
6. bivariate PDF comparison against direct empirical samples;
7. spatial-shard merge equivalence before very large production runs.

A successful process exit is not, by itself, a scientific acceptance criterion.

---

## 16. Scientific limitations

### v7.5 empirical marginals

The empirical approach removes the v6 Gaussian/Monte-Carlo approximation for the primary marginal moisture products, but it does not remove uncertainty from ERA5-Land itself, missing-data processes, temporal dependence, spatial representativeness, or thermodynamic formulation choices.

### Bivariate reference model

The current bivariate PDF is Gaussian. It is intentionally labeled a reference model because the true empirical joint distribution may be skewed, heavy-tailed, multimodal, or tail-dependent.

### Calendar convention

Pooling Feb 28 and Feb 29 is a deliberate modeling choice and should be reported in publications.

### Higher moments

Skewness and especially kurtosis are substantially more sensitive to sample size and tail behavior than the mean. Stability must be demonstrated rather than inferred visually.

---

## 17. Recommended scientific roadmap

```text
v6
Historical / educational daily-statistical Gaussian-MC baseline

        |
        v

v7.5
Hourly empirical production

        |
        +--> empirical marginal moments
        |
        +--> day-resolved bivariate reference PDFs
        |
        v

v7.5+
Empirical 2-D KDE / copula / Beta-marginal candidate models
        |
        v
per-DOY/per-cell model selection and tail diagnostics
```

The long-term goal is not to force one probability family everywhere. It is to retain the observed day-specific dependence structure and choose a compact probabilistic representation only when it demonstrably describes the empirical sample.

---

## 18. Configuration defaults

| Parameter | Default |
|---|---|
| `START_YEAR` | 1981 |
| `END_YEAR` | 2020 |
| `DOY_COUNT` | 366 |
| `MAX_WORKERS` | 2 |
| `CHUNK_LAT / CHUNK_LON` | 8192 |
| `BIVARIATE_PAIRS` | `(rh, q), (rh, r)` |
| `SCHEMA_VERSION` | `7.3` |
| `CHECKPOINT_VERSION` | `7.3` |

The executed configuration is authoritative. Scientific defaults must never be changed silently.

---

## 19. Installation and execution

The production runtime requires at minimum the packages used by the engine, including:

```text
numpy
xarray
netCDF4
tqdm
```

Then run:

```bash
python moisture_climatology_v7_5.py
```

For a production workstation, first run the built-in tests and a small real-data pilot before launching 1981-2020.

---

## 20. Release principle

A release should be treated as one scientific unit:

```text
source code
+ README
+ scientific PDF
+ tests
+ configuration contract
+ provenance manifest
+ output schema
+ changelog
```

When the formula, calendar policy, variable set, probability model, or output schema changes, the README and PDF must be regenerated in the same release.

---

## 21. License

MIT, when the complete public repository is assembled.## 12. Explicit progress and remaining-work accounting

v7.5 reports progress in calculation units rather than only files or years.

```text
total_units = years × 365 valid climatological slots/year × spatial_chunk_count
completed_units = committed checkpoint tiles
remaining_units = total_units - completed_units
percent = 100 × completed_units / total_units
ETA = remaining_units / current committed-unit rate
```

The runtime reports year-level and global progress including percent complete, completed/total units, remaining units, processing rate, and ETA. A typical message is:

```text
GLOBAL PROGRESS | 37.42% | 2,184,320/5,838,080 units | remaining 3,653,760 | active years 2
```

This progress is backed by lightweight JSON metadata, while the NetCDF completion bitmap remains the authoritative restart state.



## Empirical two-dimensional probability distribution

v7.5 does not assume that the joint distribution is bivariate normal. The primary two-dimensional product is an **empirical, piecewise-constant joint PDF** computed directly from hourly paired observations for each climatological DOY and grid cell.

For the production configuration the supported empirical pair is:

```text
(RH, q)
```

with physical ranges:

```text
RH : 0 ... 100 %
q  : 0 ... 1
```

The probability mass in each 2-D bin is observed directly:

```text
P(bin i,j) = count(i,j) / N_valid
```

and the corresponding piecewise-constant density is:

```text
f(x,y) = count(i,j) / (N_valid * Delta_x * Delta_y)
```

where `(x,y)` lies inside bin `(i,j)`. This is therefore an empirical joint distribution rather than a fitted Gaussian surface. The exact hourly counts are stored as integer histogram counts so the probabilities remain auditable.

The empirical PDF pass has its own transactional restart state. A small `next_year[doy, y_chunk, x_chunk]` bitmap is written in the same NetCDF transaction as the histogram counts. A power failure therefore cannot permanently record a count update without the corresponding progress state.

### Parametric distributions are candidates, not assumptions

The code also retains a bivariate Gaussian evaluator for comparison and future model-selection work. It is explicitly labeled as a reference candidate. For bounded marginals such as RH and q, Beta distributions can be evaluated as marginal candidates, followed by a copula model for dependence. Gaussian copula, t-copula, or other copula families should be selected only after goodness-of-fit and tail diagnostics.

The scientific rule is:

> **Empirical joint PDF first; parametric PDF only when validated against the data.**

The empirical 2-D output is generated in a second pass after the annual hourly checkpoints are complete. This keeps the annual restart checkpoints compact while giving the bivariate layer its own power-failure-safe transaction mechanism.

### Output

For each DOY and grid cell the empirical file stores:

```text
count[doy, lat, lon, x_bin, y_bin]
n_valid[doy, lat, lon]
x_bin_left / x_bin_right
y_bin_left / y_bin_right
next_year[doy, y_chunk, x_chunk]
```

The `count` field is the primary scientific record. PDF values are derived from counts and bin area; no fitted distribution is required to reconstruct the empirical surface.



## 10A. Five-day centred fitting window

The distribution-selection layer uses a **five-day centred window** for every target calendar date:

```text
D-2  D-1  D  D+1  D+2
             ^
          target DOY
```

All available **hourly observations** in this interval are passed to the fitting layer. The window is applied at the raw time-series extraction layer, not after daily summarisation. The target day remains the output label; the neighbouring four days only enlarge the fitting sample.

For the 1981-2020 target climatology, the temporal padding is:

- `1980-12-30` and `1980-12-31` are read from the dedicated combined edge file when needed for the first target days.
- Data outside 1981-2020 are padding only; they are not counted as independent target climatology days.
- Data through 2021-06 are available, so the end of 2020 has a complete centred window.

The loader uses the **actual datetime coordinate inside NetCDF**. `valid_time` is accepted explicitly, and other datetime-like coordinates are detected rather than inferred from filenames.

### Window coverage diagnostics

Each window records at least:

```text
window_start
window_end
expected_hour_count
available_hour_count
completeness_fraction
```

This is important at the first/last boundaries and whenever one physical variable has different coverage from another.

## 10B. Distribution candidate family

The fitting engine does not assume one universal distribution. For each target day and extracted series it can compare:

```text
Normal
Skew-Normal (project-style moment fit)
Skew-Normal (MLE refinement)
Pearson Type III
Beta (bounded variables only)
Bimodal Normal / 2-component Gaussian mixture
```

The Skew-Normal candidate includes the parameterisation used in the `ClimateProcessingEngine` plugin and a separate MLE-refined candidate. The Bimodal Normal candidate follows the project's five-parameter representation:

```text
w1, mu1, sigma1, mu2, sigma2
```

with `w2 = 1 - w1`, EM fitting, deterministic multi-start, regularisation, and component ordering by mean. The model-selection record also retains multimodality diagnostics such as separation and overlap.

Candidate ranking uses:

```text
log-likelihood
AIC
AICc
BIC
```

The full candidate table is retained, so the selected model is auditable rather than silently replacing the alternatives.

## 10C. Two-dimensional probability model

The primary two-dimensional product remains **empirical-first**. A histogram/PDF based on actual paired hourly observations is the reference surface. No bivariate normal assumption is imposed.

The parametric layer is separate:

```text
marginal fit
    +
dependence fit
    -> joint model
```

For bounded variables such as `RH/100` and `q`, Beta is a candidate marginal. The dependence layer can use copula models; the current lightweight fitting layer includes a Gaussian-copula estimator on rank/pseudo-observations. Alternative copulas should be compared only after diagnostics, especially for tail dependence.

Therefore:

> **Bivariate Gaussian is a reference candidate, not the definition of the joint distribution.**

## 10D. Station/grid-cell query architecture

The five-day layer is designed for the same operational idea used in `ClimateProcessingEngine`: load only the required spatial block and cache it instead of repeatedly reading the entire grid. The current v7.5 extraction API is intended for direct station/grid-cell queries and can be embedded into the full block orchestrator.

The source project separates monthly reading, block selection, and disk/in-memory caching; v7.5 follows the same separation of concerns while preserving the ERA5-Land timestamp as the source of truth.

## 10F. Bivariate dominance report generator

The release includes `bivariate_dominance_report.py` as a **diagnostic and visual reporting tool** for the day-resolved bivariate model-selection product. It is deliberately separated from the fitting engine: the fitter decides the winning model; the report generator explains where and when each model dominates.

Expected selection product:

```text
best_model_code(doy, latitude, longitude)
```

with a global `model_names` mapping, for example:

```json
{"1":"Empirical-2D","2":"Gaussian-Copula","3":"t-Copula","4":"Beta-Copula","5":"Bimodal-Copula"}
```

The generator produces a single publication-oriented PDF containing:

- model-family legend and provenance;
- monthly dominance-share chart;
- DOY x model heatmap showing the spatial fraction selected for each climatological day;
- one spatial dominance map for each month;
- exported PNG figures for downstream papers, presentations, and QA archives.

### Scientific interpretation

The report is a **diagnostic of model selection**, not a claim that the winning family is physically true. In particular:

```text
Empirical 2-D PDF
        |
        +----> reference / non-parametric surface
        |
        +----> candidate parametric joint models
                  |
                  +--> marginal family
                  +--> dependence family / copula
                  +--> fit statistics + tail diagnostics
                           |
                           +--> best_model_code
                                    |
                                    +--> bivariate_dominance_report.py
```

No Gaussian assumption is introduced by the reporting layer. A Gaussian copula, t-copula, Beta-copula construction, or bimodal construction appears only when it is selected by the upstream model-selection contract.

### What the maps answer

For each month, the spatial map answers:

> **Which bivariate probability model most often represents the paired moisture distribution at this location during this month?**

The DOY heatmap answers:

> **On which climatological days does each bivariate model family dominate across the spatial domain?**

The monthly share chart answers:

> **How does the prevalence of each candidate family change through the annual cycle?**

### Reproducibility and QA

The report must preserve the selection-product filename, model mapping, climatological calendar, software/version metadata, and the SHA-256 of the input selection file where available. The report generator must never refit distributions or silently change model codes. It consumes the frozen selection result and visualizes it deterministically.

Example command:

```bash
python bivariate_dominance_report.py \
    moisture_bivariate_model_selection_1981_2020_v7_5.nc \
    --output-pdf bivariate_distribution_dominance_report.pdf \
    --output-dir bivariate_dominance_figures
```

This report belongs in the publication/release bundle alongside the main NetCDF, diagnostics, bivariate parameters, configuration, source checksum, and validation log.

## 10E. Why v6 is retained

v6 is retained as a historical and educational baseline. Its daily-statistical Gaussian/Monte-Carlo architecture is useful for explaining statistical propagation, covariance handling, PSD checks, and Monte-Carlo uncertainty. It is **not** the preferred production path for direct hourly moisture distributions.
