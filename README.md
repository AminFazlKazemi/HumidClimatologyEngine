# HumidClimatologyEngine

> **Research-grade multivariate moisture climatology framework for ERA5-Land**
>
> Statistical state: **(T, Td, ln P)** → Monte Carlo propagation → **RH, vapor pressure, mixing ratio, specific humidity**.

[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Data](https://img.shields.io/badge/data-ERA5--Land-1F4E79)](https://cds.climate.copernicus.eu/)
[![Output](https://img.shields.io/badge/output-NetCDF4-4B8BBE)](https://www.unidata.ucar.edu/software/netcdf/)
[![License](https://img.shields.io/badge/license-MIT-2E7D32)](LICENSE)

> **Documentation contract:** This README is the human-readable front door to the same scientific contract described in the long-form PDF. The implementation reference is `moisture_climatology_v6.py`; generated documentation must not invent functionality that is not present in the implementation.

> **Implementation reference:** `K:\kazemi\papers\temperature_interpolation\HumidClimatologyEngine\moisture_climatology_v6.py`. This absolute Windows path is retained as the project source-of-truth identifier; it is not copied into the public package as a required runtime path.
>
> **Release note:** the distributable package contains the upgraded README, the long-form scientific PDF, and the hardened reset utility. The documentation deliberately distinguishes documented implementation behavior from recommended extensions and does not invent absent modules or tests.

---

## 1. Executive summary

HumidClimatologyEngine is designed around a simple scientific principle:

> **Do not climatologize nonlinear moisture variables independently when the physically relevant atmospheric state is multivariate.**

The intended workflow models the coupled state

$$
X = (T, T_d, \ln P)
$$

for each climatological day and grid cell, retains the dependence structure through a covariance matrix, samples the joint state with deterministic Monte Carlo, and only then evaluates nonlinear moisture diagnostics.

The principal products are:

- relative humidity (RH),
- vapor pressure (e),
- mixing ratio (r),
- specific humidity (q),
- distributional moments, and
- explicit Monte Carlo / numerical diagnostics.

This architecture is particularly valuable when the research question depends on **relationships between temperature, dew point, and pressure**, rather than on marginal averages alone.

---

## 2. What is actually in this archive

```text
HumidClimatologyEngine/
├── moisture_climatology_reset.py
├── README.md
└── HumidClimatologyEngine_Professional_Guide.pdf
```

The reset utility has been hardened rather than expanded into a fictional package structure. The PDF is a regenerated, publication-oriented technical guide derived from the supplied code and documentation. It now supports:

- explicit `--base` selection,
- `--dry-run`,
- `--yes` for non-interactive execution,
- structured logging,
- path-safety checks,
- symlink-safe deletion,
- clear exit codes, and
- a single explicit allow-list of artifacts.

### What this archive does **not** contain

The original README described a larger repository containing an engine module, tests, packaging metadata, configuration files, notebooks, examples, and changelog/citation/security documents. Those files are not present in the supplied ZIP, so they are not represented here as if they shipped with the archive.

That distinction is deliberate: a publication-quality README should be **accurate before it is impressive**.

---

## 3. Scientific objective

The model state is

$$
X =
\begin{bmatrix}
T\\
T_d\\
\ln P
\end{bmatrix}.
$$

For a given climatological day and grid cell, the intended statistical representation is

$$
X \sim \mathcal{N}(\mu,\Sigma),
$$

subject to empirical diagnostics, covariance validity checks, and the limitations of the Gaussian assumption.

The transformed pressure variable is used to represent multiplicative pressure variability more naturally and to guarantee positive pressure after exponentiation.

The model is **not** a claim that the real atmosphere is Gaussian. It is a tractable local approximation whose adequacy must be validated against the intended scientific use.

---

## 4. Why joint modeling matters

Moisture diagnostics are nonlinear functions of the atmospheric state. For example,

$$
RH = 100\,\frac{e_s(T_d)}{e_s(T)},
$$

while

$$
r = \frac{\epsilon e}{P-e}, \qquad q = \frac{r}{1+r}.
$$

Averaging T, Td and P separately and then inserting those means into a nonlinear formula is generally not equivalent to propagating the observed joint distribution through the formula.

The engine therefore follows the order:

```text
joint atmospheric state
        ↓
statistical representation
        ↓
Monte Carlo realization
        ↓
thermodynamic transformation
        ↓
distributional summary
```

rather than:

```text
mean T + mean Td + mean P
        ↓
one deterministic humidity calculation
```

This distinction is central to the scientific design.

---

## 5. Thermodynamic transformations

### 5.1 Relative humidity

$$
RH = 100\,\frac{e_s(T_d)}{e_s(T)}.
$$

The precise saturation-vapor-pressure formulation is part of the model definition and therefore must be preserved with the code/version used for any published result.

### 5.2 Vapor pressure

$$
e=e_s(T_d).
$$

### 5.3 Mixing ratio

$$
r=\frac{0.622e}{P-e}.
$$

The denominator must remain positive for the intended physical domain.

### 5.4 Specific humidity

$$
q=\frac{r}{1+r}.
$$

Because these transformations are nonlinear, distributional shape is expected to change even when the underlying atmospheric state is approximately Gaussian.

---

## 6. Calendar convention

The supplied specification uses a **366-slot climatological calendar** with February 28 and February 29 pooled into a single climatological day.

| Slot | Meaning |
|---:|---|
| 1–58 | January 1 through February 27 |
| 59 | Reserved slot in the supplied specification |
| 60 | February 28 + February 29 composite |
| 61–366 | March 1 through December 31 |

Operationally:

```text
Leap year:
  Feb 28 → 60
  Feb 29 → 60
  Mar 01 → 61

Non-leap year:
  Feb 28 → 60
  Mar 01 → 61
```

This convention is **not** the standard one-dimensional day-of-year mapping and therefore must be treated as a formal model contract. Any downstream consumer must use the same mapping.

> **Important:** if a production implementation uses a different slot convention, the implementation is authoritative and this section should be updated to match it exactly before publication.

---

## 7. Statistical accumulation

The intended accumulation architecture is online rather than raw-history based.

### Paired-valid sampling

A sample contributes only when all required state variables are simultaneously valid:

```text
T   finite
Td  finite
P   finite
P   > 0
```

This is important because a covariance matrix must describe a **single joint population** rather than a set of separately filtered marginals.

### Welford-style online statistics

The design uses numerically stable sufficient-statistic accumulation for counts, means, variances, and covariances. This avoids retaining a full multidecadal sample history just to estimate first- and second-order moments.

### Mergeable states

Year-level states can be merged into climatological states. Conceptually:

```text
1981 ─┐
1982 ─┤
...   ├──→ mergeable sufficient statistics ──→ 1981–2020 climatology
2020 ─┘
```

This makes restart and parallel execution much more practical than a monolithic historical array.

### Higher moments

The intended moisture products include:

- mean,
- sample standard deviation,
- bias-corrected skewness, and
- Fisher excess kurtosis.

For publication, the exact estimator definitions and finite-sample behavior should be frozen in a versioned methodology document rather than inferred from variable names alone.

---

## 8. Covariance discipline

A covariance matrix is not merely an internal numerical object. It is part of the scientific model.

The intended safeguards are:

1. validate the number of valid paired observations;
2. inspect covariance conditioning and eigenvalues;
3. test positive-semidefiniteness within an explicit tolerance;
4. test Cholesky feasibility where Cholesky factorization is used;
5. preserve diagnostics when a cell/day is numerically problematic;
6. **never silently replace an invalid covariance with an unrelated fallback model**.

A silent identity fallback would destroy the empirical dependence structure and produce a numerically convenient but scientifically different model.

---

## 9. Monte Carlo architecture

The design explicitly avoids materializing a giant sample cube such as:

```text
N_SAMPLES × latitude × longitude
```

for an entire climatological day.

Instead, the intended processing hierarchy is:

```text
day
  → cell chunk
      → sample batch
          → transform
              → online moments
```

The dominant temporary memory footprint therefore depends primarily on:

```text
CELL_CHUNK_SIZE × SAMPLE_BATCH_SIZE
```

rather than on the complete global grid.

### Nominal configuration from the supplied specification

```text
N_SAMPLES          = 5000
CELL_CHUNK_SIZE    = 1024
SAMPLE_BATCH_SIZE  = 256
MAX_WORKERS        = 2
RANDOM_SEED        = 20260821
```

These are **configuration values, not scientific truths**. A sample-count sensitivity study should compare at least a low, middle, and high setting and verify that the scientific conclusions are stable.

A useful convergence ladder is:

```text
500 → 1000 → 2000 → 5000 → 10000
```

For each level, compare means, spread, skewness, kurtosis, and Monte Carlo standard errors.

---

## 10. Monte Carlo uncertainty

For a Monte Carlo mean,

$$
SE(\bar{x}) \approx \frac{s}{\sqrt{N}},
$$

where `s` is the realization-level standard deviation and `N` is the number of valid realizations.

This allows the final product to separate two conceptually different quantities:

```text
physical/statistical variability
                 vs.
Monte Carlo sampling error
```

Representative diagnostics described by the supplied specification are:

```text
mc_se_mean_rh
mc_se_mean_e
mc_se_mean_r
mc_se_mean_q
```

A production implementation should also preserve the **valid sample count** associated with each uncertainty estimate.

---

## 11. Checkpoint and restart philosophy

Long climate computations fail for ordinary operational reasons: process termination, workstation restarts, filesystem problems, or accidental interruption. Restartability is therefore an engineering requirement, not a luxury.

The supplied specification describes two checkpoint scales.

### Annual checkpoints

Conceptually:

```text
year_YYYY_<config_hash>.*
```

with metadata such as:

- year,
- schema version,
- grid shape,
- configuration hash,
- checksum,
- creation metadata.

### Daily Monte Carlo checkpoints

Conceptually:

```text
day_DDD_<config_hash>.*
```

This permits a long climatological run to restart from its last valid daily state.

### Atomicity

A robust checkpoint sequence is:

```text
compute
  ↓
write temporary artifact
  ↓
flush / synchronize
  ↓
atomic replacement
  ↓
write/verify metadata and checksum
```

A checkpoint that is incomplete or incompatible must be treated as invalid rather than reused.

---

## 12. Reproducibility contract

A scientifically reproducible result requires more than a software version.

Preserve at minimum:

```text
software version / commit
configuration hash
source checksum
input-data version
preprocessing recipe
random seed
sample count
calendar convention
statistical estimator definitions
environment versions
final output checksum
```

For publication, the acquisition recipe for ERA5-Land should also be preserved when the raw data cannot legally or practically be redistributed.

---

## 13. ERA5-Land input contract

The intended input variables are:

| Physical quantity | Typical ERA5-Land name | Role |
|---|---|---|
| 2 m air temperature | `t2m` | T |
| 2 m dew-point temperature | `d2m` | Td |
| Surface pressure | `sp` | P |

The scientific pipeline assumes daily fields or a precisely documented daily aggregation from a higher-frequency source.

Before production, validate:

- units,
- time coverage,
- time monotonicity,
- latitude ordering,
- longitude convention,
- grid compatibility,
- missing-value semantics,
- variable naming,
- monthly completeness.

The rule should be **fail closed**: reject ambiguous or inconsistent inputs instead of guessing.

---

## 14. Data-model requirements for production

Every input dataset should have an explicit data contract covering:

```text
dimensions
coordinates
units
valid range
missing-value convention
calendar/time zone convention
spatial reference
aggregation definition
source provenance
```

A useful preflight report should answer, before expensive computation begins:

```text
What files will be read?
What exact dates do they cover?
Are the three variables on identical grids?
Are units compatible with the physics layer?
How many observations are expected?
Are any months missing or duplicated?
```

The cost of this preflight is negligible compared with the cost of discovering an input problem after many hours of Monte Carlo execution.

---

## 15. Output contract

The intended main product is a NetCDF climatology such as:

```text
moisture_climatology_1981_2020.nc
```

Representative fields are:

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

A companion diagnostic product should contain enough information to distinguish:

```text
missing input
bad pressure
invalid covariance
invalid transformed state
insufficient samples
Monte Carlo uncertainty
```

from one another.

---

## 16. Recommended metadata

A publication-grade NetCDF should preserve metadata such as:

```text
title
summary
institution
source
references
Conventions
history
creation_time
software_name
software_version
git_commit
configuration_hash
input_data_identifier
calendar_convention
random_seed
mc_samples
```

Variables should also carry:

```text
long_name
units
standard_name (where applicable)
valid_min / valid_max (where meaningful)
```

The metadata should describe the **actual executed configuration**, not merely the default configuration in a README.

---

## 17. Validation hierarchy

The intended validation strategy should proceed from the smallest scientific unit to the complete product.

### Level 1 — scalar physics reference

Test the thermodynamic equations against an independent scalar implementation.

### Level 2 — vectorized physics

Confirm that batch/vector results match the scalar reference within an explicitly chosen numerical tolerance.

### Level 3 — statistics

Compare online estimators with a trusted offline reference on small synthetic samples.

### Level 4 — covariance

Check symmetry, eigenvalues, Cholesky behavior, and pathological small-sample cases.

### Level 5 — Monte Carlo

Demonstrate convergence with increasing `N_SAMPLES` and verify deterministic behavior for fixed seeds.

### Level 6 — checkpoint/restart

Run a small job, interrupt it, resume it, and verify that the restarted output matches the uninterrupted reference to the declared tolerance.

### Level 7 — final NetCDF

Verify dimensions, coordinates, units, finite-value behavior, physical bounds, diagnostic completeness, and checksums.

The key idea is **validation before scale**.

---

## 18. Physical plausibility checks

The transformation layer should monitor at least:

- negative or non-positive pressure,
- non-finite saturation-vapor-pressure results,
- `e >= P` for the mixing-ratio denominator,
- non-finite transformed outputs,
- suspicious RH excursions,
- unexpectedly extreme transformed tails.

A failed physical check should be classified and counted rather than hidden in a generic `NaN` bucket.

---

## 19. Performance philosophy

The correct optimization target is not “fastest possible run.” It is:

> **minimum wall-clock time subject to numerical correctness, memory stability, and reproducibility.**

The main tuning controls are:

```text
MAX_WORKERS
CELL_CHUNK_SIZE
SAMPLE_BATCH_SIZE
N_SAMPLES
```

Tune one at a time and benchmark on representative days/cells.

On Windows, process-based parallelism may increase memory pressure because workers can have substantial independent state. More workers are therefore not automatically faster.

---

## 20. Reset utility

The archive includes `moisture_climatology_reset.py`, which removes only these explicitly allowed artifacts:

```text
checkpoints_moisture_v6/
moisture_climatology_1981_2020.nc
moisture_climatology_diagnostics_1981_2020.nc
```

### Preview

```bash
python moisture_climatology_reset.py --dry-run
```

### Interactive delete

```bash
python moisture_climatology_reset.py
```

### Non-interactive delete

```bash
python moisture_climatology_reset.py --yes
```

### Alternate working directory

```bash
python moisture_climatology_reset.py --base C:\\c --yes
```

### Why this is safer than the original reset script

The original utility contained duplicated deletion logic. The upgraded version centralizes the target allow-list and adds:

- no implicit deletion outside the selected base directory,
- dry-run support,
- explicit confirmation unless `--yes` is provided,
- symlink-safe handling,
- structured exit codes,
- clearer operational messages.

The utility is intentionally **destructive only for the exact three named artifacts**. It does not recursively clean arbitrary files from the working directory.

---

## 21. Production operating procedure

A robust production run should follow this order:

```text
1. Freeze the software version
2. Freeze the configuration
3. Validate input inventories
4. Validate units and coordinates
5. Run scalar/reference physics tests
6. Run a small synthetic smoke test
7. Run a small real-data pilot
8. Inspect covariance diagnostics
9. Establish Monte Carlo convergence
10. Start full climatology
11. Preserve checkpoints
12. Validate final NetCDF
13. Archive checksums and provenance
```

Never jump directly from “the code starts” to “the climatology is publishable.”

---

## 22. Publication-grade reproducibility bundle

A strong release should contain, at minimum:

```text
software source
configuration
input-data identifier
preprocessing recipe
environment lock / package versions
validation report
main NetCDF
diagnostic NetCDF
checksums
methodology document
release tag / commit
```

If the raw ERA5-Land files cannot be redistributed, publish their acquisition instructions and immutable identifiers/checksums instead.

---

## 23. Scientific limitations

### Gaussian state model

A multivariate normal model may not represent skewness, multimodality, tails, or regime mixtures adequately. This is especially important for extremes.

### Reanalysis uncertainty

ERA5-Land is a model/reanalysis product; the resulting climatology inherits source-data limitations and does not represent observational truth without qualification.

### Monte Carlo error

Finite `N_SAMPLES` produces sampling uncertainty even when the upstream covariance estimate is exact.

### Thermodynamic formulation

Different saturation-vapor-pressure formulations and unit conversions can produce measurable differences. The exact implementation must therefore be versioned.

### Calendar compression

Pooling February 28 and February 29 is a deliberate modeling choice. It changes the effective climatological calendar and should be documented in publications.

### Independence assumptions

A per-day Monte Carlo state model summarizes the distribution conditional on each climatological day. It does not automatically preserve temporal persistence or event sequencing across days.

---

## 24. Recommended sensitivity matrix

For a serious scientific analysis, vary at least:

| Factor | Example levels |
|---|---|
| Monte Carlo samples | 500, 1000, 2000, 5000, 10000 |
| workers | 1, 2, 4 |
| cell chunk | 512, 1024, 2048 |
| sample batch | 128, 256, 512 |
| covariance tolerance | documented low/central/high values |
| representative regions | dry / humid / cold / warm / mountainous |

The performance factors should be judged on computational behavior; the scientific factors should be judged on output stability.

---

## 25. Reproducibility acceptance criteria

A run should be considered reproducible only when:

```text
same inputs
+ same configuration
+ same code
+ same random seed
        ↓
consistent output
```

within the predeclared numerical tolerance.

For parallel calculations, deterministic results may additionally depend on reduction order and random-stream design. Therefore, “same seed” is not by itself a guarantee of bit-for-bit identity unless the implementation explicitly guarantees it.

This distinction should be stated in any reproducibility claim.

---

## 26. Quality gates before publication

Use a signed or versioned release gate such as:

- [ ] input inventory complete;
- [ ] no duplicate months;
- [ ] units verified;
- [ ] coordinates verified;
- [ ] leap-day mapping confirmed;
- [ ] scalar physics reference passes;
- [ ] online statistics agree with an offline reference;
- [ ] covariance diagnostics acceptable;
- [ ] Monte Carlo convergence demonstrated;
- [ ] restart test passes;
- [ ] NetCDF structure verified;
- [ ] physical plausibility checks passed;
- [ ] configuration archived;
- [ ] source checksum archived;
- [ ] input-data provenance archived;
- [ ] final output checksum archived.

A checklist is not a substitute for tests, but it is a powerful final defense against procedural omissions.

---

## 27. Minimal verification of the supplied archive

### Syntax check

```bash
python -m py_compile moisture_climatology_reset.py
```

### Preview reset behavior

```bash
python moisture_climatology_reset.py --dry-run
```

### Inspect the exact target allow-list

```bash
python -c "import moisture_climatology_reset as m; print(m.TARGETS)"
```

These checks verify the **artifact that is actually present** in this archive. They do not claim that an absent climatology engine module has been executed.

---

## 28. Recommended repository structure for the next full release

The scientific README in the original archive points toward a larger architecture. A robust future repository can formalize it as:

```text
HumidClimatologyEngine/
├── src/humidclimatology/
│   ├── calendar.py
│   ├── physics.py
│   ├── statistics.py
│   ├── covariance.py
│   ├── monte_carlo.py
│   ├── checkpoints.py
│   ├── io.py
│   ├── diagnostics.py
│   ├── validation.py
│   └── cli.py
├── tests/
├── configs/
├── docs/
├── examples/
├── scripts/
├── pyproject.toml
├── CITATION.cff
├── CHANGELOG.md
├── CONTRIBUTING.md
├── SECURITY.md
└── README.md
```

This is a **recommended target architecture**, not a claim that these files exist in the supplied ZIP.

---

## 29. Engineering principles

The project should preserve these rules as non-negotiable invariants:

> **Fail closed on ambiguous inputs.**

> **Never hide a scientific numerical failure.**

> **Never destroy dependence structure with an unrelated fallback.**

> **Never optimize by changing the model silently.**

> **Never call a result reproducible without preserving the actual configuration.**

> **Never treat a successful process exit as proof of scientific validity.**

These principles are more important than any particular package or implementation detail.

---

## 30. Maintainer note

**Maintainer:** Amin Fazlkazemi

Before public release, replace any placeholder repository URLs and ensure the software version, DOI/release identifier, and citation metadata match the exact code used for the published result.

---

## License

MIT. See `LICENSE` when the full repository is assembled.

---

# 18. Deep Scientific Reference

This section expands the operational README into a reference specification. It is intentionally explicit about what is a model definition, what is a numerical implementation choice, and what is a recommended validation activity.

## 18.1 Scientific contract

The project is best understood as a sequence of deterministic transformations applied to a long atmospheric archive:

```text
ERA5-Land state variables
        |
        +-- acquisition + provenance
        |
        +-- unit normalization
        |
        +-- climatological-day mapping
        |
        +-- paired-valid filtering
        |
        +-- online joint-state statistics
        |
        +-- annual checkpoint states
        |
        +-- mergeable multi-year state
        |
        +-- covariance / correlation validation
        |
        +-- multivariate Monte Carlo propagation
        |
        +-- thermodynamic transformations
        |
        +-- higher-moment accumulation
        |
        +-- diagnostics + NetCDF serialization
```

A scientific run should therefore be reproducible not only from the final NetCDF files but from the exact chain of inputs, calendar rules, numerical tolerances, random seed, configuration, software content, and checkpoint lineage that produced them.

## 18.2 Source-of-truth rule

The scientific implementation reference for the project is `moisture_climatology_v6.py`, identified in the project documentation as the main v6 engine. Any wrapper, reset utility, example, README section, or user-facing command should remain subordinate to that implementation. A documentation change must never silently redefine the numerical algorithm.

For release engineering, the safest rule is:

> Change the implementation first, validate it, then regenerate documentation from the validated implementation contract.

This prevents documentation drift, especially when a parameter name, checkpoint schema, default value, or calendar rule changes.

## 18.3 State vector and units

The core state is

$$
X = (T, T_d, \ln P)^T.
$$

The working units documented by the project are degrees Celsius for temperature, hectopascals for pressure after conversion, and the natural logarithm of pressure for the third state variable. The expected ERA5-Land source units are Kelvin for temperature and dew point and Pascal for surface pressure.

The transformation is therefore

$$
T_C = T_K - 273.15,
$$
$$
T_{d,C} = T_{d,K} - 273.15,
$$
$$
P_{hPa} = P_{Pa}/100,
$$
$$
L = \ln(P_{hPa}).
$$

The explicit pressure-domain transformation is important because the sampler must not generate physically impossible negative surface pressure values.

## 18.4 Paired-valid data rule

An observation enters the joint-state accumulator only when all three state variables are finite and pressure is strictly positive. This rule is stronger than three independent marginal masks. The latter would create mismatched sample sets and corrupt covariance estimates.

For a valid row,

$$
X_i = \begin{bmatrix}T_i \\ T_{d,i} \\ \ln P_i\end{bmatrix}.
$$

For $n$ paired-valid states,

$$
\mu = \frac{1}{n}\sum_{i=1}^{n}X_i,
$$

and the sample covariance is

$$
\Sigma = \frac{1}{n-1}\sum_{i=1}^{n}(X_i-\mu)(X_i-\mu)^T.
$$

The paired-valid count is itself a scientific diagnostic and should be mapped or summarized before interpreting downstream moisture statistics.

## 18.5 Why Welford-style accumulation is the correct engineering choice

The archive can span decades and a large spatial grid. A naive implementation that loads the full record and repeatedly computes sums and sums of squares is both memory-hungry and more vulnerable to cancellation. Online accumulation keeps only sufficient state for the configured statistic order.

For the scalar mean/M2 form:

$$
\delta = x-\mu,
$$
$$
\mu' = \mu + \delta/n',
$$
$$
M_2' = M_2 + \delta(x-\mu').
$$

The sample variance is $M_2/(n-1)$ for $n\ge2$.

The project extends the same principle to paired cross-products. This lets annual states be merged without replaying raw observations.

## 18.6 Mergeability

For two states with counts $n_1,n_2$, means $m_1,m_2$, and second-order accumulators $M_{2,1},M_{2,2}$, define

$$
\delta=m_2-m_1,
$$
$$
n=n_1+n_2,
$$
$$
m=m_1+\delta\frac{n_2}{n}.
$$

Then

$$
M_2=M_{2,1}+M_{2,2}+\delta^2\frac{n_1n_2}{n}.
$$

The covariance cross-term has the compatible correction

$$
C_{xy}=C_{xy,1}+C_{xy,2}+(m_{x,2}-m_{x,1})(m_{y,2}-m_{y,1})\frac{n_1n_2}{n}.
$$

This is why annual checkpoint states are scientifically meaningful: they are mergeable sufficient states for the modeled second-order dependence structure.

## 18.7 Climatological calendar contract

The project uses a 366-slot coordinate system with an explicit leap-day convention:

| Calendar date | Gregorian DOY | Climatological slot |
|---|---:|---:|
| Feb 28 in leap year | 59 | 60 |
| Feb 29 in leap year | 60 | 60 |
| Feb 28 in non-leap year | 59 | 60 |
| Mar 1 in leap year | 61 | 61 |
| Mar 1 in non-leap year | 60 | 61 |
| Mar 2 | 62/61 | 62 |

The practical consequence is that **February 28 and February 29 are deliberately pooled into the same climatological slot**. Slot 59 remains reserved under the project's convention. This is a model decision, not a trivial indexing shortcut.

Any downstream comparison with a conventional 365-day climatology must therefore document its own mapping back to a 365-slot representation.

## 18.8 Thermodynamic layer

The documented water-phase saturation expression is

$$
e_{s,w}(T)=6.112\exp\left(\frac{17.67T}{T+243.5}\right),
$$

while the documented ice-phase expression is

$$
e_{s,i}(T)=6.112\exp\left(\frac{22.46T}{T+272.62}\right).
$$

The implementation uses water for $T\ge0^\circ C$ and ice for $T<0^\circ C$.

Actual vapor pressure is obtained from dew point:

$$
e=e_s(T_d).
$$

Relative humidity is formed as

$$
RH_{raw}=100\frac{e_s(T_d)}{e_s(T)},
$$

then bounded to the physical reporting range while the supersaturation fraction is retained as a diagnostic.

The mixing ratio is

$$
r=\frac{\epsilon e}{P-e},\quad \epsilon=0.622,
$$

subject to $e>0$, $P>0$, and $e<P$. Specific humidity follows from

$$
q=\frac{r}{1+r}.
$$

The distinction between physical validity filtering and display clipping is important: clipping RH to 100% does not make an invalid pressure partition valid.

## 18.9 Joint Gaussian approximation

For each climatological day and grid cell the project represents

$$
X\sim\mathcal N(\mu,\Sigma).
$$

The benefit is preservation of dependence between temperature, dew point and log-pressure during stochastic propagation. The cost is a structural assumption that may be weak in multimodal or regime-switching climates.

Potentially vulnerable situations include:

- strongly bimodal synoptic regimes,
- coastal and mountainous cells with mixed air masses,
- very cold regimes near phase-transition boundaries,
- small observational sample sizes,
- tail-sensitive applications.

For these cases, a Gaussian copula, mixture model, or regime-conditioned model can be evaluated as a sensitivity extension. Such extensions should be treated as alternative methods, not invisible changes to the baseline engine.

## 18.10 Covariance validity and repair policy

A correlation or covariance matrix must be positive semidefinite. Small negative eigenvalues can arise from finite precision, especially after correlation clipping and reconstruction.

The documented policy is to inspect the minimum eigenvalue, allow only a controlled repair window, perform a nearest-correlation style repair where appropriate, and reject a cell/day when the matrix remains invalid. This is scientifically safer than an unconditional jitter that silently changes every matrix.

The three most important audit products are:

1. minimum eigenvalue,
2. fraction of cells requiring repair or rejection,
3. spatial distribution of covariance failures.

A modeler should review these maps before accepting the final climatology.

## 18.11 Monte Carlo propagation

Let

$$
Z\sim\mathcal N(0,I_3),
$$

and let $L$ be a valid Cholesky factor of the covariance model. Then

$$
X=\mu+LZ.
$$

The pressure state is transformed back through

$$
P=\exp(X_3).
$$

The engine processes samples in batches and spatial cells in chunks. This avoids allocating the full tensor with dimensions approximately $N_{samples}\times N_{lat}\times N_{lon}$.

The documented production defaults include 5,000 Monte Carlo samples, a cell chunk of 1,024, a sample batch of 256, a conservative two-worker default, and a fixed random seed. These are configuration defaults, not universal scientific constants.

## 18.12 Monte Carlo convergence

For a simple mean under independent sampling,

$$
SE_{MC}\propto N^{-1/2}.
$$

The practical implication is that quadrupling the sample count approximately halves idealized mean sampling error, while also increasing computational work by about a factor of four.

A serious release should compare a sequence such as 500, 1,000, 2,000, 5,000, and 10,000 samples. Means usually converge faster than kurtosis; therefore convergence must be evaluated statistic by statistic.

## 18.13 Higher moments

The engine retains more than location and spread. For transformed moisture variables it reports mean, sample standard deviation, bias-corrected skewness, and Fisher excess kurtosis.

The bias-corrected skewness follows the SciPy-style correction described in the project documentation. Fisher excess kurtosis is zero for a Gaussian population, positive for heavier tails, and negative for lighter tails.

Higher moments are much more sample-sensitive than means. A visually stable mean map is not evidence that the fourth moment is stable.

## 18.14 Output contract

The main NetCDF product uses dimensions conceptually equivalent to

```text
doy × latitude × longitude
```

and includes the moisture fields:

| Variable family | Statistics |
|---|---|
| relative humidity | mean, std, skew, kurtosis |
| vapor pressure | mean, std, skew, kurtosis |
| mixing ratio | mean, std, skew, kurtosis |
| specific humidity | mean, std, skew, kurtosis |

The diagnostic product includes supersaturation fraction, invalid $e/P$ fraction, covariance validity indicators, minimum eigenvalue, valid Monte Carlo sample count, three pairwise correlations, and paired-valid historical observation count.

The documented fill-value convention uses `-9999` where appropriate. Downstream software should honor `_FillValue` metadata rather than interpret that sentinel as a physical value.

## 18.15 Checkpoint and restart contract

The architecture separates three persistence layers:

```text
annual state
   ↓ merge
multi-year state
   ↓ daily propagation
transformed moment state
   ↓ finalize
NetCDF product
```

A restart is scientifically valid only when the checkpoint is compatible with the same schema and configuration contract. Required provenance includes the exact source content, configuration, input manifest and hashes, environment specification, checkpoint hashes, final output hashes, validation logs, and Monte Carlo convergence results.

## 18.16 Reproducibility record

A publication-grade run should be reconstructible from a run manifest containing at least:

```yaml
project: HumidClimatologyEngine
implementation: moisture_climatology_v6.py
schema_version: 6.0
period: 1981-2020
calendar_slots: 366
leap_day_policy: pool_feb28_feb29_into_slot_60
state: [T, Td, lnP]
random_seed: 20260821
mc_samples: 5000
cell_chunk_size: 1024
sample_batch_size: 256
max_workers: 2
input_manifest: SHA256...
script_sha256: SHA256...
configuration_hash: SHA256...
created_utc: 2026-08-22T...
```

The precise values must be taken from the executed configuration, not copied from an example.

## 18.17 QA release gate

A release candidate should pass all of the following conceptual gates:

| Gate | Acceptance question |
|---|---|
| calendar | Do known leap/non-leap dates map exactly to the documented slots? |
| statistics | Do online and batch statistics agree to numerical tolerance? |
| covariance | Are symmetry, eigenvalues, and factorization diagnostics valid? |
| physics | Do independent formula checks reproduce expected values? |
| Monte Carlo | Does increasing N produce convergent summaries? |
| restart | Does checkpoint + restart reproduce an uninterrupted run within tolerance? |
| serialization | Are dimensions, units, fill values, and metadata consistent? |
| provenance | Are configuration, source, and input hashes archived? |
| operations | Does the workflow fail loudly on missing or inconsistent inputs? |

## 18.18 Memory engineering

Peak memory is driven primarily by the number of cells in a chunk, the Monte Carlo batch size, and the number of simultaneously resident arrays. The correct tuning direction on a constrained workstation is to reduce chunk or batch size before changing the scientific model.

When RAM pressure is high:

1. decrease `CELL_CHUNK_SIZE`,
2. decrease `SAMPLE_BATCH_SIZE`,
3. reduce process concurrency,
4. avoid hidden copies caused by dtype conversion,
5. process fewer years per worker only if the checkpoint design supports it.

When CPU utilization is low, inspect I/O waits and serialization before simply increasing worker count. The documented two-worker default is intentionally conservative.

## 18.19 Failure taxonomy

Most failures should fall into one of a small number of recognizable classes:

**Input/provenance failure** - missing files, wrong units, date gaps, duplicate files, or mixed time conventions.

**Calendar failure** - a post-February day shift or inconsistent leap-year handling.

**Statistical failure** - insufficient paired observations, non-finite moments, or inconsistent merge states.

**Covariance failure** - materially non-PSD correlation/covariance structure.

**Physical failure** - invalid pressure partition, non-finite saturation pressure, or excessive supersaturation.

**Monte Carlo failure** - insufficient valid draws or non-converged higher moments.

**Persistence failure** - checksum mismatch, incompatible checkpoint schema, or incomplete NetCDF output.

A good production log should make the failure class obvious without reading the entire trace.

## 18.20 Scientific interpretation of the four products

### Relative humidity

RH is dimensionless after conversion to percent and describes vapor pressure relative to saturation at temperature. It is intuitive, but because it depends on both $T$ and $T_d$, its distribution is strongly nonlinear.

### Vapor pressure

Vapor pressure is directly related to the actual water-vapor partial pressure and is often useful for thermodynamic diagnostics where relative humidity's temperature dependence is undesirable.

### Mixing ratio

Mixing ratio measures water-vapor mass relative to dry-air mass under the pressure-partitioning formulation used by the project. The nonlinear denominator makes it sensitive to invalid or near-invalid $e/P$ states.

### Specific humidity

Specific humidity is water-vapor mass per total moist-air mass under the adopted formulation. It is bounded and often convenient for transport and moisture-budget applications.

The four products should not be interpreted as interchangeable. They answer different scientific questions and inherit different nonlinearities from the same joint atmospheric state.

## 18.21 Uncertainty taxonomy

Monte Carlo propagation quantifies uncertainty arising from the adopted stochastic representation of the state and its nonlinear transformation. It does **not**, by itself, quantify all scientific uncertainty.

Separate uncertainty families include:

- parameter estimation uncertainty,
- Gaussian-model structural uncertainty,
- reanalysis uncertainty,
- temporal dependence effects,
- spatial representativeness,
- calendar-convention sensitivity,
- thermodynamic formula sensitivity,
- Monte Carlo sampling error,
- data-acquisition and preprocessing uncertainty.

The distinction should appear explicitly in publications and reports.

## 18.22 Sensitivity design

A robust sensitivity campaign can vary one scientific decision at a time while preserving a frozen baseline. High-value dimensions include Monte Carlo sample count, PSD repair tolerance, Gaussian versus alternative dependence model, calendar convention, minimum observation threshold, and vapor-pressure formulation.

For each sensitivity axis, record:

```text
baseline configuration
alternative configuration
maximum absolute difference
median absolute difference
relative difference where defined
spatial pattern of change
scientific conclusion changed? yes/no
```

The purpose is not to make every choice uncertainty-free. It is to show whether the scientific conclusion is robust to reasonable modeling decisions.

## 18.23 Publication-ready methods statement

A compact methods statement derived from the documented design is:

> Daily ERA5-Land 2-m air temperature, 2-m dew-point temperature, and surface pressure were normalized to a 366-slot climatological calendar in which February 28 and February 29 were pooled. For each climatological day and grid cell, paired-valid observations were summarized in the joint state $(T,T_d,\ln P)$ using numerically stable online and mergeable moment accumulation. The resulting dependence model was propagated through phase-aware saturation vapor-pressure equations using streaming multivariate Monte Carlo, and relative humidity, vapor pressure, mixing ratio, and specific humidity were summarized by mean, sample standard deviation, bias-corrected skewness, and Fisher excess kurtosis. Checkpoints, configuration hashes, source hashes, diagnostics, and NetCDF metadata were retained to support restartability and reproducibility.

This wording should be adapted to the exact executed configuration before publication.

## 18.24 Reviewer checklist

Before accepting a result, a reviewer should ask:

- Is the exact version of `moisture_climatology_v6.py` identified?
- Is the ERA5-Land retrieval request archived?
- Are source units and internal units explicit?
- Is the leap-day convention reported?
- Is paired-valid screening documented?
- Is the Gaussian assumption justified or sensitivity-tested?
- Are covariance repairs quantified?
- Is Monte Carlo convergence demonstrated for higher moments?
- Are invalid $e/P$ realizations reported?
- Are restart and uninterrupted runs compared?
- Are the final NetCDF files checksummed?
- Can another researcher reconstruct the environment and configuration?

## 18.25 Maintainer rules

1. Treat `moisture_climatology_v6.py` as scientific source of truth.
2. Never change a default without recording it in the release notes.
3. Never change the leap-day convention silently.
4. Never hide covariance repairs.
5. Never treat Monte Carlo convergence as proof of total physical uncertainty.
6. Never commit CDS credentials.
7. Never report undocumented modules as shipped components.
8. Regenerate README and PDF when the scientific contract changes.
9. Run the complete QA gate before tagging a release.
10. Preserve hashes and provenance with final NetCDF artifacts.

# 19. Extended Configuration Reference

| Parameter | Meaning | Scientific consequence | Operational consequence |
|---|---|---|---|
| `PERIOD` | Historical climatology interval | Defines sample population | Controls input volume |
| `N_SAMPLES` | MC draws per day/cell | Controls propagation sampling error | Main CPU driver |
| `CELL_CHUNK_SIZE` | Spatial chunk | No direct scientific change | Peak RAM |
| `SAMPLE_BATCH_SIZE` | MC batch size | No direct scientific change | Peak RAM + vectorization |
| `MAX_WORKERS` | Process concurrency | No direct scientific change | CPU/RAM/I/O contention |
| `RANDOM_SEED` | Deterministic random stream | Fixes reproducibility of sampling | Enables exact/near-exact reruns |
| `MIN_OBS` | Minimum historical paired observations | Controls where a daily/cell model is considered estimable | Missing-data coverage |
| `PSD_REPAIR_TOL` | Repair threshold | Controls acceptable numerical correction | Failure/rejection frequency |
| `SCHEMA_VERSION` | Checkpoint/output contract | Defines compatibility | Restart eligibility |

Treat values shown in the current project documentation as defaults. The executed run configuration is authoritative.

# 20. Data Lineage Template

A practical lineage record should contain:

```text
Dataset name:
Product/version:
Source URL or catalogue identifier:
Variables:
Years:
Months:
Daily statistic:
Sub-daily frequency:
UTC/local-day convention:
Spatial subset:
Retrieval timestamp:
Input files:
Input SHA-256:
Preprocessing script:
Preprocessing SHA-256:
Engine script:
Engine SHA-256:
Configuration hash:
Checkpoint set:
Final NetCDF SHA-256:
Validation log:
Analyst/release identifier:
```

# 21. Reproducible Experiment Template

```yaml
experiment:
  id: HCE-YYYYMMDD-HHMM
  objective: "daily moisture climatology"
  period: "1981-2020"
  dataset: "ERA5-Land"
  calendar:
    slots: 366
    feb28_feb29_pool: true
    reserved_slot: 59
  state:
    variables: [T_C, Td_C, ln_P_hPa]
  propagation:
    method: multivariate_normal_streaming_mc
    samples: 5000
    seed: 20260821
  memory:
    cell_chunk: 1024
    sample_batch: 256
    workers: 2
  outputs:
    format: NetCDF4
    diagnostics: true
  integrity:
    script_sha256: "..."
    config_hash: "..."
```

# 22. Minimal Operational Playbook

### Before starting

Confirm the input archive is complete, units are correct, date/time conventions are frozen, and enough paired observations exist for the intended climatological period.

### During accumulation

Watch paired-valid counts, annual state sizes, checksum generation, and memory utilization. A successful process that silently drops large portions of the record is not a successful scientific run.

### During Monte Carlo

Watch valid sample counts, supersaturation, invalid $e/P$, covariance validity, and convergence indicators. Do not accept a run only because it completed without an exception.

### Before publication

Archive final NetCDF files, diagnostics, configuration, source checksum, input manifest, environment specification, validation summary, and the exact README/PDF corresponding to the release.

# 23. Versioning and Documentation Synchronization

A recommended release unit is:

```text
source code
+ configuration schema
+ README
+ scientific PDF
+ tests
+ provenance manifest
+ changelog
+ output schema
```

The README is the first-line human interface. The PDF is the long-form scientific reference. Both should describe the same versioned contract. When a formula, variable, default, schema field, or calendar rule changes, both documents should be regenerated in the same release.

# 24. Final Scientific Assessment

The architecture is strong for a reproducible, memory-aware climatological workflow because it separates the empirical state-estimation problem from the nonlinear moisture transformation problem. Online/mergeable statistics, explicit joint dependence, controlled covariance validation, streaming Monte Carlo, and checkpoint-aware persistence make the design scalable and auditable.

The main scientific caution is not numerical; it is model adequacy. A Gaussian state model can be a useful approximation, but its validity depends on season, terrain, synoptic regime, sample size, and the statistic being reported. Higher moments are especially fragile. The correct interpretation of the system is therefore:

> a transparent stochastic propagation framework whose credibility comes from validation and sensitivity analysis, not from the Gaussian assumption alone.

# Appendix A - Formula Sheet

| Quantity | Formula |
|---|---|
| Celsius temperature | $T_C=T_K-273.15$ |
| Celsius dew point | $T_{d,C}=T_{d,K}-273.15$ |
| Pressure hPa | $P_{hPa}=P_{Pa}/100$ |
| Log pressure | $L=\ln(P_{hPa})$ |
| Water saturation pressure | $e_{s,w}(T)=6.112\exp(17.67T/(T+243.5))$ |
| Ice saturation pressure | $e_{s,i}(T)=6.112\exp(22.46T/(T+272.62))$ |
| Vapor pressure | $e=e_s(T_d)$ |
| Relative humidity | $RH=100e_s(T_d)/e_s(T)$, bounded for reporting |
| Mixing ratio | $r=0.622e/(P-e)$ |
| Specific humidity | $q=r/(1+r)$ |
| Sample variance | $s^2=M_2/(n-1)$ |
| MC standard error scaling | $SE\propto N^{-1/2}$ |

# Appendix B - NetCDF Naming Contract

The documented principal naming families are:

```text
mean_rh
std_rh
skew_rh
kurt_rh
mean_vapor_pressure
std_vapor_pressure
skew_vapor_pressure
kurt_vapor_pressure
mean_mixing_ratio
std_mixing_ratio
skew_mixing_ratio
kurt_mixing_ratio
mean_specific_humidity
std_specific_humidity
skew_specific_humidity
kurt_specific_humidity
```

Diagnostics:

```text
supersaturation_fraction
invalid_e_over_p_fraction
invalid_covariance_fraction
min_eigenvalue
valid_sample_count
corr_T_Td
corr_T_logP
corr_Td_logP
valid_observation_count
```

Coordinates/metadata:

```text
doy
latitude
longitude
month
day
reserved_day
```

# Appendix C - Quality Gate Checklist

- [ ] Source checksum recorded.
- [ ] Input manifest recorded.
- [ ] Configuration hash recorded.
- [ ] Calendar tests pass.
- [ ] Paired-valid logic tested.
- [ ] Online versus batch statistics agree.
- [ ] Covariance symmetry validated.
- [ ] Eigenvalue diagnostics reviewed.
- [ ] Physics ground-truth tests pass.
- [ ] MC convergence reviewed.
- [ ] Higher-moment stability reviewed.
- [ ] Restart equivalence checked.
- [ ] NetCDF metadata audited.
- [ ] SHA-256 of final outputs recorded.
- [ ] README and PDF version synchronized.

# Appendix D - README/PDF Synchronization Rule

The README and PDF in a release should always carry the same:

- software version,
- schema version,
- climatology period,
- calendar convention,
- default Monte Carlo settings,
- physical formula set,
- output naming contract,
- limitations statement,
- and provenance expectations.

When one changes, rebuild both.
