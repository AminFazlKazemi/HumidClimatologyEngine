# HumidClimatologyEngine v11.5

> **Public release:** HumidClimatologyEngine v11.5  
> **Primary role:** reproducible empirical moisture climatology from hourly ERA5-Land data  
> **Scientific baseline:** v11.5  
> **Release status:** public release documentation

---

## 1. What this project is

HumidClimatologyEngine v11.5 is a research-oriented climate-data processing engine for constructing long-term empirical climatologies of atmospheric moisture-related variables from hourly ERA5-Land inputs.

The engine is designed around a simple scientific principle: **retain the information carried by the hourly observations for as long as practical, derive moisture variables at the observation level, and accumulate auditable statistical state without requiring the complete multidecadal raw archive to reside in memory.**

The primary input families are:

- 2-m air temperature (`T2m` / `t2m`)
- 2-m dew-point temperature (`D2m` / `d2m`)
- surface pressure (`SP` / `sp`)

From these inputs, the engine derives and accumulates moisture variables including:

- Relative humidity (`RH` / `rh`)
- Vapor pressure (`e`)
- Mixing ratio (`r`)
- Specific humidity (`q`)

The system is not a forecast model. It does not attempt to infer causal climate mechanisms merely from statistical association. It is an empirical climatology engine whose main output is a reproducible statistical description of the supplied reanalysis product.

---

## 2. Why v11.5 exists

The project evolved from a compact hourly climatology workflow into a much more explicit statistical and operational system. The central objective of v11.5 is not to make a visually larger program; it is to preserve scientific information, make statistical assumptions explicit, make long runs restartable, and make outputs auditable.

The architecture provides:

1. Explicit input/grid validation.
2. Timestamp-based handling of the monthly ERA5-Land input triplets.
3. Spatial alignment of the three input variables.
4. Observation-level psychrometric transformation.
5. Multiple temporal resolutions in one statistical state model.
6. Mergeable online central moments through fourth order.
7. Exact minimum and maximum tracking.
8. Scalar threshold counters.
9. Pair-specific dependence statistics.
10. An empirical RH-q joint histogram layer.
11. Decade-level products and a full 1981-2020 accumulation.
12. Transactional checkpoint writing.
13. Before-image/journal based recovery semantics.
14. Explicit audit and merge-audit operations.
15. Structured progress reporting during long calculations.
16. A bounded real-data pilot path for safe pre-production testing.

---

## 3. Scientific scope and boundaries

### 3.1 What the engine can answer

The accumulated products can support questions such as:

- What is the climatological mean RH, `e`, `r`, or `q` for a given climatological day and grid cell?
- How variable is moisture at that time of year?
- What are the skewness and kurtosis diagnostics where sample support permits them?
- What are the minimum and maximum observed/valid values represented in the configured state?
- How often do selected threshold exceedances occur?
- How does the moisture regime differ by time-of-day class?
- How does a distribution or threshold frequency differ between decades?
- How strong is the empirical or moment-based dependence between configured variable pairs?
- How does the joint RH-q population occupy the configured histogram support?

### 3.2 What the engine does not claim

The engine should not be interpreted as:

- a weather forecast system;
- a causal attribution framework;
- a substitute for station observations;
- a universal distribution selector;
- a guarantee that every statistical threshold corresponds to an impact event;
- a guarantee that the reanalysis values represent direct observational truth at every location.

ERA5-Land is a reanalysis/model product. The resulting climatology is therefore a climatology **of the source product**.

---

## 4. Relationship to historical versions

v8 is preserved as an historical and methodological baseline. v11.5 is the public release documented by this repository.

The important evolution is architectural and informational. The project moved from a compact direct-processing workflow to an explicitly documented engine that retains more temporal structure, statistical state, dependence information, checkpoint state, provenance and audit evidence.

Historical versions should not be overwritten merely to make the repository look cleaner. Keeping the lineage visible is part of reproducibility.

The public release identity for current project communication is:

**HumidClimatologyEngine v11.5**

---

## 5. Primary input contract

The standard production input contract consists of one aligned monthly file triplet per year and month:

```text
T2m + D2m + Surface Pressure
```

The default locations configured in the public source are:

```text
F:\Kazemi\era5\land\T2m
F:\Kazemi\era5\land\Dew_Point_Temperature
F:\Kazemi\era5\land\Surface_Pressure
```

The default output root is:

```text
C:\c\HumidClimatologyEngine_v11.5
```

These are defaults, not a claim that every installation must use the same Windows drive letters.

### 5.1 File identity

Filenames are useful inventory aids, but the authoritative interpretation is derived from the dataset coordinates and metadata validated by the engine.

For each requested year/month, the engine expects a unique and compatible T2m/D2m/SP triplet. Missing or duplicate monthly coverage is treated as an input integrity failure rather than silently substituted.

### 5.2 Units

The input validation layer examines variables and their declared units before processing. The scientific conversion layer is implemented to normalize supported physical representations before moisture derivation.

A production run should never be started solely because files exist. The inventory and unit checks are part of the run contract.

---

## 6. Grid alignment

T2m is treated as the spatial reference for the aligned triplet. When D2m or SP uses an equivalent grid but a reversed latitude axis, the engine normalizes the orientation before scientific processing.

The log message:

```text
GRID ALIGN | D2m | reversed latitude axis to T2m reference
```

therefore represents an explicit coordinate normalization step rather than an error.

The important scientific rule is that alignment is explicit and auditable. The engine should not rely on accidental array ordering.

---

## 7. Observation-level moisture derivation

At the hourly observation level, the engine derives a moisture state from temperature, dew point and pressure rather than first collapsing the raw physical drivers into an early daily statistic.

The derived quantities are represented by:

```text
rh = relative humidity
 e = vapor pressure
 r = mixing ratio
 q = specific humidity
```

The derivation layer also records validity conditions relevant to the thermodynamic transformation, including physically invalid denominator situations and supersaturation diagnostics.

This structure is important because nonlinear moisture diagnostics can behave differently from statistics of their source variables. The engine therefore keeps the nonlinear transformation close to the observation level.

---

## 8. Temporal model

One of the defining characteristics of v11.5 is the three-level temporal representation.

### L1

L1 is the pooled daily climatological state.

It contains one bin:

```text
L1 = 1 bin
```

### L2

L2 separates the day into eight three-hour periods:

```text
L2 = 8 bins
```

### L3

L3 retains all 24 hourly positions:

```text
L3 = 24 bins
```

### Combined state

The persisted temporal dimension therefore contains:

```text
1 + 8 + 24 = 33 bins
```

This distinction is fundamental to downstream analysis. The 33-bin state is the temporal statistical state and is **not** the same thing as the histogram binning used by the RH-q empirical joint product.

---

## 9. Scalar online statistical state

For the configured scalar variables, the engine maintains mergeable statistical state rather than retaining all raw historical samples in RAM.

The conceptual state includes:

```text
n
mean
M2
M3
M4
min
max
missing counts
threshold counts
```

This provides several advantages:

- controlled memory usage;
- exact accumulation of sample count;
- deterministic mergeability;
- higher-order shape diagnostics;
- explicit support accounting;
- compatibility with decade-to-FULL aggregation.

The statistical state is finalized only where the available sample support is sufficient for the requested diagnostic.

---

## 10. Central-moment mathematics

The engine uses mergeable central-moment accumulation. The stored moments are not merely descriptive metadata; they are part of the statistical state contract.

A principal reason to use mergeable moments is that the climatology is partitioned operationally while remaining mathematically combinable. A decade state can therefore be combined with another compatible state without replaying the raw archive.

Variance is derived from the accumulated second central moment under the configured sample convention. Higher-order quantities are only finalized where the relevant sample-size and numerical validity conditions are satisfied.

The implementation should not be modified to calculate a different statistic while retaining the same schema label. A change in statistical meaning is a release-level change.

---

## 11. Extremes and thresholds

Minimum and maximum are first-class outputs.

The engine also maintains configured scalar threshold counters for RH, vapor pressure, mixing ratio and specific humidity.

The threshold counters are intentionally stored as exact integer counts rather than only as percentages. This keeps the denominator and valid-sample definition visible for audit and downstream recomputation.

Joint threshold counters are also supported for selected variable combinations.

---

## 12. Pair-specific dependence

The pair state is calculated from pair-valid samples, not from assumptions about identical marginal masks.

The configured pair set in the v11.5 scientific state is:

```text
RH – q
RH – e
q  – e
r  – q
```

For each pair, the engine maintains a pair-specific state that supports pair count, paired means, paired second-moment terms, cross-product accumulation, covariance and correlation.

This matters because two variables can have different missing-value patterns. Reusing a marginal statistic derived from a different support set can produce a mathematically inconsistent correlation.

---

## 13. Empirical RH-q joint product

v11.5 maintains an empirical joint product for RH and q.

The default histogram support is described by configured physical ranges and an 8 x 8 two-dimensional resolution. The histogram is a separate statistical layer from the 33 temporal bins.

The current configured histogram levels are:

```text
L1
L2
```

The distinction is therefore:

```text
Temporal statistical state: 33 bins
Empirical RH-q histogram: 8 x 8 cells
```

The histogram is intended as an empirical reference representation of the paired sample support. It should not be silently replaced by a parametric normality assumption.

---

## 14. Missing values, support and validity

The engine distinguishes between:

- an actually missing/invalid input;
- a valid observation that happens to lie outside a configured histogram support;
- a physically invalid thermodynamic transformation;
- a valid but supersaturated state that is retained according to the configured scientific rules.

Support accounting is central to interpretation. A percentile or threshold count without a clear valid-sample definition can be misleading.

The appropriate downstream workflow is therefore to inspect sample support before interpreting changes as climate signals.

---

## 15. Water, land and all-NaN regions

Spatial locations with no valid ERA5-Land land-surface samples can naturally produce all-NaN slices during intermediate min/max candidate calculations.

An `All-NaN slice encountered` warning by itself does not establish that the scientific computation failed. Such warnings are expected at locations that have no valid support for a given operation.

A release-quality workflow should, however, distinguish harmless support-driven warnings from unexpected widespread data loss. Final validation and audit should therefore inspect support counts and output completeness rather than treating every warning as an exception.

---

## 16. Decade routing and FULL product

The default climatology period is:

```text
1981-2020
```

The engine routes each source year into its corresponding decade product and into the full climatology state.

The standard decade products are:

```text
DECADE_1981_1990
DECADE_1991_2000
DECADE_2001_2010
DECADE_2011_2020
```

The full product is:

```text
FULL_1981_2020
```

The purpose of this routing design is to allow decade comparison and a full-period reference while preserving mergeability and auditability.

---

## 17. Spatial chunking

Long climatology runs operate on spatial blocks rather than loading the complete grid into memory at once.

The public defaults are:

```text
chunk_lat = 64
chunk_lon = 128
```

Chunk sizes control the trade-off between memory use, I/O granularity and processing efficiency. They are operational parameters, not scientific parameters, provided the exact executed configuration is retained in provenance.

A production installation should benchmark representative blocks before committing to a multidecadal run if hardware or storage characteristics differ materially from the development environment.

---

## 18. Compression and storage

The default NetCDF configuration uses:

```text
compression level = 4
shuffle           = True
```

Compression settings are intended to reduce storage while keeping the checkpoint/output representation manageable. They should not be changed silently in the middle of a production series.

---

## 19. Checkpointing

Long-running climatology jobs require durable state.

The engine writes checkpoint state for spatial blocks and statistical periods. The purpose is to make the calculation resumable without restarting the entire climatology from raw data after every interruption.

The important operational rule is:

> **Durable state must become authoritative before progress is treated as committed.**

The checkpoint layer is paired with a journal/transaction concept so that partially completed updates can be distinguished from committed state.

The current checkpoint implementation also carries its own technical checkpoint identifier. That identifier is an internal storage contract and should not be casually rewritten simply because the public software release is v11.5.

---

## 20. Transactions and COMMITTED truth

The write lifecycle follows the conceptual pattern:

```text
prepare
  ↓
verify
  ↓
commit
  ↓
COMMITTED
```

The `COMMITTED` state is the authoritative operational truth for whether a work unit is complete.

If a write is interrupted before commit, the engine should treat that work unit as incomplete and recover/recompute it according to the journal and checkpoint rules.

This is deliberately stronger than inferring completion from the presence of a partially written output file.

---

## 21. Before-image recovery

The transaction layer maintains enough information to reconstruct the previous state of an affected output region when an interruption occurs during a transactional update.

The practical objective is not zero recomputation at all costs. The objective is **bounded, controlled recovery without accepting corrupted statistical state as truth**.

A power-loss or process interruption must therefore be treated as an operational event, not a reason to manually patch numerical cells.

---

## 22. Progress reporting

Long calculations should produce explicit progress messages. The public engine provides structured logging so that the operator can see where the process is and whether an individual work unit completed.

The pilot path goes further and can report each processed day.

Typical pilot messages include:

```text
PILOT 1/9 | loading grid
PILOT 2/9 | opening YYYY-MM
PILOT 3/9 | validating input triplet
PILOT 4/9 | processing first spatial block
PILOT DAY 01/31 | YYYY-MM-01 | start
PILOT DAY 01/31 | YYYY-MM-01 | PASS
...
PILOT 9/9 | PASS
```

These messages are useful for human monitoring as well as archived provenance.

---

## 23. Pilot mode

Pilot mode is the recommended first real-data test before a long production run.

The v11.5 pilot:

- loads the configured grid;
- opens one real T2m/D2m/SP monthly triplet;
- validates the triplet;
- processes the first spatial block;
- when no explicit day is supplied, processes the complete requested month;
- emits progress for each day;
- flushes the checkpoint state;
- reopens the checkpoint for structural verification;
- reports an explicit PASS/FAIL result.

This makes pilot mode much more informative than a synthetic self-test while remaining far smaller than a multidecade production run.

### Example: one-day pilot

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "pilot --year 2011 --month 1 --day 1 --verbose"
```

### Example: one-month pilot

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "pilot --year 2011 --month 1 --verbose"
```

When the day argument is omitted or zero, the monthly pilot processes every day in that month.

---

## 24. Real-data validation milestone

A real January 2011 pilot was completed successfully on one spatial block.

Observed result:

```text
31/31 days PASS
Checkpoint flush PASS
Checkpoint reopen PASS
Overall pilot PASS
```

The run demonstrated that the real ERA5-Land processing path, daily transactions, state accumulation and checkpoint reopening could complete for a full month without the earlier state-shape failures encountered during development.

This is evidence for the tested execution path. It is not, by itself, a substitute for a complete multidecadal production run and final audit.

---

## 25. Command-line interface

The public command set is:

```text
selftest
validate-input
pilot
run
audit
merge-audit
report
benchmark
```

### 25.1 selftest

Runs deterministic internal tests for calendar handling, moment accumulation, merging, physical derivation, state updates and transaction behavior.

Example:

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "selftest --verbose"
```

### 25.2 validate-input

Validates the twelve monthly input triplets for a target year. It is intentionally read-only and does not perform climatological accumulation.

Example:

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "validate-input --year 2011 --verbose"
```

### 25.3 pilot

Runs a bounded real-data pilot against the first spatial block.

Example:

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "pilot --year 2011 --month 1 --verbose"
```

### 25.4 run

Runs the configured production accumulation.

Example:

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "run --verbose"
```

The production command is guarded by a workspace lock so two independent production processes do not mutate the same production workspace simultaneously.

### 25.5 audit

Audits the generated grid/state products according to the engine's audit rules.

Example:

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "audit"
```

### 25.6 merge-audit

Checks merge consistency between compatible product states.

Example:

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "merge-audit"
```

### 25.7 report

Generates the configured reporting outputs from the accumulated product.

Example:

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "report"
```

### 25.8 benchmark

Measures representative processing characteristics so storage, CPU, RAM and chunk settings can be evaluated before a long run.

Example:

```python
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "benchmark"
```

---

## 26. Python environment

The project is intended for a modern Python environment with the scientific dependencies listed in the repository environment files.

The runtime used during v11.5 validation included Python 3.12.x and a NumPy/netCDF4-based scientific stack.

The exact environment used for a publishable result should be archived together with:

- source code;
- environment definition;
- dependency versions;
- executed configuration;
- input inventory;
- script SHA-256;
- output hashes.

---

## 27. Recommended production workflow

A reliable production workflow is:

### Phase A — source preparation

1. Freeze the source file identified as v11.5.
2. Do not mix old checkpoints from incompatible schema revisions.
3. Confirm the intended output root.
4. Preserve the release documentation alongside the code.

### Phase B — environment

5. Activate the intended conda/environment installation.
6. Confirm Python and key scientific package versions.
7. Run `selftest`.

### Phase C — inputs

8. Confirm T2m, D2m and SP inventories.
9. Confirm monthly completeness for the requested years.
10. Validate time axes and grids.
11. Confirm units.

### Phase D — real-data pilot

12. Run a one-day pilot.
13. Run a full-month pilot.
14. Inspect logs for grid alignment, warnings and support.
15. Confirm checkpoint reopen PASS.

### Phase E — production

16. Start the production `run` only after the pilot is satisfactory.
17. Keep the terminal log.
18. Monitor progress and output growth.
19. Do not manually edit checkpoint NetCDF values.
20. Allow the transaction system to determine committed state.

### Phase F — post-processing

21. Run `audit`.
22. Run `merge-audit`.
23. Generate `report` outputs.
24. Archive manifests and hashes.
25. Preserve the exact code and configuration used for the run.

---

## 28. Parallel execution guidance

The scientific state is spatially and temporally chunked, but operators should not assume that arbitrary simultaneous commands are safe merely because the machine has multiple CPU cores.

The production workspace is a stateful scientific artifact. Two processes must not independently mutate the same production checkpoint tree at the same time.

The engine therefore distinguishes between:

- **independent pilot workspaces**, which can be isolated by year/month/run identifier;
- **the shared production workspace**, which is protected by an atomic workspace lock.

The safe mental model is:

```text
Pilot A → isolated pilot workspace
Pilot B → isolated pilot workspace

Production A → production workspace lock
Production B → waits/fails rather than corrupting shared state
```

Parallelism should be introduced only where the output contract and checkpoint design explicitly support it.

---

## 29. Directory expectations

A typical installation contains:

```text
HumidClimatologyEngine/
├── humid_climatology_engine_v11.5.py
├── README.md
├── CHANGELOG.md
├── CITATION.cff
├── CONTRIBUTING.md
├── environment.yml
├── LICENSE
├── MANIFEST.txt
├── PROJECT_CONTENTS.txt
├── pyproject.toml
├── RELEASE_MANIFEST.txt
├── requirements.txt
├── SECURITY.md
├── .gitignore
└── docs/
    ├── README_DOCS.md
    ├── HumidClimatologyEngine_v11.5_Comprehensive_Scientific_Engineering_Reference_FINAL.docx
    ├── HumidClimatologyEngine_v11.5_Comprehensive_Scientific_Engineering_Reference_FINAL.pdf
    ├── HumidClimatologyEngine_v8_vs_v11.5_Complete_Comparison.docx
    ├── HumidClimatologyEngine_v11.5_User_Guide_and_Production_Runbook.docx
    ├── HumidClimatologyEngine_v11.5_Who_Should_Use_It_and_Applications.docx
    └── HumidClimatologyEngine_v11.5_Analytical_Graphical_Comparison_Toolkit_Specification.docx
```

The repository may also contain tests, development utilities, notebooks, manifests and archived historical material.

---

## 30. Configuration reference

The principal operational defaults in the public source are conceptually:

```text
start_year = 1981
end_year   = 2020
chunk_lat  = 64
chunk_lon  = 128
compression = 4
shuffle = True
```

The input directories and output root are described earlier in this document.

The histogram layer is configured from the declared histogram levels and pair set. The configuration object is validated before long processing begins.

---

## 31. Output interpretation

The output is not one undifferentiated table. It is a collection of statistical products with different semantics.

When interpreting a field, ask:

1. Which period does it belong to?
2. Which temporal level does it represent?
3. Which variable or pair does it represent?
4. What is its valid sample support?
5. Which configuration and schema produced it?
6. Is it a scalar state, a paired state or an empirical histogram?
7. Is the value a raw count, a moment, a finalized statistic or a derived diagnostic?

This discipline prevents accidental comparison of quantities with different support or temporal meaning.

---

## 32. Reproducibility contract

A published result should be reconstructable from the frozen source and its provenance, not from memory.

At minimum, archive:

```text
software version
schema version
configuration values
configuration fingerprint
source script hash
input inventory
calendar convention
chunk configuration
checkpoint lineage
output inventory
final output hashes
run timestamps
```

The exact executed configuration is more important than a generic README statement that a parameter "usually" has a certain value.

---

## 33. Statistical mergeability

A key architectural property of v11.5 is that the statistical state can be merged when the two states are compatible.

This enables:

- separate decade accumulation;
- spatial block processing;
- controlled recovery;
- reconstruction of full-period states from component products;
- validation of merge equivalence.

Mergeability is not a license to combine incompatible schemas. Compatibility must include the meaning of the variables, bins, thresholds, calendar rules and statistical state layout.

---

## 34. Calendar handling

The calendar model is part of the scientific contract.

A climatological day index must be interpreted consistently across leap and non-leap years. The engine's calendar tests exist to prevent off-by-one and leap-day errors from entering a multidecadal accumulation.

Calendar correctness should therefore be verified in `selftest` and not inferred from successful file I/O.

---

## 35. Diagnostics and numerical warnings

Numerical warnings should be classified.

### Expected/support-driven examples

- all-NaN reductions at locations with no valid source support;
- undefined higher-order statistics where sample count is insufficient;
- empty histogram support for a location/time combination.

### Potentially serious examples

- unexpected broadcasting errors;
- shape mismatches between state arrays and input blocks;
- checkpoint schema mismatch during a new release run;
- corruption detected during checkpoint reopen;
- transaction or journal disagreement;
- input grids with non-equivalent coordinate sets.

The correct reaction to a potentially serious warning is to stop and inspect the first failing work unit, not to delete the evidence and restart.

---

## 36. Troubleshooting patterns

### Error: schema mismatch

Cause: an existing checkpoint was produced under a different schema contract.

Response: do not mix the checkpoint with the new schema. Use the correct release/checkpoint lineage or start a clean compatible workspace.

### Error: index out of bounds in state arrays

Cause: a mismatch between the number of state bins and the bin indices being addressed.

Response: treat this as an implementation bug, not an input-data peculiarity. Verify that the 33 temporal bins and any histogram dimensions are being handled separately.

### Error: broadcasting failure

Cause: input data block shape and validity mask shape disagree.

Response: inspect the first failing variable and print/record its shape at the block boundary before changing the statistical code.

### Error: checkpoint already locked

Cause: another production process is using the same shared workspace.

Response: do not start a second production mutator against that workspace.

### Warning: all-NaN slice

Cause: no valid data at the affected location/support.

Response: inspect support counts and land/sea coverage. An isolated warning can be expected; widespread unexpected warnings require investigation.

---

## 37. What should never be done manually

Do not:

- edit checkpoint NetCDF values by hand;
- delete journal records to force completion;
- rename an incompatible checkpoint into a new version folder;
- copy a partial output file into a production state directory;
- change a threshold while preserving the old configuration fingerprint;
- alter bin definitions without updating the schema/documentation;
- combine outputs from different scientific contracts merely because variable names look identical.

Scientific reproducibility is more important than making a run appear complete.

---

## 38. Analytical companion tools

The core engine should remain separate from expensive exploratory analysis whenever practical.

The core accumulation produces frozen statistical products. Analytical tools can then consume those outputs for:

- decadal comparisons;
- spatial maps;
- annual and diurnal summaries;
- threshold diagnostics;
- distribution comparison;
- model-selection summaries;
- publication graphics.

This separation allows analytical methods to evolve without silently changing the underlying climatological calculation.

---

## 39. Recommended research workflow

A researcher using v11.5 should generally proceed as follows:

```text
Define scientific question
        ↓
Confirm suitable ERA5-Land period
        ↓
Validate input inventory
        ↓
Run selftest
        ↓
Run one-day real-data pilot
        ↓
Run one-month real-data pilot
        ↓
Inspect support and sanity
        ↓
Run production climatology
        ↓
Audit
        ↓
Merge-audit
        ↓
Generate reports
        ↓
Interpret with support/provenance
        ↓
Archive code + config + hashes
```

The staged approach is deliberate. It is cheaper to discover a structural problem after one day than after several weeks of multidecadal processing.

---

## 40. Example production checklist

Before starting:

- [ ] Correct v11.5 source selected.
- [ ] Python environment confirmed.
- [ ] Self-test PASS.
- [ ] Target year inventory complete.
- [ ] T2m/D2m/SP triplets present.
- [ ] Time axes valid.
- [ ] Grids compatible.
- [ ] Units valid.
- [ ] Output root clean or intentionally resumable.
- [ ] Old incompatible checkpoints excluded.
- [ ] Pilot PASS.

During the run:

- [ ] Progress is advancing.
- [ ] No unexplained shape errors.
- [ ] Checkpoints are being updated.
- [ ] Journal/commit state remains coherent.
- [ ] Disk space is sufficient.
- [ ] Logs are retained.

After completion:

- [ ] Production run completed.
- [ ] Audit PASS.
- [ ] Merge-audit PASS.
- [ ] Report generated.
- [ ] Output inventory archived.
- [ ] Source hash archived.
- [ ] Configuration archived.
- [ ] Input inventory archived.
- [ ] Final hashes archived.

---

## 41. Release documentation set

The public documentation package should contain:

1. The comprehensive scientific and engineering reference.
2. The v8 versus v11.5 comparison.
3. The user guide and production runbook.
4. The audience/applications guide.
5. The analytical and graphical comparison toolkit specification.
6. The long-form README.
7. CHANGELOG and citation metadata.
8. Release and project manifests.

The documentation index should identify v11.5 as the current public release identity and distinguish historical material from active release documentation.

---

## 42. Citation and provenance

Any scientific publication using v11.5 should cite the software release and document the exact period, configuration and source data used.

At minimum, a reproducible methods statement should identify:

- HumidClimatologyEngine v11.5;
- ERA5-Land;
- T2m, D2m and surface pressure inputs;
- target climatological period;
- temporal state definition;
- histogram configuration if used;
- threshold configuration if used;
- output schema;
- validation/pilot status;
- analysis software used downstream.

---

## 43. Interpretation of v11.5 relative to station data

ERA5-Land provides spatially complete model/reanalysis fields, while station records provide local measurements with different coverage and error characteristics.

When a study has station observations, the engine's climatology should not be presented as an automatic replacement for station analysis. Instead, station and reanalysis products can be compared as complementary evidence.

Disagreement should be investigated rather than assumed to be a software error.

---

## 44. Interpretation of dependence

Correlation and covariance are properties of a specified support set and time aggregation.

A change in pair correlation can result from:

- a change in the underlying dependence structure;
- a change in marginal support;
- a change in valid-pair membership;
- a change in temporal aggregation;
- sampling limitations.

Therefore a correlation map should always be accompanied by sufficient support information and a clear statement of the temporal level used.

---

## 45. Interpretation of thresholds

Threshold frequencies should be stated with their denominator definition.

For example, a high-RH count should be interpreted as a count among samples that were valid for the RH statistic at the specified temporal support, not necessarily among every nominal timestamp in the source archive.

Joint thresholds require the intersection of the pair-valid masks unless the product explicitly defines another support convention.

---

## 46. Interpretation of minimum and maximum

Minimum and maximum values are sensitive to sample coverage, valid support and rare events. They should therefore be interpreted alongside sample count and period definition.

A larger maximum in one decade is not automatically evidence of a physically stronger extreme unless the support and comparability of the samples are understood.

---

## 47. Why the 33-bin model matters

If the daily climatology were the only state retained, diurnal structure could be erased.

v11.5 therefore keeps:

```text
1 daily pooled bin
8 three-hour bins
24 hourly bins
```

This lets downstream analysis answer questions at progressively finer temporal resolution without requiring the raw multidecadal archive to be replayed merely to recover standard temporal summaries.

---

## 48. Why the empirical histogram is separate

The empirical RH-q histogram answers a different question from the temporal state dimension.

The 33 temporal bins say **when** the sample belongs in the statistical accumulation.

The 8 x 8 RH-q histogram says **where within the configured joint physical support** the paired samples occur.

Keeping these dimensions separate avoids the common error of treating every numerical bin count as if it represented the same statistical object.

---

## 49. Performance expectations

The dominant cost of a multidecadal ERA5-Land calculation is expected to come from:

- repeated NetCDF access;
- spatial block processing;
- nonlinear moisture transformations;
- statistical accumulation;
- checkpoint writes;
- joint counting.

Wall-clock performance is therefore a systems property as much as a CPU property.

A representative pilot provides a much better basis for estimating a production run than a synthetic microbenchmark alone.

The January 2011 one-block pilot completed 31 days in roughly 17.5 minutes in the validation environment, with daily processing stabilizing around the mid-30-second range. This is an empirical benchmark for that environment and block size, not a universal promise for other hardware.

---

## 50. Memory philosophy

The engine's statistical state is designed so that the entire raw multidecadal archive does not need to remain in RAM.

Memory pressure is controlled primarily by:

- spatial chunk size;
- temporary transformed arrays;
- checkpoint block size;
- concurrent process count.

Reducing chunk dimensions can lower peak memory at the cost of more I/O overhead. Increasing them may improve throughput until memory or cache behavior becomes limiting.

---

## 51. Release synchronization principle

The code, schema and documentation are one scientific release unit.

Whenever any of the following changes materially:

- a formula;
- a variable definition;
- a threshold;
- a bin definition;
- a calendar rule;
- an output field;
- a validity rule;
- a checkpoint schema;

the corresponding documentation and release metadata must be regenerated together.

A version number is not just a filename. It is a statement about the reproducibility contract.

---

## 52. v11.5 public release statement

HumidClimatologyEngine v11.5 is the public release represented by this README and its accompanying documentation package.

The release is intended to provide a reproducible, empirical, hourly-input moisture climatology architecture with:

- explicit grid alignment;
- observation-level moisture derivation;
- L1/L2/L3 temporal states;
- mergeable statistical moments;
- extrema and threshold counts;
- pair-specific dependence state;
- empirical RH-q joint support;
- decade and full-period routing;
- checkpointed long-running execution;
- transactional completion truth;
- audit and merge-audit tools;
- structured pilot validation;
- detailed provenance and documentation.

---

## 53. Historical and release discipline

The project history contains earlier versions and intermediate engineering changes. Those records are valuable for development and scientific reproducibility.

The current public project identity is nevertheless unambiguous:

> **HumidClimatologyEngine v11.5**

Historical implementation notes must not be mistaken for a different public scientific method when they describe engineering maintenance around the release.

---

## 54. Final operational rule

The safest production principle is:

> **Validate first. Pilot on real data. Run in a controlled workspace. Trust committed state. Audit after completion. Archive the exact provenance.**

That rule is more important than any single command-line shortcut.

---

## 55. Quick command summary

```text
# Deterministic tests
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "selftest --verbose"

# Validate one real year
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "validate-input --year 2011 --verbose"

# One-day real-data pilot
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "pilot --year 2011 --month 1 --day 1 --verbose"

# One-month real-data pilot
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "pilot --year 2011 --month 1 --verbose"

# Full production run
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "run --verbose"

# Post-run audit
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "audit"

# Merge audit
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "merge-audit"

# Reporting
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "report"

# Benchmark
%runfile "C:/Users/AminFazlKazemi/Downloads/humid_climatology_engine_v11.5.py" --args "benchmark"
```

---

## 56. Documentation map

For the full release documentation, use:

- `HumidClimatologyEngine_v11.5_Comprehensive_Scientific_Engineering_Reference_FINAL.docx`
- `HumidClimatologyEngine_v11.5_Comprehensive_Scientific_Engineering_Reference_FINAL.pdf`
- `HumidClimatologyEngine_v8_vs_v11.5_Complete_Comparison.docx`
- `HumidClimatologyEngine_v11.5_User_Guide_and_Production_Runbook.docx`
- `HumidClimatologyEngine_v11.5_Who_Should_Use_It_and_Applications.docx`
- `HumidClimatologyEngine_v11.5_Analytical_Graphical_Comparison_Toolkit_Specification.docx`
- `README_DOCS.md`

Together these files cover the scientific contract, engineering design, historical comparison, operational runbook, intended users, applications, analytical layer and release navigation.

---

## 57. Closing note

HumidClimatologyEngine v11.5 is best understood as a reproducible statistical infrastructure for climate-moisture analysis rather than as a single formula or a single report.

Its scientific value depends on keeping four things together:

```text
observations
+ transformations
+ statistical state
+ provenance
```

The release is strongest when those four remain synchronized from the initial input inventory to the final published figure.

**HumidClimatologyEngine v11.5 — public release.**
