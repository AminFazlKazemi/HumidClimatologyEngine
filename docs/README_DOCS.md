# HumidClimatologyEngine v11.5 — Documentation Index

## Public release baseline

**HumidClimatologyEngine v11.5** is the public release baseline for this documentation package.

The documentation describes the production scientific/engineering contract, the ERA5-Land inputs,
temporal aggregation model, statistical state, paired dependence, empirical RH–q histogram,
checkpointing, recovery, audits, and downstream analytical products.

## Documentation set

1. `HumidClimatologyEngine_v11.5_Comprehensive_Scientific_Engineering_Reference_FINAL.docx`
   — Primary scientific and engineering reference.

2. `HumidClimatologyEngine_v11.5_Comprehensive_Scientific_Engineering_Reference_FINAL.pdf`
   — PDF release of the primary reference.

3. `HumidClimatologyEngine_v8_vs_v11.5_Complete_Comparison.docx`
   — Historical v8 versus v11.5 comparison and migration/regression guidance.

4. `HumidClimatologyEngine_v11.5_Who_Should_Use_It_and_Applications.docx`
   — Intended users, applications, research questions, interpretation boundaries, and non-goals.

5. `HumidClimatologyEngine_v11.5_User_Guide_and_Production_Runbook.docx`
   — Installation, input validation, pilot execution, production run, checkpoint/recovery,
   monitoring, audit, merge-audit, and troubleshooting.

6. `HumidClimatologyEngine_v11.5_Analytical_Graphical_Comparison_Toolkit_Specification.docx`
   — Downstream analytical, graphical, comparison, and reporting specification.

## Scientific baseline

The v11.5 state model contains three temporal levels:

- **L1** — one daily bin.
- **L2** — eight three-hour bins.
- **L3** — twenty-four hourly bins.

Together these form **33 temporal state bins**. The empirical RH–q histogram uses a separate
configured histogram state and must not be confused with the 33 temporal state bins.

The core variables are RH, vapor pressure `e`, mixing ratio `r`, and specific humidity `q`.
Configured paired products retain pair-specific valid masks and dependence state.

## Validation baseline

The documented v11.5 release path was exercised on real ERA5-Land data. The January 2011
single-block pilot processed **31/31 days successfully**, including checkpoint flush and
checkpoint reopen verification. These bounded tests validate the execution path but do not
replace the final multidecade production run and final audit.

## Legacy documentation

v10-named documents in Version_8_1/v10_audit/ are historical/legacy artifacts and should not be used as the current
public release reference after the v11.5 documentation set is installed.
