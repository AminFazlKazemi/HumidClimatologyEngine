# HumidClimatologyEngine — Scientific Methodology

## 1. Objective

Estimate a 366-day, grid-cell moisture climatology from ERA5-Land using a
joint model of `(T, Td, logP)` and nonlinear Monte Carlo propagation.

## 2. Input sample

A sample is paired-valid only when T, Td and P are finite and P > 0.

## 3. Calendar

Feb 28 and Feb 29 are pooled into DOY 60; DOY 59 is reserved; Mar 1 maps to
DOY 61.

## 4. Statistics

The annual engine stores Welford means, M2 states and cross moments. Annual
states are merged using the standard parallel/merge formulas.

## 5. Distribution

The daily joint model is a trivariate normal representation of `(T, Td, logP)`
with empirical means, standard deviations and correlations.

## 6. Physical layer

Phase-aware saturation vapor pressure is used for T and Td. RH is calculated
as `100*es(Td)/es(T)`, vapor pressure is `es(Td)`, mixing ratio is
`0.622*e/(P-e)`, and specific humidity is `r/(1+r)`.

## 7. Monte Carlo

Sampling is chunked by grid cells and batched by sample count. The full
sample tensor is never materialized.

## 8. Higher moments

Mean, sample standard deviation, bias-corrected skewness and Fisher excess
kurtosis are retained for the transformed moisture variables.

## 9. Validation

Use leap-day tests, independent physics Ground Truth, synthetic MVN tests,
Welford/Pébay comparisons, physical bounds, convergence experiments, and
checkpoint checksums.

## 10. Limitations

The central limitations are the multivariate-normal assumption, finite Monte
Carlo error, reanalysis uncertainty, thermodynamic approximation choice, and
the custom leap-day convention.

## 11. Recommended sensitivity matrix

| Dimension | Alternatives |
|---|---|
| Monte Carlo | 500 / 1000 / 2000 / 5000 / 10000 |
| Pressure model | P / logP |
| Saturation | water-only / phase-aware |
| Calendar | pooled / alternative Feb treatment |
| Joint model | Gaussian / mixture / copula |

The final method should be selected based on quantitative validation, not on
runtime alone.
