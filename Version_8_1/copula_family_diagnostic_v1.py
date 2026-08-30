#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
copula_family_diagnostic_v1.py
================================

Scientific diagnostic for deciding whether one copula family is sufficient
for the moisture (RH, q) dependence structure.

Design
------
- Uses the SAME ERA5-Land folders and the SAME v8.2 window extractor.
- Uses one fixed decade (default 1991-2000).
- Samples a reproducible spatial grid of representative points.
- Tests 12 monthly climatological slots.
- For each point/month:
      1) extracts pooled 5-day observations over the decade,
      2) converts data to pseudo-observations,
      3) fits several copula families,
      4) compares log-likelihood, AIC and BIC.
- Produces per-case CSV + family winner summary + JSON report.

Copulas
-------
Gaussian
Student-t
Clayton
Gumbel
Frank

Interpretation
--------------
A single-family model is defensible when one family dominates most cases
and no seasonal/geographic subgroup systematically prefers another family.
A mixed copula strategy becomes scientifically justified when different
families repeatedly win in different regions/seasons with meaningful
AIC/BIC margins.

IMPORTANT:
This script is diagnostic only. It does NOT change the production runner
and does NOT silently replace the existing Gaussian-copula production logic.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln
from scipy.stats import multivariate_t, norm, rankdata, t as student_t

ENGINE_DIR = Path(r"K:\kazemi\papers\temperature_interpolation\HumidClimatologyEngine")
sys.path.insert(0, str(ENGINE_DIR))

# Import the v8.2 runner so the diagnostic uses exactly its ERA5 extractor,
# calendar implementation and physical RH/q derivation.
import moisture_copula_production_v82 as runner  # noqa: E402


# =============================================================================
# Configuration
# =============================================================================

DEFAULT_YEARS = list(range(1991, 2001))

# One representative climatological slot per month.
# These avoid reserved slot 59 and test the whole seasonal cycle.
MONTH_SLOTS = {
    1: 15,    # Jan-15
    2: 46,    # Feb-15
    3: 75,    # Mar-15
    4: 106,   # Apr-15
    5: 136,   # May-15
    6: 167,   # Jun-15
    7: 197,   # Jul-15
    8: 228,   # Aug-15
    9: 259,   # Sep-15
    10: 289,  # Oct-15
    11: 320,  # Nov-15
    12: 350,  # Dec-15
}

MIN_OBS = 200
EPS = 1e-10


# =============================================================================
# Stable pseudo-observations
# =============================================================================

def pseudo_obs(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if x.size < MIN_OBS:
        return np.empty(0, dtype=float)

    # Average ranks avoid duplicate-value problems.
    r = rankdata(x, method="average")
    return r / (x.size + 1.0)


def paired_pobs(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if x.size < MIN_OBS:
        return np.empty(0), np.empty(0)

    # Same row-wise mask is critical: copula dependence is paired.
    u = rankdata(x, method="average") / (x.size + 1.0)
    v = rankdata(y, method="average") / (y.size + 1.0)

    return np.clip(u, EPS, 1.0 - EPS), np.clip(v, EPS, 1.0 - EPS)


# =============================================================================
# Common metrics
# =============================================================================

def information_criteria(loglik: float, k: int, n: int) -> tuple[float, float]:
    aic = 2.0 * k - 2.0 * loglik
    bic = math.log(max(n, 1)) * k - 2.0 * loglik
    return float(aic), float(bic)


def safe_nll(value: float) -> float:
    if not np.isfinite(value):
        return 1e100
    return float(value)


# =============================================================================
# Gaussian copula
# =============================================================================

def gaussian_fit(u: np.ndarray, v: np.ndarray) -> dict:
    z1 = norm.ppf(u)
    z2 = norm.ppf(v)

    rho = float(np.corrcoef(z1, z2)[0, 1])
    rho = float(np.clip(rho, -0.999, 0.999))

    one_minus = max(1.0 - rho * rho, EPS)

    log_c = (
        -0.5 * np.log(one_minus)
        - (rho * rho * (z1 * z1 + z2 * z2) - 2.0 * rho * z1 * z2)
        / (2.0 * one_minus)
    )

    ll = float(np.sum(log_c))
    aic, bic = information_criteria(ll, 1, len(u))

    return {
        "family": "Gaussian",
        "parameters": {"rho": rho},
        "loglik": ll,
        "aic": aic,
        "bic": bic,
        "k": 1,
        "success": True,
    }


# =============================================================================
# Student-t copula
# =============================================================================

def t_copula_loglik(theta: np.ndarray, u: np.ndarray, v: np.ndarray) -> float:
    rho = float(np.tanh(theta[0]))
    df = float(2.0 + np.exp(theta[1]))

    x = student_t.ppf(u, df)
    y = student_t.ppf(v, df)

    cov = np.array([[1.0, rho], [rho, 1.0]], dtype=float)
    pts = np.column_stack([x, y])

    try:
        joint = multivariate_t.logpdf(pts, loc=np.zeros(2), shape=cov, df=df)
        marg1 = student_t.logpdf(x, df)
        marg2 = student_t.logpdf(y, df)
        ll = np.sum(joint - marg1 - marg2)
        if not np.isfinite(ll):
            return 1e100
        return float(-ll)
    except Exception:
        return 1e100


def t_copula_fit(u: np.ndarray, v: np.ndarray) -> dict:
    z1 = norm.ppf(u)
    z2 = norm.ppf(v)
    rho0 = float(np.clip(np.corrcoef(z1, z2)[0, 1], -0.95, 0.95))

    starts = [
        [np.arctanh(rho0), np.log(8.0 - 2.0)],
        [np.arctanh(rho0), np.log(20.0 - 2.0)],
        [0.0, np.log(10.0 - 2.0)],
    ]

    best = None
    for x0 in starts:
        res = minimize(
            t_copula_loglik,
            np.asarray(x0, dtype=float),
            args=(u, v),
            method="Nelder-Mead",
            options={"maxiter": 500, "xatol": 1e-5, "fatol": 1e-5},
        )

        if best is None or res.fun < best.fun:
            best = res

    if best is None or not np.isfinite(best.fun):
        return {
            "family": "Student-t",
            "parameters": {},
            "loglik": np.nan,
            "aic": np.nan,
            "bic": np.nan,
            "k": 2,
            "success": False,
        }

    rho = float(np.tanh(best.x[0]))
    df = float(2.0 + np.exp(best.x[1]))
    ll = -float(best.fun)
    aic, bic = information_criteria(ll, 2, len(u))

    return {
        "family": "Student-t",
        "parameters": {"rho": rho, "df": df},
        "loglik": ll,
        "aic": aic,
        "bic": bic,
        "k": 2,
        "success": bool(best.success),
    }


# =============================================================================
# Clayton copula
# =============================================================================

def clayton_loglik_theta(theta: float, u: np.ndarray, v: np.ndarray) -> float:
    if theta <= 0.0 or not np.isfinite(theta):
        return 1e100

    a = np.power(u, -theta)
    b = np.power(v, -theta)
    s = a + b - 1.0

    if np.any(s <= 0.0) or not np.all(np.isfinite(s)):
        return 1e100

    log_c = (
        np.log1p(theta)
        + (-theta - 1.0) * (np.log(u) + np.log(v))
        + (-2.0 - 1.0 / theta) * np.log(s)
    )

    if not np.all(np.isfinite(log_c)):
        return 1e100

    return float(-np.sum(log_c))


def clayton_fit(u: np.ndarray, v: np.ndarray) -> dict:
    res = minimize_scalar(
        clayton_loglik_theta,
        bounds=(1e-4, 50.0),
        args=(u, v),
        method="bounded",
        options={"xatol": 1e-6, "maxiter": 500},
    )

    if not np.isfinite(res.fun):
        return {
            "family": "Clayton",
            "parameters": {},
            "loglik": np.nan,
            "aic": np.nan,
            "bic": np.nan,
            "k": 1,
            "success": False,
        }

    theta = float(res.x)
    ll = -float(res.fun)
    aic, bic = information_criteria(ll, 1, len(u))

    return {
        "family": "Clayton",
        "parameters": {"theta": theta},
        "loglik": ll,
        "aic": aic,
        "bic": bic,
        "k": 1,
        "success": bool(res.success),
    }


# =============================================================================
# Gumbel copula
# =============================================================================

def gumbel_loglik_theta(theta: float, u: np.ndarray, v: np.ndarray) -> float:
    if theta < 1.0 or not np.isfinite(theta):
        return 1e100

    x = -np.log(u)
    y = -np.log(v)

    xt = np.power(x, theta)
    yt = np.power(y, theta)
    s = xt + yt

    if np.any(s <= 0.0) or not np.all(np.isfinite(s)):
        return 1e100

    s1 = np.power(s, 1.0 / theta)

    # c(u,v) =
    # exp(-s^(1/theta))
    # * x^(theta-1)y^(theta-1)
    # * s^(2/theta-2)
    # * (theta-1+s^(1/theta))/(uv)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        log_c = (
            -s1
            + np.log(theta - 1.0 + s1)
            + (theta - 1.0) * (np.log(x) + np.log(y))
            + (2.0 / theta - 2.0) * np.log(s)
            - np.log(u)
            - np.log(v)
        )

    if not np.all(np.isfinite(log_c)):
        return 1e100

    return float(-np.sum(log_c))


def gumbel_fit(u: np.ndarray, v: np.ndarray) -> dict:
    res = minimize_scalar(
        gumbel_loglik_theta,
        bounds=(1.0, 20.0),
        args=(u, v),
        method="bounded",
        options={"xatol": 1e-5, "maxiter": 500},
    )

    if not np.isfinite(res.fun):
        return {
            "family": "Gumbel",
            "parameters": {},
            "loglik": np.nan,
            "aic": np.nan,
            "bic": np.nan,
            "k": 1,
            "success": False,
        }

    theta = float(res.x)
    ll = -float(res.fun)
    aic, bic = information_criteria(ll, 1, len(u))

    return {
        "family": "Gumbel",
        "parameters": {"theta": theta},
        "loglik": ll,
        "aic": aic,
        "bic": bic,
        "k": 1,
        "success": bool(res.success),
    }


# =============================================================================
# Frank copula
# =============================================================================

def frank_loglik_theta(theta: float, u: np.ndarray, v: np.ndarray) -> float:
    if abs(theta) < 1e-5 or not np.isfinite(theta):
        return 1e100

    # Stable enough for the moderate dependence range used for diagnosis.
    a = np.exp(-theta * u)
    b = np.exp(-theta * v)
    e = math.exp(-theta)

    denom = 1.0 - e - (1.0 - a) * (1.0 - b)

    if np.any(denom <= 0.0) or not np.all(np.isfinite(denom)):
        return 1e100

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        log_c = (
            np.log(abs(theta))
            + np.log(abs(1.0 - e))
            - theta * (u + v)
            - 2.0 * np.log(denom)
        )

    if not np.all(np.isfinite(log_c)):
        return 1e100

    return float(-np.sum(log_c))


def frank_fit(u: np.ndarray, v: np.ndarray) -> dict:
    candidates = []

    for lo, hi in [
        (-30.0, -0.05),
        (0.05, 30.0),
    ]:
        res = minimize_scalar(
            frank_loglik_theta,
            bounds=(lo, hi),
            args=(u, v),
            method="bounded",
            options={"xatol": 1e-5, "maxiter": 500},
        )
        if np.isfinite(res.fun):
            candidates.append(res)

    if not candidates:
        return {
            "family": "Frank",
            "parameters": {},
            "loglik": np.nan,
            "aic": np.nan,
            "bic": np.nan,
            "k": 1,
            "success": False,
        }

    res = min(candidates, key=lambda r: r.fun)
    theta = float(res.x)
    ll = -float(res.fun)
    aic, bic = information_criteria(ll, 1, len(u))

    return {
        "family": "Frank",
        "parameters": {"theta": theta},
        "loglik": ll,
        "aic": aic,
        "bic": bic,
        "k": 1,
        "success": bool(res.success),
    }


# =============================================================================
# Family comparison
# =============================================================================

FITTERS = [
    gaussian_fit,
    t_copula_fit,
    clayton_fit,
    gumbel_fit,
    frank_fit,
]


def compare_copulas(rh: np.ndarray, q: np.ndarray) -> dict:
    u, v = paired_pobs(rh, q)

    if len(u) < MIN_OBS:
        return {
            "n": len(u),
            "results": [],
            "winner_bic": None,
            "winner_aic": None,
        }

    results = []

    for fitter in FITTERS:
        try:
            r = fitter(u, v)
        except Exception as exc:
            r = {
                "family": fitter.__name__,
                "parameters": {},
                "loglik": np.nan,
                "aic": np.nan,
                "bic": np.nan,
                "k": np.nan,
                "success": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
        results.append(r)

    valid = [r for r in results if np.isfinite(r.get("bic", np.nan))]

    winner_bic = min(valid, key=lambda r: r["bic"])["family"] if valid else None
    winner_aic = min(valid, key=lambda r: r["aic"])["family"] if valid else None

    if valid:
        best_bic = min(r["bic"] for r in valid)
        for r in valid:
            r["delta_bic"] = float(r["bic"] - best_bic)
    else:
        for r in results:
            r["delta_bic"] = np.nan

    return {
        "n": len(u),
        "results": results,
        "winner_bic": winner_bic,
        "winner_aic": winner_aic,
    }


# =============================================================================
# Spatial sample generation
# =============================================================================

def make_sample_points(
    lat: np.ndarray,
    lon: np.ndarray,
    n_lat: int = 4,
    n_lon: int = 4,
) -> list[dict]:
    if n_lat < 1 or n_lon < 1:
        raise ValueError("n_lat and n_lon must be >= 1")

    lat_positions = np.unique(
        np.round(np.linspace(0, len(lat) - 1, n_lat)).astype(int)
    )
    lon_positions = np.unique(
        np.round(np.linspace(0, len(lon) - 1, n_lon)).astype(int)
    )

    points = []
    pid = 0

    for j in lat_positions:
        for i in lon_positions:
            points.append({
                "point_id": pid,
                "lat_index": int(j),
                "lon_index": int(i),
                "latitude": float(lat[j]),
                "longitude": float(lon[i]),
            })
            pid += 1

    return points


# =============================================================================
# ERA5 extraction
# =============================================================================

def extract_point_decade(
    point: dict,
    slot: int,
    years: list[int],
    file_index,
) -> tuple[np.ndarray, np.ndarray]:

    # Configure the shared v8.2 worker state to one exact grid point.
    lat_idx = np.asarray([point["lat_index"]], dtype=np.int32)
    lon_idx = np.asarray([point["lon_index"]], dtype=np.int32)

    runner.worker_initializer(
        file_index,
        lat_idx,
        lon_idx,
        np.asarray([point["latitude"]], dtype=np.float32),
        np.asarray([point["longitude"]], dtype=np.float32),
        runner.ENGINE_CONFIG_HASH,
    )

    rh_parts = []
    q_parts = []

    try:
        for y in years:
            for target in runner.target_dates_for_slot(slot, [y]):
                rh, q = runner.extract_target_window_block(
                    target,
                    0, 1,
                    0, 1,
                )

                if rh.size:
                    rh_parts.append(rh[:, 0])
                    q_parts.append(q[:, 0])

    finally:
        runner.close_worker_cache()

    if not rh_parts:
        return np.empty(0), np.empty(0)

    rh_all = np.concatenate(rh_parts)
    q_all = np.concatenate(q_parts)

    return rh_all, q_all


# =============================================================================
# Main diagnostic
# =============================================================================

def run_diagnostic(
    years: list[int],
    n_lat: int,
    n_lon: int,
    output_dir: Path,
) -> None:

    output_dir.mkdir(parents=True, exist_ok=True)

    if len(years) < 5:
        raise ValueError("Use at least five years for a meaningful decade-scale diagnostic.")

    lat, lon = runner.get_grid()
    points = make_sample_points(lat, lon, n_lat=n_lat, n_lon=n_lon)

    file_index = runner.build_all_file_indices(years)

    cases = []
    detailed_rows = []

    total_cases = len(points) * len(MONTH_SLOTS)
    done_cases = 0

    print("=" * 88)
    print("COPULA FAMILY DIAGNOSTIC")
    print(f"Years       : {years[0]}-{years[-1]} ({len(years)} years)")
    print(f"Points      : {len(points)}")
    print(f"Months      : {len(MONTH_SLOTS)}")
    print(f"Cases       : {total_cases}")
    print("Families    : Gaussian | Student-t | Clayton | Gumbel | Frank")
    print("=" * 88)

    for point in points:
        for month, slot in MONTH_SLOTS.items():
            rh, q = extract_point_decade(
                point,
                slot,
                years,
                file_index,
            )

            comparison = compare_copulas(rh, q)

            case = {
                "point_id": point["point_id"],
                "latitude": point["latitude"],
                "longitude": point["longitude"],
                "lat_index": point["lat_index"],
                "lon_index": point["lon_index"],
                "month": month,
                "slot": slot,
                "n_obs": int(comparison["n"]),
                "winner_bic": comparison["winner_bic"],
                "winner_aic": comparison["winner_aic"],
            }
            cases.append(case)

            for r in comparison["results"]:
                detailed_rows.append({
                    **case,
                    "family": r.get("family"),
                    "loglik": r.get("loglik"),
                    "aic": r.get("aic"),
                    "bic": r.get("bic"),
                    "delta_bic": r.get("delta_bic"),
                    "k": r.get("k"),
                    "success": r.get("success"),
                    "parameters": json.dumps(
                        r.get("parameters", {}),
                        ensure_ascii=False,
                        sort_keys=True,
                    ),
                })

            done_cases += 1
            print(
                f"[{done_cases:4d}/{total_cases}] "
                f"point={point['point_id']:02d} "
                f"lat={point['latitude']:8.3f} "
                f"lon={point['longitude']:8.3f} "
                f"month={month:02d} "
                f"n={comparison['n']:4d} "
                f"BIC winner={comparison['winner_bic']}"
            )

            del rh, q, comparison
            gc.collect()

    case_df = pd.DataFrame(cases)
    detail_df = pd.DataFrame(detailed_rows)

    case_csv = output_dir / "copula_case_winners.csv"
    detail_csv = output_dir / "copula_family_scores.csv"

    case_df.to_csv(case_csv, index=False, encoding="utf-8-sig")
    detail_df.to_csv(detail_csv, index=False, encoding="utf-8-sig")

    # -------------------------------------------------------------------------
    # Winner summaries
    # -------------------------------------------------------------------------

    valid_cases = case_df.dropna(subset=["winner_bic"]).copy()

    winner_counts = (
        valid_cases["winner_bic"]
        .value_counts()
        .rename_axis("family")
        .reset_index(name="wins")
    )

    winner_counts["fraction"] = (
        winner_counts["wins"] / max(len(valid_cases), 1)
    )

    winner_counts_csv = output_dir / "copula_winner_summary.csv"
    winner_counts.to_csv(
        winner_counts_csv,
        index=False,
        encoding="utf-8-sig",
    )

    monthly = pd.crosstab(
        valid_cases["month"],
        valid_cases["winner_bic"],
    )

    monthly_csv = output_dir / "copula_winners_by_month.csv"
    monthly.to_csv(monthly_csv, encoding="utf-8-sig")

    # Spatial summary
    spatial = pd.crosstab(
        valid_cases["point_id"],
        valid_cases["winner_bic"],
    )

    spatial_csv = output_dir / "copula_winners_by_point.csv"
    spatial.to_csv(spatial_csv, encoding="utf-8-sig")

    # BIC margin: how often is the winner strongly preferred?
    # ΔBIC > 10 is conventionally interpreted as strong evidence relative to
    # the best competing model; here it is reported descriptively, not used as
    # an automatic scientific threshold for production.
    best_margin_rows = []

    for (pid, month), grp in detail_df.groupby(["point_id", "month"]):
        g = grp[np.isfinite(grp["bic"])].copy()

        if len(g) < 2:
            continue

        g = g.sort_values("bic")
        best = float(g.iloc[0]["bic"])
        second = float(g.iloc[1]["bic"])

        best_margin_rows.append({
            "point_id": int(pid),
            "month": int(month),
            "winner": str(g.iloc[0]["family"]),
            "best_bic": best,
            "second_bic": second,
            "delta_bic_2nd_minus_best": second - best,
        })

    margins = pd.DataFrame(best_margin_rows)
    margins_csv = output_dir / "copula_bic_margins.csv"
    margins.to_csv(margins_csv, index=False)

    strong = (
        int((margins["delta_bic_2nd_minus_best"] >= 10).sum())
        if not margins.empty else 0
    )

    report = {
        "diagnostic_version": "1.0",
        "years": years,
        "n_years": len(years),
        "n_points": len(points),
        "n_months": len(MONTH_SLOTS),
        "n_cases": int(len(case_df)),
        "n_valid_cases": int(len(valid_cases)),
        "copulas": [
            "Gaussian",
            "Student-t",
            "Clayton",
            "Gumbel",
            "Frank",
        ],
        "winner_counts_bic": {
            str(k): int(v)
            for k, v in valid_cases["winner_bic"].value_counts().items()
        },
        "winner_fractions_bic": {
            str(k): float(v)
            for k, v in winner_counts.set_index("family")["fraction"].items()
        },
        "cases_with_delta_bic_ge_10": strong,
        "total_margin_cases": int(len(margins)),
        "files": {
            "case_winners": str(case_csv),
            "family_scores": str(detail_csv),
            "winner_summary": str(winner_counts_csv),
            "winners_by_month": str(monthly_csv),
            "winners_by_point": str(spatial_csv),
            "bic_margins": str(margins_csv),
        },
    }

    report_path = output_dir / "copula_diagnostic_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("\n" + "=" * 88)
    print("FINAL DIAGNOSTIC SUMMARY")
    print("=" * 88)

    if winner_counts.empty:
        print("No valid cases were fitted.")
    else:
        print(winner_counts.to_string(index=False))

    print("\nBIC winner by month:")
    print(monthly.to_string())

    print(f"\nCases with second-best ΔBIC >= 10: {strong}/{len(margins)}")
    print(f"\nReport: {report_path}")
    print(f"Scores : {detail_csv}")
    print(f"Winners: {case_csv}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Diagnose whether multiple copula families are needed."
    )

    p.add_argument(
        "--years",
        type=str,
        default="1991-2000",
        help="Inclusive year range, e.g. 1991-2000",
    )

    p.add_argument(
        "--n-lat",
        type=int,
        default=4,
        help="Number of latitude sample positions.",
    )

    p.add_argument(
        "--n-lon",
        type=int,
        default=4,
        help="Number of longitude sample positions.",
    )

    p.add_argument(
        "--output",
        type=Path,
        default=Path(r"C:\c\copula_diagnostic_v1"),
        help="Diagnostic output directory.",
    )

    return p.parse_args()


def parse_year_range(spec: str) -> list[int]:
    if "-" not in spec:
        y = int(spec)
        return [y]

    a, b = spec.split("-", 1)
    y0 = int(a)
    y1 = int(b)

    if y0 > y1:
        y0, y1 = y1, y0

    return list(range(y0, y1 + 1))


def main() -> None:
    args = parse_args()
    years = parse_year_range(args.years)

    runner.ensure_env_for_worker()

    run_diagnostic(
        years=years,
        n_lat=args.n_lat,
        n_lon=args.n_lon,
        output_dir=args.output,
    )


if __name__ == "__main__":
    main()
