#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Visual reporting for per-DOY x grid-cell bivariate model selection.

Expected NetCDF variables:
  best_model_code(doy, latitude, longitude): int16 codes
  model_names as a global attribute, JSON object {"1":"...", ...}
  optional pair, month, day coordinates

Produces:
  - PNG spatial dominance maps by month
  - PNG monthly dominance share
  - PNG DOY x model occurrence heatmap
  - PDF report combining the figures and QA summary

This reporter is intentionally agnostic to the statistical family. It does not
assume Gaussian copulas or any fixed bivariate family.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import xarray as xr


def _month_from_doy(ds: xr.Dataset) -> np.ndarray:
    if "month" in ds:
        return np.asarray(ds["month"].values, dtype=int)
    # fallback for the project's 366-slot calendar
    month = np.full(ds.dims["doy"], -1, dtype=int)
    d = np.datetime64("2000-01-01")
    # Slot 59 is reserved; slot 60 is Feb-28/29 composite.
    for i in range(58):
        month[i] = int(str(d)[5:7])
        d = d + np.timedelta64(1, "D")
    month[58] = -1
    month[59] = 2
    d = np.datetime64("2000-03-01")
    for i in range(60, 366):
        month[i] = int(str(d)[5:7])
        d = d + np.timedelta64(1, "D")
    return month


def load_selection(path: Path) -> tuple[xr.Dataset, dict[int, str]]:
    ds = xr.open_dataset(path)
    if "best_model_code" not in ds:
        raise ValueError("Expected variable 'best_model_code'.")
    raw = ds.attrs.get("model_names", "{}")
    if isinstance(raw, str):
        names = json.loads(raw)
    else:
        names = raw
    mapping = {int(k): str(v) for k, v in names.items()}
    return ds, mapping


def build_report(path: Path, output_pdf: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    ds, model_names = load_selection(path)
    codes = np.asarray(ds["best_model_code"].values)
    lat = np.asarray(ds["latitude"].values)
    lon = np.asarray(ds["longitude"].values)
    months = _month_from_doy(ds)
    unique_codes = sorted(int(x) for x in np.unique(codes[np.isfinite(codes)]))

    generated = []
    with PdfPages(output_pdf) as pdf:
        fig = plt.figure(figsize=(11.7, 8.3))
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.text(0.02, 0.92, "Bivariate Distribution Dominance Report", fontsize=22, fontweight="bold")
        ax.text(0.02, 0.85, f"Selection product: {path.name}", fontsize=11)
        ax.text(0.02, 0.80, "Dominance is computed from the selected model code for each DOY x grid cell.", fontsize=11)
        ax.text(0.02, 0.74, "No distribution family is assumed a priori; the report visualizes the fitted selection result.", fontsize=11)
        y = 0.63
        for code in unique_codes:
            ax.text(0.04, y, f"{code}: {model_names.get(code, 'UNKNOWN')}", fontsize=12)
            y -= 0.045
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Global monthly shares
        share_rows = []
        for m in range(1, 13):
            mask = months == m
            vals = codes[mask]
            vals = vals[np.isfinite(vals)]
            row = {code: float(np.mean(vals == code)) if vals.size else np.nan for code in unique_codes}
            share_rows.append(row)
        fig, ax = plt.subplots(figsize=(11.7, 8.3))
        x = np.arange(1, 13)
        bottom = np.zeros(12)
        for code in unique_codes:
            vals = np.array([r[code] for r in share_rows], dtype=float)
            ax.bar(x, vals, bottom=bottom, label=model_names.get(code, str(code)))
            bottom += np.nan_to_num(vals)
        ax.set_xlabel("Month")
        ax.set_ylabel("Fraction of selected grid-day cells")
        ax.set_title("Monthly dominance of bivariate model families")
        ax.set_xticks(x)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=8, ncol=2)
        pdf.savefig(fig, bbox_inches="tight")
        fig.savefig(output_dir / "monthly_dominance_share.png", dpi=180, bbox_inches="tight")
        generated.append(output_dir / "monthly_dominance_share.png")
        plt.close(fig)

        # DOY x model occurrence heatmap
        fig, ax = plt.subplots(figsize=(11.7, 8.3))
        occ = np.zeros((codes.shape[0], len(unique_codes)), dtype=float)
        for i, code in enumerate(unique_codes):
            occ[:, i] = np.mean(codes == code, axis=(1, 2))
        im = ax.imshow(occ.T, aspect="auto", interpolation="nearest", origin="lower")
        ax.set_xlabel("Climatological DOY")
        ax.set_ylabel("Model")
        ax.set_yticks(np.arange(len(unique_codes)))
        ax.set_yticklabels([model_names.get(c, str(c)) for c in unique_codes])
        ax.set_title("DOY-wise spatial dominance fraction")
        fig.colorbar(im, ax=ax, label="Fraction of grid cells")
        pdf.savefig(fig, bbox_inches="tight")
        fig.savefig(output_dir / "doy_model_dominance.png", dpi=180, bbox_inches="tight")
        generated.append(output_dir / "doy_model_dominance.png")
        plt.close(fig)

        # Monthly spatial maps: winner for each cell
        for month in range(1, 13):
            doy_idx = np.flatnonzero(months == month)
            if doy_idx.size == 0:
                continue
            sub = codes[doy_idx]
            # mode over the DOYs belonging to the month
            winner = np.full(sub.shape[1:], -1, dtype=np.int16)
            for j in range(sub.shape[1]):
                for i in range(sub.shape[2]):
                    vals = sub[:, j, i]
                    vals = vals[np.isfinite(vals)]
                    if vals.size:
                        counts = [(int(np.sum(vals == c)), c) for c in unique_codes]
                        winner[j, i] = max(counts)[1]
            fig, ax = plt.subplots(figsize=(11.7, 8.3))
            im = ax.imshow(winner, origin="upper", aspect="equal",
                           extent=(lon.min(), lon.max(), lat.min(), lat.max()))
            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.set_title(f"Dominant bivariate model by location - month {month:02d}")
            fig.colorbar(im, ax=ax, label="Model code")
            pdf.savefig(fig, bbox_inches="tight")
            out = output_dir / f"dominant_model_month_{month:02d}.png"
            fig.savefig(out, dpi=180, bbox_inches="tight")
            generated.append(out)
            plt.close(fig)

    ds.close()
    print(f"Report written: {output_pdf}")
    for p in generated:
        print(p)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("selection_netcdf", type=Path)
    p.add_argument("--output-pdf", type=Path, default=Path("bivariate_distribution_dominance_report.pdf"))
    p.add_argument("--output-dir", type=Path, default=Path("bivariate_dominance_figures"))
    args = p.parse_args()
    build_report(args.selection_netcdf, args.output_pdf, args.output_dir)


if __name__ == "__main__":
    main()
