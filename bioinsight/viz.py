"""Visualization stage — Manhattan and QQ plots (implemented), plus planned plots.

Plots take a normalized results frame (snp, chr, bp, p) and optional threshold
lines. matplotlib is imported lazily so importing this module stays cheap.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from .registry import Stage, method, planned
from .stats import chrom_sort_key, genomic_inflation

_PALETTE = ["#2c7fb8", "#7fcdbb"]


def _lazy_plt():
    import matplotlib

    matplotlib.use("Agg")  # headless / no display required
    import matplotlib.pyplot as plt

    return plt


@method(
    Stage.VIZ,
    "manhattan",
    "Manhattan plot",
    description="Genome-wide -log10(p) scatter with chromosomes laid end to end "
    "and significance threshold lines.",
)
def manhattan_plot(df: pd.DataFrame, out_path, lines=None, title="GWAS Manhattan Plot") -> Path:
    plt = _lazy_plt()
    df = df.copy()
    df["neg_log10_p"] = -df["p"].apply(math.log10)
    df["_key"] = df["chr"].apply(chrom_sort_key)
    df = df.sort_values(["_key", "bp"]).reset_index(drop=True)

    chroms = sorted(df["chr"].unique(), key=chrom_sort_key)
    x, colors, ticks, tick_labels = [], [], [], []
    offset = 0.0
    for i, chrom in enumerate(chroms):
        sub = df[df["chr"] == chrom]
        xs = sub["bp"].to_numpy(dtype=float) - sub["bp"].min() + offset
        x.extend(xs)
        colors.extend([_PALETTE[i % 2]] * len(sub))
        if len(xs):
            ticks.append((xs.min() + xs.max()) / 2)
            offset = xs.max() + 1
        tick_labels.append(str(chrom))

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.scatter(x, df["neg_log10_p"], c=colors, s=10, alpha=0.8)
    for p_line, label in lines or []:
        if p_line and p_line > 0:
            ax.axhline(-math.log10(p_line), linestyle="--", linewidth=1, label=label)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title(title)
    if lines:
        ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    out_path = Path(out_path)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


@method(
    Stage.VIZ,
    "qq",
    "QQ plot",
    description="Quantile-quantile plot of observed vs. expected -log10(p) with "
    "the genomic inflation factor (λ) annotated.",
    aliases=("quantile",),
)
def qq_plot(df: pd.DataFrame, out_path, title="GWAS QQ Plot") -> Path:
    plt = _lazy_plt()
    p = np.sort(np.asarray(df["p"], dtype=float))
    m = p.size
    expected = -np.log10((np.arange(1, m + 1) - 0.5) / m)
    observed = -np.log10(p)
    lam = genomic_inflation(df["p"])

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=10, alpha=0.8, color="#2c7fb8")
    lim = max(expected.max(), observed.max())
    ax.plot([0, lim], [0, lim], color="red", linestyle="--", linewidth=1)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    ax.set_title(f"{title}  (λ = {lam:.3f})")
    fig.tight_layout()
    out_path = Path(out_path)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# --- Planned plots ---------------------------------------------------------

planned(
    Stage.VIZ,
    "regional",
    "Regional association plot",
    description="Zoomed locus plot with LD coloring and gene track (LocusZoom-style).",
)
planned(
    Stage.VIZ,
    "ld_heatmap",
    "LD heatmap",
    description="Pairwise linkage-disequilibrium heatmap for a region.",
)
planned(
    Stage.VIZ,
    "effects",
    "Effect-size / forest plot",
    description="Effect sizes and confidence intervals for top markers.",
)
