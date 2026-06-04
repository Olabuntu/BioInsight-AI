"""Significance / multiple-testing correction stage.

These methods operate directly on a vector of association p-values, so they
work today on any results table. Pick one (or compare several) to decide which
markers are significant.

Each handler takes a 1-D array-like of p-values plus keyword params and returns
a :class:`ThresholdResult`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .registry import Param, Stage, method


@dataclass
class ThresholdResult:
    """Outcome of applying a significance method to a set of p-values."""

    method: str
    n_tests: int
    n_significant: int
    significant: np.ndarray  # bool mask aligned to the input order
    adjusted_p: Optional[np.ndarray]  # adjusted p-values, or None if undefined
    line_p: Optional[float]  # representative raw-p cutoff for plotting
    label: str  # human-readable legend label
    alpha: Optional[float] = None


def _clean(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    return p


@method(
    Stage.THRESHOLD,
    "genomewide",
    "Genome-wide fixed threshold",
    description="Declare significance at a fixed p-value cutoff (default 5e-8, "
    "the standard human GWAS threshold).",
    aliases=("gw", "fixed", "5e-8"),
    params=(
        Param("threshold", "float", 5e-8, help="Absolute p-value cutoff."),
    ),
    reference="Dudbridge & Gusnanto 2008, Genet Epidemiol 32:227-234; "
    "Pe'er et al. 2008, Genet Epidemiol 32:381-385",
)
def genomewide(pvals, threshold: float = 5e-8, **_) -> ThresholdResult:
    p = _clean(pvals)
    mask = p <= threshold
    return ThresholdResult(
        method="genomewide",
        n_tests=p.size,
        n_significant=int(mask.sum()),
        significant=mask,
        adjusted_p=None,
        line_p=threshold,
        label=f"genome-wide (p ≤ {threshold:g})",
    )


@method(
    Stage.THRESHOLD,
    "suggestive",
    "Suggestive fixed threshold",
    description="A looser fixed cutoff (default 1e-5) for hypothesis generation.",
    params=(
        Param("threshold", "float", 1e-5, help="Absolute p-value cutoff."),
    ),
)
def suggestive(pvals, threshold: float = 1e-5, **_) -> ThresholdResult:
    p = _clean(pvals)
    mask = p <= threshold
    return ThresholdResult(
        method="suggestive",
        n_tests=p.size,
        n_significant=int(mask.sum()),
        significant=mask,
        adjusted_p=None,
        line_p=threshold,
        label=f"suggestive (p ≤ {threshold:g})",
    )


@method(
    Stage.THRESHOLD,
    "bonferroni",
    "Bonferroni correction",
    description="Family-wise error control: significant if p ≤ alpha / n_tests.",
    aliases=("bonf",),
    params=(
        Param("alpha", "float", 0.05, help="Family-wise error rate."),
    ),
    reference="Bonferroni 1936; Dunn 1961, J Am Stat Assoc 56:52-64",
)
def bonferroni(pvals, alpha: float = 0.05, **_) -> ThresholdResult:
    p = _clean(pvals)
    m = p.size
    cutoff = alpha / m
    adj = np.minimum(p * m, 1.0)
    mask = p <= cutoff
    return ThresholdResult(
        method="bonferroni",
        n_tests=m,
        n_significant=int(mask.sum()),
        significant=mask,
        adjusted_p=adj,
        line_p=cutoff,
        label=f"Bonferroni (α={alpha:g}, p ≤ {cutoff:.2e})",
        alpha=alpha,
    )


@method(
    Stage.THRESHOLD,
    "sidak",
    "Šidák correction",
    description="Family-wise error control assuming independence: "
    "significant if p ≤ 1 - (1 - alpha)^(1/n_tests).",
    params=(
        Param("alpha", "float", 0.05, help="Family-wise error rate."),
    ),
    reference="Šidák 1967, J Am Stat Assoc 62:626-633",
)
def sidak(pvals, alpha: float = 0.05, **_) -> ThresholdResult:
    p = _clean(pvals)
    m = p.size
    cutoff = 1.0 - (1.0 - alpha) ** (1.0 / m)
    adj = np.minimum(1.0 - np.power(1.0 - p, m), 1.0)
    mask = p <= cutoff
    return ThresholdResult(
        method="sidak",
        n_tests=m,
        n_significant=int(mask.sum()),
        significant=mask,
        adjusted_p=adj,
        line_p=cutoff,
        label=f"Šidák (α={alpha:g}, p ≤ {cutoff:.2e})",
        alpha=alpha,
    )


@method(
    Stage.THRESHOLD,
    "holm",
    "Holm step-down",
    description="Sequentially-rejective family-wise control; uniformly more "
    "powerful than Bonferroni.",
    aliases=("holm-bonferroni",),
    params=(
        Param("alpha", "float", 0.05, help="Family-wise error rate."),
    ),
    reference="Holm 1979, Scand J Stat 6:65-70",
)
def holm(pvals, alpha: float = 0.05, **_) -> ThresholdResult:
    p = _clean(pvals)
    m = p.size
    order = np.argsort(p)
    sorted_p = p[order]
    # Step-down thresholds: alpha / (m - i) for i = 0..m-1.
    thresholds = alpha / (m - np.arange(m))
    passed = sorted_p <= thresholds
    # Reject up to (but not including) the first failure.
    if passed.all():
        n_rej = m
    else:
        n_rej = int(np.argmin(passed))  # index of first False
    reject_sorted = np.zeros(m, dtype=bool)
    reject_sorted[:n_rej] = True
    # Holm-adjusted p-values: cumulative max of (m - i) * p_(i), capped at 1.
    raw = (m - np.arange(m)) * sorted_p
    adj_sorted = np.minimum(np.maximum.accumulate(raw), 1.0)
    mask = np.zeros(m, dtype=bool)
    adj = np.empty(m)
    mask[order] = reject_sorted
    adj[order] = adj_sorted
    line_p = float(sorted_p[n_rej - 1]) if n_rej > 0 else None
    return ThresholdResult(
        method="holm",
        n_tests=m,
        n_significant=int(mask.sum()),
        significant=mask,
        adjusted_p=adj,
        line_p=line_p,
        label=f"Holm (α={alpha:g})",
        alpha=alpha,
    )


def _step_up(pvals, alpha, c_m, name, label_prefix):
    """Shared Benjamini-Hochberg / Benjamini-Yekutieli step-up procedure.

    c_m is the BH/BY constant (1 for BH; sum(1/i) for BY).
    """
    p = _clean(pvals)
    m = p.size
    order = np.argsort(p)
    sorted_p = p[order]
    ranks = np.arange(1, m + 1)
    crit = ranks / (m * c_m) * alpha
    below = sorted_p <= crit
    if below.any():
        k = np.max(np.flatnonzero(below))  # largest passing rank (0-indexed)
        n_rej = k + 1
        line_p = float(sorted_p[k])
    else:
        n_rej = 0
        line_p = None
    reject_sorted = np.zeros(m, dtype=bool)
    reject_sorted[:n_rej] = True
    # Adjusted (q) values via cumulative minimum from the largest p.
    adj_sorted = np.minimum.accumulate((m * c_m / ranks * sorted_p)[::-1])[::-1]
    adj_sorted = np.minimum(adj_sorted, 1.0)
    mask = np.zeros(m, dtype=bool)
    adj = np.empty(m)
    mask[order] = reject_sorted
    adj[order] = adj_sorted
    return ThresholdResult(
        method=name,
        n_tests=m,
        n_significant=int(mask.sum()),
        significant=mask,
        adjusted_p=adj,
        line_p=line_p,
        label=f"{label_prefix} (FDR={alpha:g})",
        alpha=alpha,
    )


@method(
    Stage.THRESHOLD,
    "fdr_bh",
    "Benjamini-Hochberg FDR",
    description="Controls the false discovery rate (assumes independence or "
    "positive dependence). More powerful than family-wise methods.",
    aliases=("fdr", "bh"),
    params=(
        Param("alpha", "float", 0.05, help="Target false discovery rate."),
    ),
    reference="Benjamini & Hochberg 1995, J R Stat Soc B 57:289-300",
)
def fdr_bh(pvals, alpha: float = 0.05, **_) -> ThresholdResult:
    return _step_up(pvals, alpha, 1.0, "fdr_bh", "Benjamini-Hochberg")


@method(
    Stage.THRESHOLD,
    "fdr_by",
    "Benjamini-Yekutieli FDR",
    description="FDR control valid under arbitrary dependence (more "
    "conservative than Benjamini-Hochberg).",
    aliases=("by",),
    params=(
        Param("alpha", "float", 0.05, help="Target false discovery rate."),
    ),
    reference="Benjamini & Yekutieli 2001, Ann Stat 29:1165-1188",
)
def fdr_by(pvals, alpha: float = 0.05, **_) -> ThresholdResult:
    m = np.asarray(pvals).size
    c_m = float(np.sum(1.0 / np.arange(1, m + 1)))
    return _step_up(pvals, alpha, c_m, "fdr_by", "Benjamini-Yekutieli")
