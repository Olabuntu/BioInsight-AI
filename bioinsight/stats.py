"""Small statistical helpers shared across the pipeline."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def neg_log10(p) -> pd.Series:
    """Return -log10(p) for a series/array of p-values."""
    arr = np.asarray(p, dtype=float)
    return pd.Series(-np.log10(arr))


def genomic_inflation(pvals) -> float:
    """Genomic inflation factor lambda_GC.

    lambda_GC = median(observed chi-square) / median(expected chi-square),
    with 1 degree of freedom. Values near 1.0 indicate well-calibrated
    statistics; values >> 1 suggest residual confounding (e.g. population
    stratification).

    Reference: Devlin & Roeder 1999, Biometrics 55(4):997-1004 (genomic control).
    """
    p = np.asarray(pvals, dtype=float)
    p = p[(p > 0) & (p <= 1)]
    if p.size == 0:
        return float("nan")
    try:
        from scipy.stats import chi2

        chisq = chi2.isf(p, df=1)
        return float(np.median(chisq) / chi2.ppf(0.5, df=1))
    except Exception:
        from statistics import NormalDist

        nd = NormalDist()
        chisq = np.array([(nd.inv_cdf(1 - x / 2)) ** 2 for x in p])
        return float(np.median(chisq) / 0.4549364231195736)


def chrom_sort_key(c) -> tuple:
    """Sort chromosomes numerically when possible, else alphabetically."""
    s = str(c).lower().replace("chr", "")
    try:
        return (0, int(s), "")
    except ValueError:
        # Handle X/Y/MT and arbitrary scaffold names.
        special = {"x": 23, "y": 24, "xy": 25, "mt": 26, "m": 26}
        if s in special:
            return (0, special[s], "")
        return (1, math.inf, s)
