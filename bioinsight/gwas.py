"""GWAS analysis: load results, find significant SNPs, plot, summarize.

A GWAS results file is expected to be a CSV (or TSV) with at least these
columns (names are matched case-insensitively, common aliases accepted):

    snp / rsid / marker   - variant identifier
    chr / chrom           - chromosome
    bp / pos / position   - base-pair position
    p / pval / p_value    - association p-value

Anything extra is preserved and ignored.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Genome-wide significance threshold commonly used in human GWAS.
GENOME_WIDE_THRESHOLD = 5e-8
# A looser "suggestive" threshold.
SUGGESTIVE_THRESHOLD = 1e-5

# Map of canonical column name -> accepted aliases (lowercase).
_COLUMN_ALIASES = {
    "snp": ["snp", "rsid", "rs", "marker", "variant", "id"],
    "chr": ["chr", "chrom", "chromosome", "#chrom"],
    "bp": ["bp", "pos", "position", "base_pair", "ps"],
    "p": ["p", "pval", "p_value", "pvalue", "p.value", "p_bonferroni"],
}


@dataclass
class GwasResult:
    """Outcome of a GWAS analysis run."""

    data: pd.DataFrame
    n_snps: int
    n_significant: int
    n_suggestive: int
    threshold: float
    top_hits: pd.DataFrame
    lambda_gc: float
    chromosomes: list = field(default_factory=list)


def _resolve_columns(df: pd.DataFrame) -> dict:
    """Return a mapping {canonical: actual_column_name} found in df."""
    lowered = {c.lower(): c for c in df.columns}
    resolved = {}
    for canonical, aliases in _COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in lowered:
                resolved[canonical] = lowered[alias]
                break
    return resolved


def load_gwas(path: str | Path) -> pd.DataFrame:
    """Load a GWAS results file (CSV or TSV) into a normalized DataFrame.

    The returned frame is guaranteed to have columns: snp, chr, bp, p.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GWAS file not found: {path}")

    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    df = pd.read_csv(path, sep=sep)

    resolved = _resolve_columns(df)
    missing = [c for c in ("p",) if c not in resolved]
    if missing:
        raise ValueError(
            "Could not find a p-value column. Looked for one of: "
            f"{_COLUMN_ALIASES['p']}. Found columns: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["p"] = pd.to_numeric(df[resolved["p"]], errors="coerce")
    out["snp"] = (
        df[resolved["snp"]].astype(str)
        if "snp" in resolved
        else [f"snp{i}" for i in range(len(df))]
    )
    out["chr"] = df[resolved["chr"]].astype(str) if "chr" in resolved else "NA"
    out["bp"] = (
        pd.to_numeric(df[resolved["bp"]], errors="coerce")
        if "bp" in resolved
        else range(len(df))
    )

    out = out.dropna(subset=["p"])
    out = out[(out["p"] > 0) & (out["p"] <= 1)]
    if out.empty:
        raise ValueError("No valid p-values (0 < p <= 1) found in the file.")
    return out.reset_index(drop=True)


def _genomic_inflation(pvals: pd.Series) -> float:
    """Compute the genomic inflation factor lambda_GC.

    lambda_GC = median(chi-square stats) / 0.456 (the median of a 1-df chi2).
    Values near 1.0 indicate well-calibrated statistics.
    """
    try:
        from scipy.stats import chi2

        chisq = chi2.isf(pvals, df=1)
        return float(chisq.median() / chi2.ppf(0.5, df=1))
    except Exception:
        # Fallback without scipy: approximate via normal quantiles.
        # chi2_1 stat = z^2 where z is the standard-normal two-sided quantile.
        from statistics import NormalDist

        nd = NormalDist()
        chisq = [(nd.inv_cdf(1 - p / 2)) ** 2 for p in pvals]
        chisq.sort()
        n = len(chisq)
        median = (
            chisq[n // 2]
            if n % 2
            else (chisq[n // 2 - 1] + chisq[n // 2]) / 2
        )
        return float(median / 0.4549)


def analyze(
    df: pd.DataFrame,
    threshold: float = GENOME_WIDE_THRESHOLD,
    top_n: int = 10,
) -> GwasResult:
    """Run the core GWAS summary analysis on a normalized DataFrame."""
    df = df.copy()
    df["neg_log10_p"] = -df["p"].apply(math.log10)

    n_sig = int((df["p"] <= threshold).sum())
    n_sugg = int((df["p"] <= SUGGESTIVE_THRESHOLD).sum())

    top = df.nsmallest(top_n, "p").reset_index(drop=True)
    chroms = sorted(df["chr"].unique(), key=_chrom_sort_key)

    return GwasResult(
        data=df,
        n_snps=len(df),
        n_significant=n_sig,
        n_suggestive=n_sugg,
        threshold=threshold,
        top_hits=top,
        lambda_gc=_genomic_inflation(df["p"]),
        chromosomes=chroms,
    )


def _chrom_sort_key(c: str):
    """Sort chromosomes numerically when possible, then alphabetically."""
    s = str(c).replace("chr", "").replace("Chr", "")
    try:
        return (0, int(s))
    except ValueError:
        return (1, s)


def manhattan_plot(result: GwasResult, out_path: str | Path) -> Path:
    """Render a Manhattan plot to ``out_path`` and return the path."""
    import matplotlib

    matplotlib.use("Agg")  # headless / no display required
    import matplotlib.pyplot as plt

    df = result.data.copy()
    df["_chrkey"] = df["chr"].apply(lambda c: _chrom_sort_key(c))
    df = df.sort_values(["_chrkey", "bp"]).reset_index(drop=True)

    # Build a cumulative x position so chromosomes lay end to end.
    x = []
    ticks = []
    tick_labels = []
    offset = 0
    colors = []
    palette = ["#2c7fb8", "#7fcdbb"]
    for i, chrom in enumerate(result.chromosomes):
        sub = df[df["chr"] == chrom]
        xs = sub["bp"].to_numpy() - sub["bp"].min() + offset
        x.extend(xs)
        colors.extend([palette[i % 2]] * len(sub))
        ticks.append(offset + (xs.max() - xs.min()) / 2 if len(xs) else offset)
        tick_labels.append(str(chrom))
        offset = (xs.max() if len(xs) else offset) + 1

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.scatter(x, df["neg_log10_p"], c=colors, s=10, alpha=0.8)
    ax.axhline(
        -math.log10(result.threshold),
        color="red",
        linestyle="--",
        linewidth=1,
        label=f"genome-wide ({result.threshold:g})",
    )
    ax.axhline(
        -math.log10(SUGGESTIVE_THRESHOLD),
        color="blue",
        linestyle=":",
        linewidth=1,
        label=f"suggestive ({SUGGESTIVE_THRESHOLD:g})",
    )
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, rotation=0, fontsize=8)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title("GWAS Manhattan Plot")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()

    out_path = Path(out_path)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
