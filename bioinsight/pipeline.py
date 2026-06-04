"""Pipeline orchestration — compose chosen components into an analysis.

The full vision (genotypes -> QC -> structure -> model -> threshold -> viz ->
report) is described in DESIGN.md. What runs end-to-end today is the
results-analysis path: load a results table, apply a chosen significance
method, render plots, and produce a report. Other stages dispatch through the
same registry and raise a clear "planned" error until implemented.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import registry
from .registry import Stage, Status
from .stats import chrom_sort_key, genomic_inflation
from .thresholds import ThresholdResult


@dataclass
class GwasSummary:
    """Everything needed to report on a results-analysis run."""

    data: pd.DataFrame  # normalized + annotated (significant, p_adj, neg_log10_p)
    n_snps: int
    n_chromosomes: int
    chromosomes: list
    lambda_gc: float
    primary: ThresholdResult
    top_hits: pd.DataFrame
    method_comparison: dict = field(default_factory=dict)  # label -> n_significant


def apply_threshold(df: pd.DataFrame, method: str, **params) -> ThresholdResult:
    """Apply a single significance method (by key/alias) to a results frame.

    Supplied params are validated and cast against the method's declared spec;
    unknown params raise a clear error.
    """
    comp = registry.resolve(Stage.THRESHOLD, method)
    kwargs = comp.resolve_params(params)
    return comp.handler(df["p"].to_numpy(), **kwargs)


def analyze_results(
    df: pd.DataFrame,
    method: str = "genomewide",
    top_n: int = 10,
    compare_all: bool = True,
    **params,
) -> GwasSummary:
    """Run the results-analysis pipeline on a normalized results frame."""
    df = df.copy()
    df["neg_log10_p"] = -df["p"].apply(math.log10)

    primary = apply_threshold(df, method, **params)
    df["significant"] = primary.significant
    df["p_adj"] = primary.adjusted_p if primary.adjusted_p is not None else pd.NA

    comparison: dict = {}
    if compare_all:
        for comp in registry.list_components(Stage.THRESHOLD):
            if comp.status is not Status.AVAILABLE:
                continue
            # Lenient: only forward params this method actually accepts.
            safe = {k: v for k, v in params.items() if comp.get_param(k)}
            try:
                r = comp.handler(df["p"].to_numpy(), **comp.resolve_params(safe))
                comparison[r.label] = r.n_significant
            except Exception:
                continue

    chroms = sorted(df["chr"].unique(), key=chrom_sort_key)
    top = df.nsmallest(top_n, "p").reset_index(drop=True)

    return GwasSummary(
        data=df,
        n_snps=len(df),
        n_chromosomes=len(chroms),
        chromosomes=chroms,
        lambda_gc=genomic_inflation(df["p"]),
        primary=primary,
        top_hits=top,
        method_comparison=comparison,
    )


def run_results_analysis(
    input_path,
    method: str = "genomewide",
    outdir="bioinsight_report",
    plots=("manhattan", "qq"),
    top_n: int = 10,
    **params,
) -> dict:
    """End-to-end: load -> threshold -> plots -> interpret -> report.

    Returns a dict with the summary, interpretation text, and output paths.
    """
    from . import ai, report

    loader = registry.resolve(Stage.IO, "results")
    df = loader.handler(input_path)

    summary = analyze_results(df, method=method, top_n=top_n, **params)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_paths: dict = {}
    for plot_key in plots:
        comp = registry.resolve(Stage.VIZ, plot_key)
        out = outdir / f"{plot_key}.png"
        if plot_key == "manhattan":
            lines = [(summary.primary.line_p, summary.primary.label)]
            comp.handler(df, out, lines=lines)
        else:
            comp.handler(df, out)
        plot_paths[comp.name] = out

    interpretation = ai.interpret(summary)
    md = report.render_markdown(summary, interpretation, plot_paths)
    report_path = outdir / "report.md"
    report_path.write_text(md)

    return {
        "summary": summary,
        "interpretation": interpretation,
        "report_path": report_path,
        "plot_paths": plot_paths,
    }
