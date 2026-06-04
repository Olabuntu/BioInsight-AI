"""End-to-end pipeline tests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from bioinsight import pipeline
from bioinsight.loaders import load_results

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "sample_gwas.csv"


def _toy():
    return pd.DataFrame(
        {
            "snp": ["a", "b", "c", "d"],
            "chr": [1, 1, 2, 2],
            "bp": [100, 200, 300, 400],
            "p": [1e-9, 0.5, 1e-6, 0.8],
        }
    )


def test_analyze_results_genomewide():
    s = pipeline.analyze_results(_toy(), method="genomewide")
    assert s.n_snps == 4
    assert s.primary.n_significant == 1
    assert s.top_hits.iloc[0]["snp"] == "a"
    assert "significant" in s.data.columns


def test_analyze_results_builds_method_comparison():
    s = pipeline.analyze_results(load_results(EXAMPLE), method="fdr_bh")
    # every available threshold method should appear in the comparison
    assert len(s.method_comparison) >= 5
    assert all(isinstance(v, int) for v in s.method_comparison.values())


def test_analyze_results_rejects_unknown_param():
    with pytest.raises(ValueError):
        pipeline.analyze_results(_toy(), method="fdr_bh", bogus=1)


def test_run_results_analysis_writes_outputs(tmp_path):
    out = pipeline.run_results_analysis(
        EXAMPLE, method="bonferroni", outdir=tmp_path, plots=("manhattan", "qq")
    )
    assert out["report_path"].exists()
    assert (tmp_path / "manhattan.png").exists()
    assert (tmp_path / "qq.png").exists()
    assert "GWAS tested" in out["interpretation"]


def test_run_with_planned_plot_raises(tmp_path):
    with pytest.raises(NotImplementedError):
        pipeline.run_results_analysis(
            EXAMPLE, outdir=tmp_path, plots=("regional",)
        )
