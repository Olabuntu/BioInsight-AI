"""Tests for the GWAS analysis core."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from bioinsight import ai, gwas

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "sample_gwas.csv"


def _toy_df():
    return pd.DataFrame(
        {
            "snp": ["a", "b", "c", "d"],
            "chr": [1, 1, 2, 2],
            "bp": [100, 200, 300, 400],
            "p": [1e-9, 0.5, 1e-6, 0.8],
        }
    )


def test_load_example_normalizes_columns():
    df = gwas.load_gwas(EXAMPLE)
    assert set(["snp", "chr", "bp", "p"]).issubset(df.columns)
    assert len(df) > 0
    assert df["p"].between(0, 1).all()


def test_load_resolves_aliases(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text("rsID,CHROM,POS,P_VALUE\nrs1,1,500,1e-8\nrs2,2,900,0.4\n")
    df = gwas.load_gwas(f)
    assert list(df.columns) == ["p", "snp", "chr", "bp"]
    assert df.iloc[0]["snp"] == "rs1"


def test_load_rejects_missing_pvalue(tmp_path):
    f = tmp_path / "bad.csv"
    f.write_text("snp,chr,bp\nrs1,1,500\n")
    with pytest.raises(ValueError):
        gwas.load_gwas(f)


def test_load_filters_invalid_pvalues(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text("snp,p\nrs1,0.01\nrs2,1.5\nrs3,-0.1\nrs4,abc\n")
    df = gwas.load_gwas(f)
    assert len(df) == 1
    assert df.iloc[0]["snp"] == "rs1"


def test_analyze_counts_significance():
    result = gwas.analyze(_toy_df(), threshold=5e-8)
    assert result.n_snps == 4
    assert result.n_significant == 1  # only p=1e-9
    assert result.n_suggestive == 2  # 1e-9 and 1e-6
    assert result.top_hits.iloc[0]["snp"] == "a"


def test_analyze_on_example_finds_hits():
    df = gwas.load_gwas(EXAMPLE)
    result = gwas.analyze(df)
    assert result.n_significant >= 2
    assert result.n_suggestive >= result.n_significant


def test_manhattan_plot_writes_file(tmp_path):
    result = gwas.analyze(gwas.load_gwas(EXAMPLE))
    out = gwas.manhattan_plot(result, tmp_path / "m.png")
    assert out.exists()
    assert out.stat().st_size > 0


def test_interpret_mentions_lead_snp():
    result = gwas.analyze(_toy_df(), threshold=5e-8)
    text = ai.interpret_gwas(result)
    assert "a" in text
    assert "significance" in text.lower()


def test_interpret_handles_no_hits():
    df = pd.DataFrame({"snp": ["x", "y"], "chr": [1, 1], "bp": [1, 2], "p": [0.4, 0.9]})
    text = ai.interpret_gwas(gwas.analyze(df))
    assert "No variant" in text
