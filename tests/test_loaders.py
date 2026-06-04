"""Tests for the data-import stage."""

from __future__ import annotations

from pathlib import Path

import pytest

from bioinsight.loaders import load_results

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "sample_gwas.csv"


def test_load_example_normalizes_columns():
    df = load_results(EXAMPLE)
    assert {"snp", "chr", "bp", "p"}.issubset(df.columns)
    assert len(df) > 0
    assert df["p"].between(0, 1).all()


def test_load_resolves_aliases(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text("rsID,CHROM,POS,P_VALUE\nrs1,1,500,1e-8\nrs2,2,900,0.4\n")
    df = load_results(f)
    assert df.iloc[0]["snp"] == "rs1"
    assert df.iloc[0]["chr"] == "1"


def test_load_tsv(tmp_path):
    f = tmp_path / "g.tsv"
    f.write_text("snp\tchr\tbp\tp\nrs1\t1\t500\t0.01\n")
    df = load_results(f)
    assert len(df) == 1


def test_load_rejects_missing_pvalue(tmp_path):
    f = tmp_path / "bad.csv"
    f.write_text("snp,chr,bp\nrs1,1,500\n")
    with pytest.raises(ValueError):
        load_results(f)


def test_load_filters_invalid_pvalues(tmp_path):
    f = tmp_path / "g.csv"
    f.write_text("snp,p\nrs1,0.01\nrs2,1.5\nrs3,-0.1\nrs4,abc\n")
    df = load_results(f)
    assert len(df) == 1
    assert df.iloc[0]["snp"] == "rs1"


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_results(tmp_path / "nope.csv")
