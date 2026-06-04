"""Tests for the significance / multiple-testing methods."""

from __future__ import annotations

import numpy as np
import pytest

from bioinsight import thresholds as T


@pytest.fixture
def pvals():
    # 1 very strong, 1 moderate, plus nulls.
    return np.array([1e-9, 1e-3, 0.2, 0.5, 0.9])


def test_genomewide_fixed_cutoff(pvals):
    r = T.genomewide(pvals, threshold=5e-8)
    assert r.n_significant == 1
    assert bool(r.significant[0]) is True
    assert r.line_p == 5e-8


def test_bonferroni_matches_definition(pvals):
    r = T.bonferroni(pvals, alpha=0.05)
    m = len(pvals)
    assert r.line_p == pytest.approx(0.05 / m)
    # adjusted p = p * m capped at 1
    assert r.adjusted_p[0] == pytest.approx(1e-9 * m)
    assert r.adjusted_p[-1] == 1.0
    assert r.n_significant == int((pvals <= 0.05 / m).sum())


def test_sidak_close_to_bonferroni(pvals):
    b = T.bonferroni(pvals, alpha=0.05)
    s = T.sidak(pvals, alpha=0.05)
    # Šidák cutoff is slightly larger (less conservative) than Bonferroni.
    assert s.line_p >= b.line_p
    assert s.n_significant >= b.n_significant


def test_holm_at_least_as_powerful_as_bonferroni(pvals):
    b = T.bonferroni(pvals, alpha=0.05)
    h = T.holm(pvals, alpha=0.05)
    assert h.n_significant >= b.n_significant


def test_fdr_bh_more_powerful_than_fwer(pvals):
    bh = T.fdr_bh(pvals, alpha=0.05)
    bonf = T.bonferroni(pvals, alpha=0.05)
    assert bh.n_significant >= bonf.n_significant
    # adjusted q-values are monotone non-decreasing in p order
    order = np.argsort(pvals)
    q = bh.adjusted_p[order]
    assert np.all(np.diff(q) >= -1e-12)


def test_fdr_by_more_conservative_than_bh(pvals):
    bh = T.fdr_bh(pvals, alpha=0.05)
    by = T.fdr_by(pvals, alpha=0.05)
    assert by.n_significant <= bh.n_significant


def test_masks_are_boolean_and_aligned(pvals):
    for fn in (T.genomewide, T.bonferroni, T.holm, T.fdr_bh, T.fdr_by, T.sidak):
        r = fn(pvals)
        assert r.significant.dtype == bool
        assert len(r.significant) == len(pvals)
        assert r.n_significant == int(r.significant.sum())
