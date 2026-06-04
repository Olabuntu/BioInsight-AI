"""Validation tests: cross-check our methods against reference implementations.

These guard the statistical correctness of the correction methods by comparing
to statsmodels (the canonical Python implementation) and confirm the genomic
inflation factor behaves as theory predicts.
"""

from __future__ import annotations

import numpy as np
import pytest

from bioinsight import thresholds as T
from bioinsight.stats import genomic_inflation

statsmodels = pytest.importorskip("statsmodels")
from statsmodels.stats.multitest import multipletests  # noqa: E402


@pytest.fixture
def pvals():
    rng = np.random.default_rng(0)
    return np.concatenate([rng.uniform(0, 1, 95), [1e-9, 1e-7, 1e-4, 1e-3, 5e-3]])


@pytest.mark.parametrize(
    "fn,sm_method",
    [
        (T.bonferroni, "bonferroni"),
        (T.sidak, "sidak"),
        (T.holm, "holm"),
        (T.fdr_bh, "fdr_bh"),
        (T.fdr_by, "fdr_by"),
    ],
)
def test_matches_statsmodels(pvals, fn, sm_method):
    alpha = 0.05
    r = fn(pvals, alpha=alpha)
    rej_sm, padj_sm, _, _ = multipletests(pvals, alpha=alpha, method=sm_method)
    assert np.array_equal(r.significant, rej_sm)
    assert np.allclose(r.adjusted_p, padj_sm, atol=1e-12)
    assert r.n_significant == int(rej_sm.sum())


def test_lambda_gc_null_is_one():
    rng = np.random.default_rng(1)
    p = rng.uniform(0, 1, 100000)
    assert genomic_inflation(p) == pytest.approx(1.0, abs=0.03)


def test_lambda_gc_detects_inflation():
    from scipy.stats import chi2

    rng = np.random.default_rng(2)
    z = rng.normal(0, 1, 100000)
    p = chi2.sf((z**2) * 1.5, df=1)  # stats inflated by 1.5x
    assert genomic_inflation(p) == pytest.approx(1.5, abs=0.05)
