"""Tests for the component registry and typed parameter system."""

from __future__ import annotations

import pytest

import bioinsight  # noqa: F401  (populates the registry on import)
from bioinsight import registry
from bioinsight.registry import Param, Stage, Status


def test_every_stage_has_components():
    for stage in Stage:
        assert registry.list_components(stage), f"no components for {stage}"


def test_available_threshold_methods_have_handlers():
    for c in registry.list_components(Stage.THRESHOLD):
        if c.status is Status.AVAILABLE:
            assert callable(c.handler)


def test_get_by_alias():
    assert registry.get(Stage.THRESHOLD, "bh").key == "fdr_bh"
    assert registry.get(Stage.THRESHOLD, "BH").key == "fdr_bh"  # case-insensitive


def test_resolve_planned_raises_notimplemented():
    with pytest.raises(NotImplementedError):
        registry.resolve(Stage.MODEL, "mlm")


def test_param_cast_types():
    assert Param("a", "int").cast("5") == 5
    assert Param("a", "float").cast("1e-6") == 1e-6
    assert Param("a", "bool").cast("yes") is True
    assert Param("a", "bool").cast("0") is False


def test_param_cast_bad_value_raises():
    with pytest.raises(ValueError):
        Param("alpha", "float").cast("abc")


def test_param_choices_enforced():
    p = Param("kinship", "str", choices=("vanraden", "ibs"))
    assert p.cast("ibs") == "ibs"
    with pytest.raises(ValueError):
        p.cast("nonsense")


def test_resolve_params_fills_defaults_and_rejects_unknown():
    comp = registry.get(Stage.THRESHOLD, "bonferroni")
    assert comp.resolve_params({}) == {"alpha": 0.05}
    assert comp.resolve_params({"alpha": "0.1"}) == {"alpha": 0.1}
    with pytest.raises(ValueError):
        comp.resolve_params({"unknown": 1})


def test_resolve_params_required_missing_raises():
    comp = registry.get(Stage.STRUCTURE, "qmatrix")  # k is required
    with pytest.raises(ValueError):
        comp.resolve_params({})
    assert comp.resolve_params({"k": "3"}) == {"k": 3}
