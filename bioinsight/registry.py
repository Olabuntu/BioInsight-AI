"""Component registry — the backbone of the modular pipeline.

Every capability in BioInsight (a file loader, a QC filter, an association
model, a significance method, a plot, a report) is a :class:`Component`
registered against a :class:`Stage`. Users compose a pipeline by *picking*
which component to use at each stage.

Components advertise a :class:`Status`: ``AVAILABLE`` ones run today;
``PLANNED`` ones are part of the published design and raise a helpful error
until implemented. This lets the architecture be complete and inspectable
(`bioinsight list`) long before every method is coded.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Optional


class Stage(str, Enum):
    """Ordered stages of a GWAS pipeline."""

    IO = "input"
    QC = "qc"
    STRUCTURE = "structure"
    MODEL = "model"
    THRESHOLD = "threshold"
    VIZ = "viz"
    REPORT = "report"

    @property
    def title(self) -> str:
        return {
            "input": "Data import",
            "qc": "Quality control",
            "structure": "Population structure",
            "model": "Association model",
            "threshold": "Significance / correction",
            "viz": "Visualization",
            "report": "Reporting",
        }[self.value]


class Status(str, Enum):
    AVAILABLE = "available"
    PLANNED = "planned"


_CASTERS = {
    "str": str,
    "int": int,
    "float": float,
    "path": lambda v: Path(v),
    "bool": lambda v: str(v).strip().lower() in {"1", "true", "yes", "y", "on"},
}


@dataclass(frozen=True)
class Param:
    """A typed, self-describing input for a component.

    The CLI and pipeline use this to parse, validate, and document the
    method-specific inputs a user must (or may) provide.
    """

    name: str
    type: str = "float"  # one of: str | int | float | path | bool
    default: Any = None
    required: bool = False
    help: str = ""
    choices: tuple = ()

    def cast(self, raw: Any) -> Any:
        """Cast a raw (string) value to this param's type and validate choices."""
        if raw is None:
            return None
        try:
            value = _CASTERS[self.type](raw)
        except KeyError:
            raise ValueError(f"Unknown param type '{self.type}' for '{self.name}'.")
        except (ValueError, TypeError):
            raise ValueError(
                f"Parameter '{self.name}' expects a {self.type}; got {raw!r}."
            )
        if self.choices and value not in self.choices:
            raise ValueError(
                f"Parameter '{self.name}' must be one of {self.choices}; got {value!r}."
            )
        return value


@dataclass(frozen=True)
class Component:
    """A single pluggable capability."""

    key: str
    name: str
    stage: Stage
    status: Status
    description: str
    handler: Optional[Callable] = None
    aliases: tuple = ()
    params: tuple = ()  # tuple[Param, ...] — the inputs this component accepts
    reference: str = ""  # citation / method origin, when relevant

    @property
    def available(self) -> bool:
        return self.status is Status.AVAILABLE and self.handler is not None

    def get_param(self, name: str) -> Optional[Param]:
        for p in self.params:
            if p.name == name:
                return p
        return None

    def resolve_params(self, supplied: Optional[dict] = None) -> dict:
        """Merge user-supplied values with declared defaults.

        Casts each supplied value to its declared type, rejects unknown keys,
        validates ``choices`` and ``required``, and fills in defaults.
        """
        supplied = supplied or {}
        known = {p.name for p in self.params}
        unknown = set(supplied) - known
        if unknown:
            raise ValueError(
                f"Unknown parameter(s) for '{self.key}': {sorted(unknown)}. "
                f"Accepted: {sorted(known) or '(none)'}."
            )
        out: dict = {}
        for p in self.params:
            if p.name in supplied and supplied[p.name] is not None:
                out[p.name] = p.cast(supplied[p.name])
            elif p.default is not None:
                out[p.name] = p.default
            elif p.required:
                raise ValueError(
                    f"Missing required parameter '{p.name}' for '{self.key}'."
                )
        return out


# stage -> {key: Component}
_REGISTRY: dict[Stage, dict[str, Component]] = {s: {} for s in Stage}
# stage -> {alias: canonical key}
_ALIASES: dict[Stage, dict[str, str]] = {s: {} for s in Stage}


def register(component: Component) -> Component:
    """Register a component. Raises if the key is already taken in its stage."""
    bucket = _REGISTRY[component.stage]
    if component.key in bucket:
        raise ValueError(
            f"Component '{component.key}' already registered for stage "
            f"'{component.stage.value}'."
        )
    bucket[component.key] = component
    for alias in component.aliases:
        _ALIASES[component.stage][alias.lower()] = component.key
    return component


def method(
    stage: Stage,
    key: str,
    name: str,
    *,
    status: Status = Status.AVAILABLE,
    description: str = "",
    aliases: tuple = (),
    params: tuple = (),
    reference: str = "",
) -> Callable:
    """Decorator registering the wrapped function as an AVAILABLE component."""

    def decorator(func: Callable) -> Callable:
        register(
            Component(
                key=key,
                name=name,
                stage=stage,
                status=status,
                description=description,
                handler=func,
                aliases=aliases,
                params=tuple(params),
                reference=reference,
            )
        )
        return func

    return decorator


def planned(
    stage: Stage,
    key: str,
    name: str,
    *,
    description: str = "",
    aliases: tuple = (),
    params: tuple = (),
    reference: str = "",
) -> Component:
    """Register a PLANNED component (no handler yet)."""
    return register(
        Component(
            key=key,
            name=name,
            stage=stage,
            status=Status.PLANNED,
            description=description,
            handler=None,
            aliases=aliases,
            params=tuple(params),
            reference=reference,
        )
    )


def get(stage: Stage, key: str) -> Component:
    """Look up a component by key or alias within a stage."""
    bucket = _REGISTRY[stage]
    k = key.lower()
    if k in bucket:
        return bucket[k]
    if k in _ALIASES[stage]:
        return bucket[_ALIASES[stage][k]]
    available_keys = ", ".join(sorted(bucket)) or "(none)"
    raise KeyError(
        f"No '{key}' component for stage '{stage.value}'. "
        f"Available: {available_keys}."
    )


def resolve(stage: Stage, key: str) -> Component:
    """Like :func:`get`, but raise NotImplementedError for PLANNED components."""
    comp = get(stage, key)
    if not comp.available:
        raise NotImplementedError(
            f"'{comp.name}' ({stage.value}:{comp.key}) is on the roadmap but "
            f"not implemented in this version. See DESIGN.md for the plan. "
            f"Run `bioinsight list --stage {stage.value}` to see what works today."
        )
    return comp


def list_components(stage: Optional[Stage] = None) -> list[Component]:
    """Return all components, optionally filtered to one stage."""
    stages = [stage] if stage else list(Stage)
    out: list[Component] = []
    for s in stages:
        out.extend(_REGISTRY[s].values())
    return out
