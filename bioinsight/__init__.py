"""BioInsight: a modular GWAS pipeline — pick the tools you need.

Importing the package registers every available and planned component into the
shared registry (see :mod:`bioinsight.registry`).
"""

from __future__ import annotations

__version__ = "0.2.0"

# Import side-effect: each module registers its components on import.
from . import loaders, qc, structure, models, thresholds, viz, ai  # noqa: E402,F401
from . import registry  # noqa: E402

__all__ = ["registry", "__version__"]
