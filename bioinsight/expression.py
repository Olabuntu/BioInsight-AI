"""Expression analysis (planned for a future release).

v0.1.0 focuses on GWAS. This module defines the intended interface so the
roadmap is visible in code, but the implementation is not yet available.
"""

from __future__ import annotations

from pathlib import Path


def analyze_expression(path: str | Path) -> None:
    """Analyze a gene-expression matrix (differential expression, etc.).

    Not yet implemented in v0.1.0.
    """
    raise NotImplementedError(
        "Expression analysis is planned for a future release. "
        "v0.1.0 supports GWAS analysis — try `bioinsight gwas <file>`."
    )
