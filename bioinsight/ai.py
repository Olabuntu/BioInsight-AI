"""AI-assisted interpretation of analysis results.

v0.1.0 ships a deterministic, rule-based interpreter so the tool works with
zero configuration and no API keys. The function signature is intentionally
LLM-shaped (``interpret_gwas`` takes a result and returns prose) so a real
language-model backend can be slotted in behind the same interface later.
"""

from __future__ import annotations

from .gwas import GENOME_WIDE_THRESHOLD, GwasResult


def interpret_gwas(result: GwasResult) -> str:
    """Return a plain-English interpretation of a GWAS result."""
    lines: list[str] = []

    lines.append(
        f"This GWAS tested {result.n_snps:,} variants across "
        f"{len(result.chromosomes)} chromosome(s)."
    )

    if result.n_significant > 0:
        lead = result.top_hits.iloc[0]
        lines.append(
            f"{result.n_significant} variant(s) reached genome-wide "
            f"significance (p ≤ {result.threshold:g}). The strongest "
            f"association is {lead['snp']} on chromosome {lead['chr']} "
            f"(p = {lead['p']:.2e}), a strong candidate for follow-up."
        )
    elif result.n_suggestive > 0:
        lines.append(
            f"No variant passed the genome-wide threshold, but "
            f"{result.n_suggestive} reached suggestive significance "
            f"(p ≤ 1e-5). These may be worth replicating in a larger "
            f"cohort before drawing conclusions."
        )
    else:
        lines.append(
            "No variant reached suggestive or genome-wide significance. "
            "This often indicates limited statistical power (small sample "
            "size) rather than the absence of any true effect."
        )

    # Genomic inflation interpretation.
    lam = result.lambda_gc
    if lam < 1.05:
        lines.append(
            f"The genomic inflation factor (λ = {lam:.3f}) is close to "
            f"1.0, suggesting well-calibrated statistics with little evidence "
            f"of population stratification or systematic bias."
        )
    elif lam < 1.2:
        lines.append(
            f"The genomic inflation factor (λ = {lam:.3f}) is mildly "
            f"elevated. Some inflation is present; consider adjusting for "
            f"principal components or relatedness."
        )
    else:
        lines.append(
            f"The genomic inflation factor (λ = {lam:.3f}) is high, "
            f"which points to substantial confounding (e.g. population "
            f"stratification). Interpret hits with caution and re-run with "
            f"appropriate covariates."
        )

    return "\n".join(lines)


def is_llm_available() -> bool:
    """Placeholder hook for a future LLM backend. Always False in v0.1.0."""
    return False
