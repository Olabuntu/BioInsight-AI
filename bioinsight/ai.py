"""Interpretation stage — turn numbers into plain-English insight.

v0 ships a deterministic, rule-based interpreter so the tool works with zero
configuration and no API keys. The function is intentionally LLM-shaped
(``interpret`` takes a summary and returns prose) so a real language-model
backend can be slotted in behind the same interface later.
"""

from __future__ import annotations

from .registry import Stage, method, planned


@method(
    Stage.REPORT,
    "interpret",
    "AI-assisted interpretation",
    description="Plain-English summary of the result (rule-based; no API key "
    "required). Pluggable with an LLM backend in future.",
    aliases=("ai", "insight"),
)
def interpret(summary) -> str:
    """Return a plain-English interpretation of a GWAS summary."""
    lines: list[str] = []
    res = summary.primary

    lines.append(
        f"This GWAS tested {summary.n_snps:,} variants across "
        f"{summary.n_chromosomes} chromosome(s), with significance called by "
        f"the {res.label}."
    )

    if res.n_significant > 0:
        lead = summary.top_hits.iloc[0]
        lines.append(
            f"{res.n_significant} variant(s) were declared significant. The "
            f"strongest association is {lead['snp']} on chromosome "
            f"{lead['chr']} (p = {lead['p']:.2e}), a strong candidate for "
            f"follow-up and replication."
        )
    else:
        lines.append(
            "No variant passed the chosen significance criterion. This often "
            "reflects limited statistical power (small sample size) rather than "
            "the absence of any true effect — consider a larger cohort or a "
            "less conservative method (e.g. FDR)."
        )

    lam = summary.lambda_gc
    if lam < 1.05:
        lines.append(
            f"The genomic inflation factor (λ = {lam:.3f}) is close to 1.0, "
            f"suggesting well-calibrated statistics with little evidence of "
            f"population stratification or systematic bias."
        )
    elif lam < 1.2:
        lines.append(
            f"The genomic inflation factor (λ = {lam:.3f}) is mildly elevated; "
            f"consider adjusting for principal components or relatedness."
        )
    else:
        lines.append(
            f"The genomic inflation factor (λ = {lam:.3f}) is high, pointing to "
            f"substantial confounding (e.g. population stratification). "
            f"Interpret hits cautiously and re-run with appropriate covariates."
        )

    return "\n".join(lines)


def is_llm_available() -> bool:
    """Placeholder hook for a future LLM backend. Always False for now."""
    return False


planned(
    Stage.REPORT,
    "pdf",
    "PDF report",
    description="Typeset PDF report bundling plots, tables, and interpretation.",
)
