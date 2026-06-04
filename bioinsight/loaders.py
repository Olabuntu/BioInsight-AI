"""Data import stage — load association data from various sources.

Implemented today:
    * results  — a precomputed GWAS results table (CSV/TSV of p-values)

Planned (see DESIGN.md):
    * plink    — PLINK1 .bed/.bim/.fam or PLINK2 .pgen genotype sets
    * vcf      — VCF / BCF variant calls (+ a phenotype file)
    * hapmap   — HapMap-format genotypes (GAPIT/TASSEL convention)
    * numeric  — numeric genotype matrix + phenotype/covariate tables
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .registry import Param, Stage, method, planned

# Canonical column name -> accepted aliases (lowercase) for a results table.
_COLUMN_ALIASES = {
    "snp": ["snp", "rsid", "rs", "marker", "variant", "id", "name"],
    "chr": ["chr", "chrom", "chromosome", "#chrom", "linkage_group"],
    "bp": ["bp", "pos", "position", "base_pair", "ps", "coordinate"],
    "p": ["p", "pval", "p_value", "pvalue", "p.value", "p-value", "pval_nominal"],
}


def _resolve_columns(df: pd.DataFrame) -> dict:
    lowered = {c.lower(): c for c in df.columns}
    resolved = {}
    for canonical, aliases in _COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in lowered:
                resolved[canonical] = lowered[alias]
                break
    return resolved


@method(
    Stage.IO,
    "results",
    "Precomputed GWAS results",
    description="Load a results table (CSV/TSV) with at least a p-value column. "
    "Column names are auto-detected (snp/chr/bp/p and common aliases).",
    aliases=("csv", "tsv", "summary", "sumstats"),
    params=(
        Param("path", "path", required=True, help="Path to the results file (.csv/.tsv)."),
    ),
)
def load_results(path: str | Path) -> pd.DataFrame:
    """Load and normalize a GWAS results table.

    The returned frame is guaranteed to have columns: snp, chr, bp, p.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GWAS file not found: {path}")

    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    df = pd.read_csv(path, sep=sep)

    resolved = _resolve_columns(df)
    if "p" not in resolved:
        raise ValueError(
            "Could not find a p-value column. Looked for one of: "
            f"{_COLUMN_ALIASES['p']}. Found columns: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["p"] = pd.to_numeric(df[resolved["p"]], errors="coerce")
    out["snp"] = (
        df[resolved["snp"]].astype(str)
        if "snp" in resolved
        else [f"snp{i}" for i in range(len(df))]
    )
    out["chr"] = df[resolved["chr"]].astype(str) if "chr" in resolved else "NA"
    out["bp"] = (
        pd.to_numeric(df[resolved["bp"]], errors="coerce")
        if "bp" in resolved
        else range(len(df))
    )

    out = out.dropna(subset=["p"])
    out = out[(out["p"] > 0) & (out["p"] <= 1)]
    if out.empty:
        raise ValueError("No valid p-values (0 < p <= 1) found in the file.")
    return out.reset_index(drop=True)


# --- Planned loaders -------------------------------------------------------

planned(
    Stage.IO,
    "plink",
    "PLINK genotypes",
    description="PLINK1 (.bed/.bim/.fam) or PLINK2 (.pgen/.pvar/.psam) genotype sets.",
    aliases=("bed", "bfile", "pgen"),
    params=(
        Param("prefix", "str", required=True, help="PLINK fileset prefix."),
        Param("pheno", "path", help="Phenotype file."),
        Param("pheno_col", "str", help="Phenotype column name."),
    ),
)
planned(
    Stage.IO,
    "vcf",
    "VCF / BCF variants",
    description="Variant calls in VCF/BCF, paired with a phenotype table.",
    aliases=("bcf",),
    params=(
        Param("path", "path", required=True, help="Path to .vcf/.vcf.gz/.bcf."),
        Param("pheno", "path", help="Phenotype file."),
        Param("pheno_col", "str", help="Phenotype column name."),
    ),
)
planned(
    Stage.IO,
    "hapmap",
    "HapMap genotypes",
    description="HapMap-format genotypes (GAPIT / TASSEL convention).",
    aliases=("hmp",),
)
planned(
    Stage.IO,
    "numeric",
    "Numeric genotype matrix",
    description="Numeric (0/1/2) genotype matrix with phenotype and covariate tables.",
)
