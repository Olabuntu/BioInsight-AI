"""Quality-control stage — filter markers and samples before association.

All registered as PLANNED (they need genotype-level input, which arrives with
the PLINK/VCF loaders). See DESIGN.md.
"""

from __future__ import annotations

from .registry import Param, Stage, planned

planned(
    Stage.QC,
    "maf",
    "Minor-allele-frequency filter",
    description="Drop markers below a MAF cutoff (default 0.05).",
    params=(
        Param("min_maf", "float", 0.05, help="Minimum minor allele frequency."),
    ),
)
planned(
    Stage.QC,
    "missingness",
    "Marker / sample missingness",
    description="Drop markers and samples with call rate below a threshold.",
    aliases=("geno", "mind", "call_rate"),
    params=(
        Param("max_marker_missing", "float", 0.1, help="Max marker missing rate."),
        Param("max_sample_missing", "float", 0.1, help="Max sample missing rate."),
    ),
)
planned(
    Stage.QC,
    "hwe",
    "Hardy-Weinberg equilibrium",
    description="Drop markers deviating from HWE below a p-value cutoff.",
    params=(
        Param("hwe_p", "float", 1e-6, help="HWE exact-test p-value cutoff."),
    ),
)
planned(
    Stage.QC,
    "ld_prune",
    "LD pruning",
    description="Prune markers in high linkage disequilibrium (window/step/r²).",
    params=(
        Param("window", "int", 50, help="Window size (markers)."),
        Param("step", "int", 5, help="Step size (markers)."),
        Param("r2", "float", 0.2, help="r² threshold."),
    ),
)
planned(
    Stage.QC,
    "relatedness",
    "Relatedness / duplicate removal",
    description="Remove duplicate or closely related samples above a kinship cutoff.",
)
planned(
    Stage.QC,
    "heterozygosity",
    "Heterozygosity outliers",
    description="Flag samples with excess or deficient genome-wide heterozygosity.",
)
