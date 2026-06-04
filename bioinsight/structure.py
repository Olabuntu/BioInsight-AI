"""Population-structure stage — covariates that control confounding.

All registered as PLANNED (they need genotype-level input). See DESIGN.md.
The outputs (PCs, kinship matrix, Q matrix) feed the association models as
fixed or random effects.
"""

from __future__ import annotations

from .registry import Param, Stage, planned

planned(
    Stage.STRUCTURE,
    "pca",
    "Principal components",
    description="Top principal components of the genotype matrix, used as fixed "
    "covariates to correct for population stratification.",
    aliases=("pc", "eigenstrat"),
    params=(
        Param("n_components", "int", 5, help="Number of PCs to compute."),
    ),
    reference="Price et al. 2006, Nat Genet 38:904-909 (EIGENSTRAT)",
)
planned(
    Stage.STRUCTURE,
    "kinship",
    "Kinship / GRM",
    description="Genomic relationship (kinship) matrix used as the random-effect "
    "covariance in mixed models.",
    aliases=("grm", "vanraden"),
    params=(
        Param(
            "method", "str", "vanraden",
            choices=("vanraden", "astle-balding", "ibs"),
            help="Kinship estimator.",
        ),
    ),
    reference="VanRaden 2008, J Dairy Sci 91:4414-4423",
)
planned(
    Stage.STRUCTURE,
    "qmatrix",
    "Q matrix (admixture)",
    description="Ancestry/admixture proportions (STRUCTURE/ADMIXTURE-style) used "
    "as fixed covariates.",
    aliases=("admixture", "q"),
    params=(
        Param("k", "int", required=True, help="Number of ancestral subpopulations."),
    ),
)
