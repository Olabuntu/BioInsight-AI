"""Association model stage — the engines that test markers against a trait.

None are implemented yet; they are registered as PLANNED so the catalogue of
intended methods is visible (`bioinsight list --stage model`) and the pipeline
can refer to them. See DESIGN.md for the build order.

Two broad families:
  * single-locus : test one marker at a time
  * multi-locus  : fit many markers jointly (better power, fewer false positives)
"""

from __future__ import annotations

from .registry import Param, Stage, planned

# --- Single-locus models ---------------------------------------------------

planned(
    Stage.MODEL,
    "glm",
    "GLM (general linear model)",
    description="Single-locus fixed-effect model; phenotype ~ marker + covariates "
    "(e.g. PCs). Linear for quantitative traits, logistic for binary.",
    aliases=("linear", "logistic", "lm"),
    params=(
        Param("trait", "str", "quantitative",
              choices=("quantitative", "binary"), help="Trait type."),
        Param("n_pcs", "int", 0, help="Number of PCs to include as covariates."),
    ),
    reference="Price et al. 2006, Nat Genet 38:904-909 (EIGENSTRAT)",
)
planned(
    Stage.MODEL,
    "mlm",
    "MLM (mixed linear model)",
    description="Single-locus model adding a random polygenic effect via a kinship "
    "matrix to control relatedness and structure (Q+K).",
    aliases=("emmax", "q+k"),
    params=(
        Param("kinship", "str", "vanraden",
              choices=("vanraden", "ibs", "identity"), help="Kinship estimator."),
        Param("n_pcs", "int", 0, help="Number of PCs as fixed covariates (Q)."),
    ),
    reference="Yu et al. 2006, Nat Genet 38:203-208; "
    "Kang et al. 2010, Nat Genet 42:348-354 (EMMAX)",
)
planned(
    Stage.MODEL,
    "cmlm",
    "Compressed MLM",
    description="MLM with individuals clustered into groups to speed up and "
    "stabilize variance-component estimation (CMLM + P3D).",
    reference="Zhang et al. 2010, Nat Genet 42:355-360 "
    "(GAPIT software: Lipka et al. 2012, Bioinformatics)",
)
planned(
    Stage.MODEL,
    "farmcpu",
    "FarmCPU",
    description="Iterates between a fixed-effect marker model and a random-effect "
    "kinship model to control false positives without over-correction.",
    params=(
        Param("max_iterations", "int", 10, help="Maximum FarmCPU iterations."),
        Param("bin_sizes", "str", "5e5,5e6,5e7",
              help="Comma-separated candidate bin sizes (bp)."),
    ),
    reference="Liu et al. 2016, PLoS Genet 12(2):e1005767",
)
planned(
    Stage.MODEL,
    "blink",
    "BLINK",
    description="Bayesian-information and LD-based iterative model; FarmCPU "
    "successor that replaces the kinship step for speed.",
    reference="Huang et al. 2019, GigaScience 8(2):giy154",
)

# --- Multi-locus models ----------------------------------------------------

planned(
    Stage.MODEL,
    "mrmlm",
    "mrMLM (multi-locus random-SNP-effect MLM)",
    description="Two-stage multi-locus method: scan, then jointly estimate "
    "selected markers' effects under a LASSO-style model.",
    aliases=("mr-mlm",),
    params=(
        Param("lod", "float", 3.0, help="LOD threshold for the second stage."),
        Param("scan_p", "float", 0.005, help="P-value cutoff for stage-1 selection."),
    ),
    reference="Wang et al. 2016, Sci Rep 6:19444",
)
planned(
    Stage.MODEL,
    "fastmremma",
    "FASTmrEMMA",
    description="Fast multi-locus random-effect EMMA; the faster mrMLM-family "
    "method for large marker sets.",
    aliases=("fastmrmlm",),
    reference="Wen et al. 2018, Brief Bioinform 19(4):700-712",
)
planned(
    Stage.MODEL,
    "plarmeb",
    "pLARmEB",
    description="Polygenic-background-controlled LARS multi-locus method with "
    "empirical-Bayes effect estimation.",
    reference="Zhang et al. 2017, Heredity 118(6):517-524",
)
planned(
    Stage.MODEL,
    "isis_em_blasso",
    "ISIS EM-BLASSO",
    description="Iterative sure-independence screening followed by EM-Bayesian "
    "LASSO for multi-locus selection.",
    reference="Tamba et al. 2017, PLoS Comput Biol 13(1):e1005357",
)
planned(
    Stage.MODEL,
    "bayes",
    "Bayesian regression (BayesA/B/C, BSLMM)",
    description="Whole-genome Bayesian models that estimate all marker effects "
    "jointly with shrinkage priors.",
    reference="BayesA/B: Meuwissen et al. 2001, Genetics 157:1819-1829; "
    "BayesC: Habier et al. 2011, BMC Bioinformatics 12:186; "
    "BSLMM: Zhou et al. 2013, PLoS Genet 9(2):e1003264",
)
