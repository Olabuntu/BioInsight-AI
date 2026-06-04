# BioInsight — Pipeline Design

This document describes the full target architecture: a **modular GWAS
pipeline** where every capability is a pluggable *component* and the user
composes an analysis by picking one component per *stage*. Some components are
implemented today; the rest are registered as **planned** so the design is
complete and inspectable (`bioinsight list`) while we build incrementally.

## Design goals

1. **Pick what you need.** Each stage offers interchangeable methods. Want a
   GLM scan with 5 PCs, Bonferroni correction, and a QQ plot? Or an MLM with a
   VanRaden kinship and FDR control? You select; the components fit together.
2. **One contract per stage.** Components in a stage share an input/output
   contract, so they are swappable without touching the rest of the pipeline.
3. **Self-describing inputs.** Every component declares its parameters
   (`Param`: name, type, default, required, choices, help). The CLI, the config
   loader, and the docs all read from the same declarations — no method-specific
   inputs are hidden.
4. **Honest status.** `available` components run now; `planned` components raise
   a clear "on the roadmap" error instead of failing obscurely.
5. **Validated correctness.** Implemented statistics are cross-checked against
   reference implementations and the primary literature (see `tests/` and the
   citations attached to each component).

## Stages

The pipeline is an ordered sequence of stages. A run uses a subset of them
depending on the input you start from.

```
 ┌──────────┐   ┌──────┐   ┌────────────┐   ┌────────────┐   ┌──────────────┐   ┌─────────┐   ┌────────┐
 │  INPUT   │──▶│  QC  │──▶│ STRUCTURE  │──▶│   MODEL    │──▶│  THRESHOLD   │──▶│   VIZ   │──▶│ REPORT │
 │ load     │   │filter│   │ PCA/kinship│   │ assoc test │   │ significance │   │ plots   │   │ md/pdf │
 └──────────┘   └──────┘   └────────────┘   └────────────┘   └──────────────┘   └─────────┘   └────────┘
        genotype path: INPUT → QC → STRUCTURE → MODEL → THRESHOLD → VIZ → REPORT
        results  path:                       INPUT(results) → THRESHOLD → VIZ → REPORT   ← works today
```

| Stage | Purpose | Output contract |
| --- | --- | --- |
| **input** | Load genotypes/phenotypes or precomputed results | normalized table / genotype+pheno handles |
| **qc** | Filter markers & samples | filtered genotype set |
| **structure** | Estimate covariates for confounding | PCs / kinship matrix / Q matrix |
| **model** | Test markers vs. trait | results table (snp, chr, bp, p, effect) |
| **threshold** | Decide significance | per-marker significant flag + adjusted p |
| **viz** | Plot | image files |
| **report** | Summaries & interpretation | Markdown / PDF |

## Component catalogue

Legend: ✓ available · planned. Run `bioinsight list -v` for parameters and
citations, or `bioinsight describe <stage> <key>` for one component.

### input
- ✓ `results` — precomputed results table (CSV/TSV), auto-detected columns
- · `plink` — PLINK1 `.bed/.bim/.fam` / PLINK2 `.pgen`
- · `vcf` — VCF/BCF + phenotype
- · `hapmap` — HapMap genotypes (GAPIT/TASSEL)
- · `numeric` — numeric genotype matrix + phenotype/covariates

### qc (all planned)
- · `maf`, `missingness`, `hwe`, `ld_prune`, `relatedness`, `heterozygosity`

### structure (all planned)
- · `pca` (Price et al. 2006) · `kinship` (VanRaden 2008) · `qmatrix` (admixture)

### model (all planned)
Single-locus: · `glm` (Price 2006) · `mlm`/EMMAX (Yu 2006; Kang 2010)
· `cmlm` (Zhang 2010) · `farmcpu` (Liu 2016) · `blink` (Huang 2019)
Multi-locus: · `mrmlm` (Wang 2016) · `fastmremma` (Wen 2018)
· `plarmeb` (Zhang 2017) · `isis_em_blasso` (Tamba 2017) · `bayes` (Meuwissen 2001; Zhou 2013)

### threshold (all available)
- ✓ `genomewide` (Dudbridge & Gusnanto 2008; Pe'er 2008)
- ✓ `suggestive`
- ✓ `bonferroni` (Dunn 1961)
- ✓ `sidak` (Šidák 1967)
- ✓ `holm` (Holm 1979)
- ✓ `fdr_bh` (Benjamini & Hochberg 1995)
- ✓ `fdr_by` (Benjamini & Yekutieli 2001)

### viz
- ✓ `manhattan` · ✓ `qq` (with λ) · · `regional` · · `ld_heatmap` · · `effects`

### report
- ✓ `interpret` — AI-assisted interpretation (rule-based now; LLM backend planned)
- · `pdf`

## Configuration (planned)

Beyond the `gwas` quick command, a full run will be describable by a YAML
config so analyses are reproducible and shareable:

```yaml
# pipeline.yaml (planned)
input:     { plink: { prefix: data/maize, pheno: data/trait.tsv } }
qc:        { maf: { min_maf: 0.05 }, missingness: { max_marker_missing: 0.1 } }
structure: { pca: { n_components: 5 }, kinship: { method: vanraden } }
model:     { farmcpu: { max_iterations: 10 } }
threshold: { fdr_bh: { alpha: 0.05 } }
viz:       [ manhattan, qq ]
report:    [ interpret, pdf ]
```

```bash
bioinsight run --config pipeline.yaml      # planned
```

## Roadmap (build order)

1. **(done)** Registry + typed params, threshold methods, Manhattan/QQ,
   results loader, AI-assisted interpretation, Markdown report.
2. PLINK/VCF loaders + numeric genotype matrix.
3. QC filters (MAF, missingness, HWE).
4. Structure: PCA + VanRaden kinship.
5. Models: GLM → MLM/EMMAX → CMLM → FarmCPU/BLINK → multi-locus family.
6. YAML config runner (`bioinsight run`) + PDF report + LLM interpretation backend.

## Adding a component

```python
from bioinsight.registry import Stage, method, Param

@method(
    Stage.THRESHOLD, "mymethod", "My correction",
    description="...",
    params=(Param("alpha", "float", 0.05, help="Error rate."),),
    reference="Author YEAR, Journal",
)
def mymethod(pvals, alpha=0.05, **_):
    ...  # return a ThresholdResult
```

It immediately appears in `bioinsight list`, `describe`, and is selectable via
`-m mymethod`. Planned components use `planned(...)` with no handler.
