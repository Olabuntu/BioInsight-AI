# BioInsight

**BioInsight** is a modular GWAS pipeline. Every capability — data loaders,
quality control, population-structure correction, association models,
significance methods, plots, reports — is a pluggable component, and you
compose an analysis by **picking the tools you need**. Implemented methods work
today; the rest of the architecture is a published, inspectable roadmap (see
[DESIGN.md](DESIGN.md)).

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)
![Status](https://img.shields.io/badge/status-alpha-orange.svg)

```bash
bioinsight list                       # see every component and what's available
bioinsight describe threshold fdr_bh  # see what a method does + what to pass
bioinsight gwas results.csv -m fdr_bh -P alpha=0.1
```

---

## Why "pick what you need"

A GWAS is a sequence of stages, and at each stage there are many valid methods.
BioInsight makes each stage a slot you fill:

```
INPUT → QC → STRUCTURE → MODEL → THRESHOLD → VIZ → REPORT
```

You choose, for example, *PLINK input → MAF filter → 5 PCs → MLM → FDR → Manhattan + QQ*,
or *a results CSV → Bonferroni → QQ plot*. Components in a stage are
interchangeable because they share one input/output contract. Each method
declares its own parameters, so the tool always tells you what to enter.

> **What runs today (v0.2.0):** the results-analysis path — load a results
> table, apply any of **7 significance/correction methods**, render Manhattan
> and QQ plots, and get a report with plain-English interpretation. Genotype
> loaders (PLINK/VCF), QC, structure, and association models (GLM, MLM, CMLM,
> FarmCPU, BLINK, mrMLM, …) are registered as **planned** — visible in
> `bioinsight list`, documented in [DESIGN.md](DESIGN.md), and built next.

## Features

- 🧩 **Modular registry** — every method is a component you can list, describe, and pick
- 🎚️ **7 significance methods, validated** — genome-wide, suggestive, Bonferroni,
  Šidák, Holm, Benjamini-Hochberg & Benjamini-Yekutieli FDR
  (cross-checked against `statsmodels` to 1e-12; see [tests](tests/test_validation.py))
- 🧮 **Typed, self-describing parameters** — `bioinsight describe <stage> <method>`
  shows exactly what each method accepts, with types, defaults, and choices
- 🌋 **Manhattan & QQ plots** — publication-ready, headless-friendly, with genomic inflation (λ)
- 🤖 **AI-assisted interpretation** — plain-English summary (rule-based today, no API key; LLM backend planned)
- 📄 **Reproducible Markdown reports** bundling plots, a top-hits table, and a method-comparison table
- 📚 **Cited methods** — every component carries its primary literature reference

## Installation

```bash
git clone https://github.com/Olabuntu/BioInsight-AI.git
cd BioInsight-AI
pip install -e .
```

## Quickstart

```bash
bioinsight gwas examples/sample_gwas.csv -m fdr_bh
```

```
================================================================
  Variants tested:      304
  Significance method:  Benjamini-Hochberg (FDR=0.05)
  Significant markers:  4
  Genomic inflation λ:  0.960
================================================================

  Significant markers by method:
       2  genome-wide (p ≤ 5e-08)
       4  Bonferroni (α=0.05, p ≤ 1.64e-04)
       4  Benjamini-Hochberg (FDR=0.05)
       ...
Report written to: bioinsight_report/report.md
Plot written to:   bioinsight_report/manhattan.png
Plot written to:   bioinsight_report/qq.png
```

## Picking methods and parameters

List everything and its status:

```bash
bioinsight list                  # all stages
bioinsight list --stage threshold -v   # with parameters + citations
```

Ask what a method needs before you run it:

```bash
$ bioinsight describe threshold bonferroni
Bonferroni correction  (threshold:bonferroni)  [available]
parameters:
  alpha   Family-wise error rate.  [type=float, default=0.05]
reference: Bonferroni 1936; Dunn 1961, J Am Stat Assoc 56:52-64
example: bioinsight gwas results.csv -m bonferroni -P alpha=0.05
```

Pass method-specific parameters with `-P key=value` (repeatable):

```bash
bioinsight gwas results.csv -m genomewide -P threshold=1e-6
bioinsight gwas results.csv -m fdr_bh     -P alpha=0.1 -p manhattan,qq
```

Unknown parameters, wrong types, and invalid choices are rejected with a clear
message — the validation comes straight from each method's declared spec.

## Input format

A CSV or TSV with at least a p-value column. Column names are auto-detected
(case-insensitive, common aliases accepted):

| Field | Accepted names |
| --- | --- |
| SNP ID | `snp`, `rsid`, `marker`, `variant`, `id` |
| Chromosome | `chr`, `chrom`, `chromosome` |
| Position | `bp`, `pos`, `position` |
| P-value | `p`, `pval`, `p_value`, `pvalue` |

See [`examples/sample_gwas.csv`](examples/sample_gwas.csv).

## Use as a library

```python
from bioinsight import pipeline

result = pipeline.run_results_analysis(
    "results.csv", method="fdr_bh", alpha=0.1, plots=("manhattan", "qq")
)
print(result["interpretation"])
print(result["summary"].method_comparison)
```

## Architecture & roadmap

The full pipeline design — every stage, every planned method, the config format,
and how to add a component — is in **[DESIGN.md](DESIGN.md)**.

## Development

```bash
pip install -e ".[dev]"
pytest        # includes validation against statsmodels
```

## License

[MIT](LICENSE) © 2026 Abhishek Sahu

---

*Contributions welcome — open an issue or PR.*
