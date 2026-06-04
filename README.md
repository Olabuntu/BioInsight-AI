# BioInsight-AI

**BioInsight-AI** is an open-source bioinformatics toolkit that transforms
sequencing and omics outputs into interpretable biological insights — with
publication-ready visualizations and plain-English, AI-assisted summaries.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)
![Status](https://img.shields.io/badge/status-alpha-orange.svg)

> **v0.1.0** focuses on **GWAS analysis**. Variant summarization and
> expression analysis are on the roadmap (see [Roadmap](#roadmap)).

---

## Features

- 📊 **GWAS analysis** — load results, flag genome-wide and suggestive hits
- 🌋 **Manhattan plots** — publication-ready, headless-friendly (no display needed)
- 🧬 **Quality metrics** — genomic inflation factor (λ) to catch confounding
- 🤖 **AI-assisted interpretation** — a plain-English summary of what your
  results mean (rule-based in v0.1.0, no API key required)
- 📄 **Reproducible reports** — a single Markdown report with plot, table, and summary
- 🔌 **Flexible input** — CSV or TSV, with automatic column-name detection

## Installation

```bash
# From source (recommended during alpha)
git clone https://github.com/abhisheksahu/BioInsight-AI.git
cd BioInsight-AI
pip install -e .
```

Once published to PyPI:

```bash
pip install bioinsight-ai
```

## Quickstart

```bash
bioinsight gwas examples/sample_gwas.csv
```

This reads the GWAS results, finds significant SNPs, renders a Manhattan plot,
and writes a summary report:

```
============================================================
  Variants tested:        304
  Genome-wide significant: 2
  Suggestive:              4
  Genomic inflation (λ):   1.0xx
============================================================

This GWAS tested 304 variants across 4 chromosome(s).
2 variant(s) reached genome-wide significance (p ≤ 5e-08). The strongest
association is rs1042 on chromosome 1 (p = 3.10e-09), a strong candidate
for follow-up.
...
Report written to: bioinsight_report/report.md
Plot written to:   bioinsight_report/manhattan.png
```

### Options

```bash
bioinsight gwas results.csv \
  --outdir my_report/ \      # where to write report.md + manhattan.png
  --threshold 5e-8 \         # genome-wide significance cutoff
  --top-n 20                 # number of top hits to tabulate
```

## Input format

A CSV or TSV with (at minimum) a p-value column. Column names are detected
case-insensitively, and common aliases are accepted:

| Field | Accepted column names |
| --- | --- |
| SNP ID | `snp`, `rsid`, `marker`, `variant`, `id` |
| Chromosome | `chr`, `chrom`, `chromosome` |
| Position | `bp`, `pos`, `position` |
| P-value | `p`, `pval`, `p_value`, `pvalue` |

See [`examples/sample_gwas.csv`](examples/sample_gwas.csv) for a template.

## Use as a library

```python
from bioinsight import gwas, ai

df = gwas.load_gwas("results.csv")
result = gwas.analyze(df, threshold=5e-8)
gwas.manhattan_plot(result, "manhattan.png")
print(ai.interpret_gwas(result))
```

## Roadmap

- [x] GWAS analysis, Manhattan plot, AI summary, Markdown report
- [ ] Variant summarization (VCF: counts, types, Ts/Tv)
- [ ] Expression analysis (differential expression)
- [ ] QQ plots and additional diagnostics
- [ ] Optional LLM backend for richer interpretation

## Development

```bash
pip install -e ".[dev]"
pytest
```

## License

[MIT](LICENSE) © 2026 Abhishek Sahu

---

*Contributions welcome — open an issue or PR.*
