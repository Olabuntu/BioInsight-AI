"""Command-line interface for BioInsight.

    bioinsight list                          # every pipeline component + status
    bioinsight list --stage threshold -v     # ... with parameters + references
    bioinsight describe threshold fdr_bh     # what a method does + what to pass
    bioinsight gwas results.csv              # analyze (default genome-wide)
    bioinsight gwas results.csv -m fdr_bh -P alpha=0.1
    bioinsight gwas results.csv -m genomewide -P threshold=1e-6 -p manhattan,qq
"""

from __future__ import annotations

import argparse
import sys

from . import __version__, pipeline, registry
from .registry import Stage, Status


def _print_params(comp, indent="      "):
    if not comp.params:
        print(f"{indent}(no parameters)")
        return
    for p in comp.params:
        bits = [f"type={p.type}"]
        if p.default is not None:
            bits.append(f"default={p.default}")
        if p.required:
            bits.append("required")
        if p.choices:
            bits.append("choices=" + "|".join(map(str, p.choices)))
        print(f"{indent}{p.name:<18} {p.help}  [{', '.join(bits)}]")


def _cmd_list(args: argparse.Namespace) -> int:
    stage = Stage(args.stage) if args.stage else None
    stages = [stage] if stage else list(Stage)
    for s in stages:
        comps = registry.list_components(s)
        if not comps:
            continue
        print(f"\n{s.title}  (stage: {s.value})")
        print("-" * 64)
        for c in sorted(comps, key=lambda c: (c.status.value, c.key)):
            mark = "✓" if c.status is Status.AVAILABLE else "·"
            tag = "" if c.status is Status.AVAILABLE else "  [planned]"
            print(f"  {mark} {c.key:<16} {c.name}{tag}")
            if args.verbose:
                if c.description:
                    print(f"      {c.description}")
                if c.params:
                    _print_params(c)
                if c.reference:
                    print(f"      ref: {c.reference}")
    print("\n✓ = available now    · = planned (see DESIGN.md)")
    return 0


def _cmd_describe(args: argparse.Namespace) -> int:
    comp = registry.get(Stage(args.stage), args.key)
    status = "available" if comp.available else "planned"
    print(f"\n{comp.name}  ({args.stage}:{comp.key})  [{status}]")
    print("-" * 64)
    if comp.description:
        print(comp.description)
    if comp.aliases:
        print(f"\naliases: {', '.join(comp.aliases)}")
    print("\nparameters:")
    _print_params(comp, indent="  ")
    if comp.reference:
        print(f"\nreference: {comp.reference}")
    if comp.params:
        example = " ".join(
            f"-P {p.name}={p.default if p.default is not None else '<value>'}"
            for p in comp.params
        )
        print(f"\nexample: bioinsight gwas results.csv -m {comp.key} {example}")
    return 0


def _parse_csv(value: str):
    return tuple(p.strip() for p in value.split(",") if p.strip())


def _parse_params(pairs):
    """Turn ['alpha=0.1', 'trait=binary'] into {'alpha': '0.1', ...}."""
    out = {}
    for item in pairs or []:
        if "=" not in item:
            raise ValueError(f"--param must be KEY=VALUE; got {item!r}.")
        key, _, val = item.partition("=")
        out[key.strip()] = val.strip()
    return out


def _cmd_gwas(args: argparse.Namespace) -> int:
    params = _parse_params(args.param)
    # Convenience flags fold into the param dict (explicit --param wins).
    if args.alpha is not None:
        params.setdefault("alpha", args.alpha)
    if args.threshold is not None:
        params.setdefault("threshold", args.threshold)

    print(f"Loading results from {args.file} ...")
    print(f"Significance method: {args.method}")
    result = pipeline.run_results_analysis(
        args.file,
        method=args.method,
        outdir=args.outdir,
        plots=_parse_csv(args.plots),
        top_n=args.top_n,
        **params,
    )

    summary = result["summary"]
    res = summary.primary
    print()
    print("=" * 64)
    print(f"  Variants tested:      {summary.n_snps:,}")
    print(f"  Significance method:  {res.label}")
    print(f"  Significant markers:  {res.n_significant}")
    print(f"  Genomic inflation λ:  {summary.lambda_gc:.3f}")
    print("=" * 64)
    if summary.method_comparison:
        print("\n  Significant markers by method:")
        for label, n in summary.method_comparison.items():
            print(f"    {n:>6}  {label}")
    print()
    print(result["interpretation"])
    print()
    print(f"Report written to: {result['report_path']}")
    for name, path in result["plot_paths"].items():
        print(f"Plot written to:   {path}  ({name})")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioinsight",
        description="A modular GWAS pipeline — pick the tools you need.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_list = sub.add_parser("list", help="List pipeline components and status.")
    p_list.add_argument(
        "--stage", choices=[s.value for s in Stage], help="Restrict to one stage."
    )
    p_list.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show descriptions, parameters, and references.",
    )
    p_list.set_defaults(func=_cmd_list)

    p_desc = sub.add_parser(
        "describe", help="Show a component's details and accepted parameters."
    )
    p_desc.add_argument("stage", choices=[s.value for s in Stage])
    p_desc.add_argument("key", help="Component key or alias.")
    p_desc.set_defaults(func=_cmd_describe)

    g = sub.add_parser("gwas", help="Analyze a GWAS results file.")
    g.add_argument("file", help="Path to GWAS results (CSV or TSV).")
    g.add_argument(
        "-m", "--method", default="genomewide",
        help="Significance method (see `bioinsight list --stage threshold`).",
    )
    g.add_argument(
        "-P", "--param", action="append", metavar="KEY=VALUE",
        help="Method-specific parameter; repeatable. See `bioinsight describe "
        "threshold <method>`.",
    )
    g.add_argument(
        "-a", "--alpha", type=float, default=None,
        help="Shortcut for -P alpha=… (correction methods).",
    )
    g.add_argument(
        "-t", "--threshold", type=float, default=None,
        help="Shortcut for -P threshold=… (fixed-threshold methods).",
    )
    g.add_argument(
        "-o", "--outdir", default="bioinsight_report", help="Output directory."
    )
    g.add_argument(
        "-p", "--plots", default="manhattan,qq",
        help="Comma-separated plots to render (e.g. manhattan,qq).",
    )
    g.add_argument(
        "-n", "--top-n", type=int, default=10, help="Number of top hits to report."
    )
    g.set_defaults(func=_cmd_gwas)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (FileNotFoundError, ValueError, KeyError, NotImplementedError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
