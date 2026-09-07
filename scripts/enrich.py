"""Gene set enrichment by the hypergeometric test.

This script is provided so that the enrichment step takes one command rather
than a programming exercise. It is short on purpose: read it before you use it.

The test asks a single question for each gene set. Given N genes tested in
total, of which K are in the gene set, and given that you selected n genes as
significant, of which k are in the gene set, how surprising is k? The p-value is
the probability of observing k or more members by chance.

P-values are adjusted across gene sets with the Benjamini-Hochberg procedure.

Usage:
    python scripts/enrich.py results/de_results.tsv
    python scripts/enrich.py results/de_results.tsv --padj 0.05 --lfc 1.0
    python scripts/enrich.py results/de_results.tsv --out results/enrichment.tsv

The input file must contain the columns gene_symbol, log2_fold_change and
padj. Column names are matched case-insensitively and common variants are
accepted, for example log2FoldChange or adj_pvalue.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_GMT = ROOT / "data" / "gene_sets.gmt"

COLUMN_ALIASES = {
    "gene_symbol": ["gene_symbol", "gene", "symbol", "genesymbol", "gene_name"],
    "log2_fold_change": [
        "log2_fold_change", "log2foldchange", "log2fc", "lfc", "logfc",
        "log2_fc", "log2fold",
    ],
    "padj": [
        "padj", "adj_pvalue", "adj_p", "fdr", "qvalue", "q_value",
        "adjusted_pvalue", "p_adj", "pvalue_adj", "bh",
    ],
}


def read_gmt(path):
    """Read a GMT file into a dictionary of set name to member list."""
    sets = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            sets[fields[0]] = [g for g in fields[2:] if g]
    return sets


def resolve_columns(frame):
    """Map the columns of the input table onto the names this script needs."""
    lookup = {c.lower().replace(" ", "_"): c for c in frame.columns}
    resolved = {}
    for wanted, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in lookup:
                resolved[wanted] = lookup[alias]
                break
        if wanted not in resolved:
            raise SystemExit(
                f"The input file has no column for '{wanted}'.\n"
                f"Columns found: {list(frame.columns)}\n"
                f"Accepted names: {COLUMN_ALIASES[wanted]}"
            )
    return resolved


def benjamini_hochberg(pvalues):
    """Return Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


def enrich(results, gene_sets, padj_threshold=0.05, lfc_threshold=1.0,
           min_set_size=10):
    """Run the hypergeometric test for every gene set."""
    cols = resolve_columns(results)
    table = results[[cols["gene_symbol"], cols["log2_fold_change"], cols["padj"]]]
    table.columns = ["gene_symbol", "log2_fold_change", "padj"]
    table = table.dropna(subset=["gene_symbol"])

    universe = set(table["gene_symbol"])
    significant = set(
        table.loc[
            (table["padj"] < padj_threshold)
            & (table["log2_fold_change"].abs() >= lfc_threshold),
            "gene_symbol",
        ]
    )

    if not significant:
        raise SystemExit(
            "No genes passed the thresholds, so there is nothing to test.\n"
            f"Thresholds used: padj < {padj_threshold}, "
            f"absolute log2 fold change >= {lfc_threshold}."
        )

    n_universe = len(universe)
    n_selected = len(significant)

    rows = []
    for name, members in gene_sets.items():
        in_universe = set(members) & universe
        if len(in_universe) < min_set_size:
            continue
        hits = in_universe & significant
        k = len(hits)
        K = len(in_universe)
        # Survival function at k-1 gives P(X >= k).
        pvalue = hypergeom.sf(k - 1, n_universe, K, n_selected)
        expected = n_selected * K / n_universe
        rows.append(
            {
                "gene_set": name,
                "set_size": K,
                "n_significant_in_set": k,
                "expected_by_chance": round(expected, 2),
                "fold_enrichment": round(k / expected, 3) if expected > 0 else np.nan,
                "pvalue": pvalue,
                "genes": ";".join(sorted(hits)),
            }
        )

    out = pd.DataFrame(rows).sort_values("pvalue").reset_index(drop=True)
    out["padj"] = benjamini_hochberg(out["pvalue"])
    out = out[
        [
            "gene_set", "set_size", "n_significant_in_set", "expected_by_chance",
            "fold_enrichment", "pvalue", "padj", "genes",
        ]
    ]
    return out, n_universe, n_selected


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("results", help="differential expression results table")
    parser.add_argument("--gmt", default=str(DEFAULT_GMT), help="gene set file")
    parser.add_argument("--padj", type=float, default=0.05,
                        help="adjusted p-value threshold for calling a gene significant")
    parser.add_argument("--lfc", type=float, default=1.0,
                        help="minimum absolute log2 fold change")
    parser.add_argument("--out", default=None, help="output file (optional)")
    parser.add_argument("--top", type=int, default=10,
                        help="number of gene sets to print")
    args = parser.parse_args()

    sep = "\t" if args.results.endswith((".tsv", ".txt")) else ","
    results = pd.read_csv(args.results, sep=sep)
    gene_sets = read_gmt(args.gmt)

    table, n_universe, n_selected = enrich(
        results, gene_sets, padj_threshold=args.padj, lfc_threshold=args.lfc
    )

    print(f"Genes tested:      {n_universe}")
    print(f"Genes selected:    {n_selected} "
          f"(padj < {args.padj}, |log2 fold change| >= {args.lfc})")
    print(f"Gene sets tested:  {len(table)}")
    print()
    shown = table.head(args.top).copy()
    shown["pvalue"] = shown["pvalue"].map(lambda v: f"{v:.2e}")
    shown["padj"] = shown["padj"].map(lambda v: f"{v:.2e}")
    print(shown.drop(columns="genes").to_string(index=False))
    print()
    print("A small p-value means the gene set contains more of your significant")
    print("genes than would be expected by chance. It does not tell you that the")
    print("biological process is active. Read data/README.md before you draw a")
    print("biological conclusion from this table.")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.out, sep="\t", index=False)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
