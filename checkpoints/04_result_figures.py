"""Checkpoint 4: the five result figures.

Run this if you are behind in Exercise 02 and want the result figures immediately.
It runs checkpoint 3 first if results/de_results.tsv is missing.

    python checkpoints/04_result_figures.py

It writes five figures to results/figures/. It does not run the enrichment
test; that is one command, given in the Exercise 02 handout.
"""

import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIG = RESULTS / "figures"
FIG.mkdir(parents=True, exist_ok=True)

GREY, ORANGE = "#999999", "#D55E00"
PADJ, LFC = 0.05, 1.0

if not (RESULTS / "de_results.tsv").exists():
    print("results/de_results.tsv not found. Running checkpoint 3 first.\n")
    subprocess.run([sys.executable, str(ROOT / "checkpoints" /
                                        "03_differential_expression.py")], check=True)
    print()

results = pd.read_csv(RESULTS / "de_results.tsv", sep="\t")
normalized = pd.read_csv(RESULTS / "normalized_counts.tsv", sep="\t")
meta = pd.read_csv(ROOT / "data" / "sample_metadata.tsv", sep="\t")
samples = list(meta["sample_id"])

called = (results["padj"] < PADJ) & (results["log2_fold_change"].abs() >= LFC)


def finish(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=11, loc="left")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)


# Figure 4: distribution of unadjusted p-values.
fig, ax = plt.subplots(figsize=(5.4, 3.6))
ax.hist(results["pvalue"], bins=40, range=(0, 1), color="#0072B2", edgecolor="white")
expected = len(results) / 40
ax.axhline(expected, color=GREY, linestyle="--", linewidth=1)
ax.text(0.30, expected * 1.35, "expected if no gene were differential",
        fontsize=8, color=GREY, bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
finish(ax, f"Distribution of unadjusted p-values ({len(results)} genes tested)",
       "p-value", "Number of genes")
fig.savefig(FIG / "04_pvalue_histogram.png", dpi=200, bbox_inches="tight")

# Figure 5: MA plot.
fig, ax = plt.subplots(figsize=(5.8, 4.2))
ax.scatter(results.loc[~called, "mean_expression"],
           results.loc[~called, "log2_fold_change"],
           s=5, c=GREY, alpha=0.4, linewidths=0)
ax.scatter(results.loc[called, "mean_expression"],
           results.loc[called, "log2_fold_change"],
           s=7, c=ORANGE, linewidths=0)
ax.axhline(0, color="black", linewidth=0.8)
finish(ax, f"MA plot. {int(called.sum())} genes pass padj < {PADJ} "
           f"and |log2 fold change| >= {LFC}",
       "Mean expression (log2 normalized counts)",
       "log2 fold change (treated / control)")
fig.savefig(FIG / "05_ma_plot.png", dpi=200, bbox_inches="tight")

# Figure 6: volcano plot.
y = -np.log10(results["pvalue"].clip(lower=1e-300))
fig, ax = plt.subplots(figsize=(5.8, 4.6))
ax.scatter(results.loc[~called, "log2_fold_change"], y[~called],
           s=5, c=GREY, alpha=0.4, linewidths=0)
ax.scatter(results.loc[called, "log2_fold_change"], y[called],
           s=8, c=ORANGE, linewidths=0)
ax.axvline(LFC, color="black", linestyle="--", linewidth=0.8)
ax.axvline(-LFC, color="black", linestyle="--", linewidth=0.8)
cutoff = results.loc[results["padj"] < PADJ, "pvalue"].max()
if pd.notna(cutoff):
    ax.axhline(-np.log10(cutoff), color="black", linestyle=":", linewidth=0.8)
    ax.text(ax.get_xlim()[0], -np.log10(cutoff) * 1.02, f"padj = {PADJ}",
            fontsize=8, va="bottom")
for direction, side in [(1, 1), (-1, -1)]:
    subset = results[called & (np.sign(results["log2_fold_change"]) == direction)]
    for rank, (_, row) in enumerate(subset.head(5).iterrows()):
        ax.annotate(row["gene_symbol"],
                    (row["log2_fold_change"], -np.log10(max(row["pvalue"], 1e-300))),
                    fontsize=7, xytext=(7 * side, 9 - 6 * rank),
                    textcoords="offset points",
                    ha="left" if side > 0 else "right")
finish(ax, "Volcano plot. Dashed lines mark the thresholds (n = 8 per condition)",
       "log2 fold change (treated / control)", "-log10 unadjusted p-value")
fig.savefig(FIG / "06_volcano.png", dpi=200, bbox_inches="tight")

# Figure 7: heatmap of the 40 genes with the smallest adjusted p-value.
top = results.head(40)
block = normalized.set_index("gene_id").loc[top["gene_id"], samples].to_numpy(dtype=float)
block = np.log2(block + 1)
z = (block - block.mean(axis=1, keepdims=True)) / block.std(axis=1, keepdims=True)
order = np.argsort(meta["condition"].to_numpy(), kind="stable")

fig, ax = plt.subplots(figsize=(7.2, 7.6))
image = ax.imshow(z[:, order], cmap="RdBu_r",
                  norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5), aspect="auto")
ax.set_yticks(range(len(top)), top["gene_symbol"], fontsize=7)
ax.set_xticks(range(len(order)),
              [f"{samples[i]} {meta['condition'][i][:4]}" for i in order],
              rotation=90, fontsize=7)
bar = fig.colorbar(image, ax=ax, shrink=0.5)
bar.set_label("Row z-score of log2 normalized counts", fontsize=9)
ax.set_title("40 genes with the smallest adjusted p-value\n"
             "Rows scaled to mean 0 and standard deviation 1. n = 8 per condition.",
             fontsize=11, loc="left")
fig.savefig(FIG / "07_heatmap_top_genes.png", dpi=200, bbox_inches="tight")

print("Wrote:")
for name in ["04_pvalue_histogram.png", "05_ma_plot.png", "06_volcano.png",
             "07_heatmap_top_genes.png"]:
    print(f"  results/figures/{name}")
print()
print("The enrichment figure is not made here. Run the enrichment test first:")
print("  python scripts/enrich.py results/de_results.tsv --out results/enrichment.tsv")
print("then ask Claude to plot the result. Check the figure against")
print("docs/figure-standards.md before you accept it.")
