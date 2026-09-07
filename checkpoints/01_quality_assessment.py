"""Checkpoint 1: data loaded and quality assessed.

Run this if you are behind and want to reach the end of step 2 of the analysis
path immediately.

    python checkpoints/01_quality_assessment.py

It loads the count matrix and the sample table, prints the quality summary, and
writes three figures to results/figures/.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE = "#0072B2", "#D55E00"

counts = pd.read_csv(ROOT / "data" / "counts.tsv", sep="\t")
meta = pd.read_csv(ROOT / "data" / "sample_metadata.tsv", sep="\t")
samples = list(meta["sample_id"])
matrix = counts[samples].to_numpy(dtype=float)

print(f"Genes: {matrix.shape[0]}")
print(f"Samples: {matrix.shape[1]}")
print(meta.groupby(["condition", "batch"]).size().to_string())
print()

library_size = matrix.sum(axis=0)
detected = (matrix > 0).sum(axis=0)
summary = pd.DataFrame({
    "sample_id": samples,
    "condition": meta["condition"],
    "batch": meta["batch"],
    "library_size": library_size.astype(int),
    "genes_detected": detected,
})
print(summary.to_string(index=False))

# Figure 1: library size per sample.
colors = [BLUE if c == "control" else ORANGE for c in meta["condition"]]
fig, ax = plt.subplots(figsize=(7, 3.4))
ax.bar(samples, library_size / 1e6, color=colors)
for i, batch in enumerate(meta["batch"]):
    ax.text(i, 0.3, batch, ha="center", fontsize=8, color="white")
ax.set_xlabel("Sample")
ax.set_ylabel("Total counts (millions)")
ax.set_title("Library size per sample (batch letter inside each bar)", loc="left")
ax.tick_params(axis="x", rotation=90)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig(FIG / "01_library_size.png", dpi=200, bbox_inches="tight")

# Log counts per million, used for the next two figures.
cpm = matrix / matrix.sum(axis=0) * 1e6
log_cpm = np.log2(cpm + 1)

# Figure 2: principal component analysis of the 1000 most variable genes.
variable = np.argsort(log_cpm.var(axis=1))[-1000:]
x = log_cpm[variable].T
x = x - x.mean(axis=0)
u, s, _ = np.linalg.svd(x, full_matrices=False)
scores = u * s
explained = s ** 2 / (s ** 2).sum() * 100

fig, ax = plt.subplots(figsize=(5.4, 4.6))
for condition, color in [("control", BLUE), ("treated", ORANGE)]:
    for batch, marker in [("A", "o"), ("B", "^")]:
        mask = ((meta["condition"] == condition) & (meta["batch"] == batch)).to_numpy()
        ax.scatter(scores[mask, 0], scores[mask, 1], c=color, marker=marker,
                   s=60, edgecolor="white", label=f"{condition}, batch {batch}")
ax.set_xlabel(f"PC1 ({explained[0]:.1f} per cent of variance)")
ax.set_ylabel(f"PC2 ({explained[1]:.1f} per cent of variance)")
ax.set_title("Principal component analysis, 1000 most variable genes\n"
             "n = 8 per condition", loc="left")
ax.legend(frameon=False, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig(FIG / "02_pca.png", dpi=200, bbox_inches="tight")

# Figure 3: sample to sample correlation.
corr = np.corrcoef(log_cpm.T)
labels = [f"{s} {c[:4]} {b}" for s, c, b in
          zip(samples, meta["condition"], meta["batch"])]
fig, ax = plt.subplots(figsize=(6.4, 5.6))
image = ax.imshow(corr, cmap="viridis")
ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=7)
ax.set_yticks(range(len(labels)), labels, fontsize=7)
bar = fig.colorbar(image, ax=ax, shrink=0.8)
bar.set_label("Pearson correlation of log2 counts per million")
ax.set_title("Sample to sample correlation", loc="left")
fig.savefig(FIG / "03_sample_correlation.png", dpi=200, bbox_inches="tight")

print()
print("Wrote results/figures/01_library_size.png")
print("Wrote results/figures/02_pca.png")
print("Wrote results/figures/03_sample_correlation.png")
print()
print("Look at the PCA figure before you continue. Which variable separates the")
print("samples along PC1, and which along PC2?")
