"""Checkpoint 2: genes filtered and counts normalized.

Run this if you are behind and want to reach the end of step 4 of the analysis
path immediately.

    python checkpoints/02_filtered_and_normalized.py

It writes results/normalized_counts.tsv, which the next checkpoint reads.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

counts = pd.read_csv(ROOT / "data" / "counts.tsv", sep="\t")
meta = pd.read_csv(ROOT / "data" / "sample_metadata.tsv", sep="\t")
samples = list(meta["sample_id"])
matrix = counts[samples].to_numpy(dtype=float)

# Step 3: filter.
# Keep a gene if it has more than one count per million in at least eight
# samples. Eight is the size of the smaller group, so a gene expressed in only
# one condition is still kept.
cpm = matrix / matrix.sum(axis=0) * 1e6
keep = (cpm > 1).sum(axis=1) >= 8
print(f"Genes before filtering: {len(keep)}")
print(f"Genes removed:          {int((~keep).sum())}")
print(f"Genes kept:             {int(keep.sum())}")

matrix = matrix[keep]
genes = counts.loc[keep, ["gene_id", "gene_symbol"]].reset_index(drop=True)

# Step 4: normalize by the median of ratios method.
# For each gene, take the ratio of its count to its geometric mean across all
# samples. The size factor for a sample is the median of those ratios. This
# method is used by DESeq2. It is not affected by a small number of very highly
# expressed genes, which is why it is preferred here over scaling by the total
# count.
with np.errstate(divide="ignore"):
    log_matrix = np.log(matrix)
usable = np.all(np.isfinite(log_matrix), axis=1)
reference = log_matrix[usable].mean(axis=1, keepdims=True)
factors = np.exp(np.median(log_matrix[usable] - reference, axis=0))

print()
print("Size factors:")
for sample, factor in zip(samples, factors):
    print(f"  {sample}  {factor:.3f}")

normalized = matrix / factors

out = pd.DataFrame(normalized, columns=samples).round(3)
out.insert(0, "gene_symbol", genes["gene_symbol"].values)
out.insert(0, "gene_id", genes["gene_id"].values)
out.to_csv(RESULTS / "normalized_counts.tsv", sep="\t", index=False)

print()
print(f"Wrote results/normalized_counts.tsv "
      f"({out.shape[0]} genes x {len(samples)} samples)")
print()
print("Compare the size factors with the library sizes from checkpoint 1.")
print("They are not proportional. Ask Claude why, and check the answer against")
print("the comment above the size factor calculation in this file.")
