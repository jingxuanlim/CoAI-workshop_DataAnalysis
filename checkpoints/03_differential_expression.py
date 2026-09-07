"""Checkpoint 3: differential expression tested, results table written.

Run this if you are behind and want to reach the end of step 6 of the analysis
path immediately. It runs checkpoint 2 first if the normalized counts are
missing.

    python checkpoints/03_differential_expression.py

It writes results/de_results.tsv, which is the file that
scripts/enrich.py and scripts/score_results.py both read.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
NORMALIZED = RESULTS / "normalized_counts.tsv"

if not NORMALIZED.exists():
    print("results/normalized_counts.tsv not found. Running checkpoint 2 first.\n")
    subprocess.run([sys.executable, str(ROOT / "checkpoints" /
                                        "02_filtered_and_normalized.py")], check=True)
    print()

table = pd.read_csv(NORMALIZED, sep="\t")
meta = pd.read_csv(ROOT / "data" / "sample_metadata.tsv", sep="\t")
samples = list(meta["sample_id"])

normalized = table[samples].to_numpy(dtype=float)
log_expression = np.log2(normalized + 1)

# Step 5: test each gene for a difference between the two conditions.
# Welch's t-test on the log transformed normalized counts. It does not assume
# that the two groups have the same variance.
is_treated = (meta["condition"] == "treated").to_numpy()
treated = log_expression[:, is_treated]
control = log_expression[:, ~is_treated]

t_statistic, pvalue = stats.ttest_ind(treated, control, axis=1, equal_var=False)
log2_fold_change = treated.mean(axis=1) - control.mean(axis=1)


def benjamini_hochberg(p):
    """Adjust p-values for multiple testing.

    Almost 4000 genes are tested. At an unadjusted threshold of 0.05 about 200
    genes would be called significant even if no gene differed at all. This
    procedure controls the expected proportion of false discoveries instead.
    """
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted


padj = benjamini_hochberg(pvalue)

# Step 6: results table.
results = pd.DataFrame({
    "gene_id": table["gene_id"],
    "gene_symbol": table["gene_symbol"],
    "mean_expression": np.round(log_expression.mean(axis=1), 4),
    "log2_fold_change": np.round(log2_fold_change, 4),
    "t_statistic": np.round(t_statistic, 4),
    "pvalue": pvalue,
    "padj": padj,
}).sort_values("padj").reset_index(drop=True)

results.to_csv(RESULTS / "de_results.tsv", sep="\t", index=False)

significant = (results["padj"] < 0.05) & (results["log2_fold_change"].abs() >= 1.0)
print(f"Genes tested:                 {len(results)}")
print(f"Adjusted p-value below 0.05:  {int((results['padj'] < 0.05).sum())}")
print(f"Also |log2 fold change| >= 1: {int(significant.sum())}")
print()
print("Ten genes with the smallest adjusted p-value:")
print(results.head(10)[
    ["gene_symbol", "mean_expression", "log2_fold_change", "padj"]
].to_string(index=False))
print()
print("Wrote results/de_results.tsv")
