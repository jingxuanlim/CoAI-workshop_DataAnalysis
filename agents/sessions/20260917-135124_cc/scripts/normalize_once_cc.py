# Execution Mode: One-shot (Single-run)
# Filter low-count genes and normalize counts using TMM.
# Writes: results/normalized_counts.tsv

import os


def load_data(counts_path, metadata_path):
    import pandas as pd
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    gene_info = counts[["gene_symbol", "gene_length"]].copy()
    matrix = counts.drop(columns=["gene_symbol", "gene_length"])
    return matrix, gene_info, metadata


def filter_low_counts(matrix, min_cpm=1.0, min_samples=8):
    """
    Keep genes with CPM >= min_cpm in at least min_samples samples.
    CPM is computed per-sample before filtering so the threshold is
    library-size-independent.
    """
    import numpy as np
    lib_sizes = matrix.sum(axis=0)
    cpm = matrix.divide(lib_sizes, axis=1) * 1e6
    passes = (cpm >= min_cpm).sum(axis=1) >= min_samples
    return matrix.loc[passes], passes


def tmm_normalize(matrix):
    """
    TMM normalization (Robinson & Oshlack 2010, edgeR method).

    Steps:
      1. Select a reference sample whose 75th-percentile CPM is closest to
         the cross-sample geometric mean of 75th-percentile CPMs.
      2. For each sample, compute log-ratio (M) and mean-log-expression (A)
         vs the reference for genes with counts > 0 in both.
      3. Trim the top/bottom 30% of M values and top/bottom 5% of A values.
      4. Compute a precision-weighted mean of the remaining M values.
         The weight for gene g is 1 / Var(log-ratio) ≈
         (N_s - CPM_s_g) / (N_s * CPM_s_g) + (N_r - CPM_r_g) / (N_r * CPM_r_g),
         where N is the library size.
      5. TMM factor = 2 ^ (weighted mean M).
      6. Effective library size = raw library size * TMM factor.
      7. Output: TMM-CPM = raw counts / effective library size * 1e6.

    Returns normalized CPM matrix and the per-sample TMM factors.
    """
    import numpy as np
    import pandas as pd

    lib_sizes = matrix.sum(axis=0).values.astype(float)
    n_samples = matrix.shape[1]
    samples = matrix.columns.tolist()
    raw = matrix.values.astype(float)

    # Step 1: choose reference sample
    cpm_all = raw / lib_sizes * 1e6
    q75 = np.percentile(cpm_all, 75, axis=0)
    geom_mean_q75 = np.exp(np.mean(np.log(q75 + 1)))
    ref_idx = int(np.argmin(np.abs(q75 - geom_mean_q75)))

    # Step 2–4: compute TMM factor for every sample vs reference
    ref_counts = raw[:, ref_idx]
    ref_lib = lib_sizes[ref_idx]
    tmm_factors = np.ones(n_samples)

    for i in range(n_samples):
        if i == ref_idx:
            continue
        s_counts = raw[:, i]
        s_lib = lib_sizes[i]

        # Only genes expressed in both
        mask = (s_counts > 0) & (ref_counts > 0)
        if mask.sum() < 10:
            continue

        s_cpm = s_counts[mask] / s_lib * 1e6
        r_cpm = ref_counts[mask] / ref_lib * 1e6

        M = np.log2(s_cpm / r_cpm)
        A = 0.5 * np.log2(s_cpm * r_cpm)

        # Trim 30% tails of M, 5% tails of A
        m_lo, m_hi = np.percentile(M, [30, 70])
        a_lo, a_hi = np.percentile(A, [5, 95])
        keep = (M >= m_lo) & (M <= m_hi) & (A >= a_lo) & (A <= a_hi)
        if keep.sum() < 3:
            continue

        M_k = M[keep]
        s_c_k = s_counts[mask][keep]
        r_c_k = ref_counts[mask][keep]

        # Precision weights: inverse variance of log-ratio
        w = 1.0 / (
            (s_lib - s_c_k) / (s_lib * s_c_k) +
            (ref_lib - r_c_k) / (ref_lib * r_c_k)
        )
        tmm_factors[i] = 2 ** (np.sum(w * M_k) / np.sum(w))

    # Step 5–7: effective library size → TMM-CPM
    eff_lib = lib_sizes * tmm_factors
    norm_cpm = raw / eff_lib * 1e6

    norm_df = pd.DataFrame(norm_cpm, index=matrix.index, columns=samples)
    factors_s = pd.Series(tmm_factors, index=samples, name="tmm_factor")
    return norm_df, factors_s, ref_idx


def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    counts_path = os.path.join(base, "data", "counts.tsv")
    metadata_path = os.path.join(base, "data", "sample_metadata.tsv")
    out_path = os.path.join(base, "results", "normalized_counts.tsv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print("Loading data...")
    matrix, gene_info, metadata = load_data(counts_path, metadata_path)
    n_genes_before = matrix.shape[0]
    print(f"  Genes before filtering: {n_genes_before}")
    print(f"  Samples: {matrix.shape[1]}")

    print("\nFiltering (CPM >= 1 in >= 8 samples)...")
    filtered, passes = filter_low_counts(matrix, min_cpm=1.0, min_samples=8)
    n_removed = n_genes_before - filtered.shape[0]
    print(f"  Genes after filtering:  {filtered.shape[0]}")
    print(f"  Genes removed:          {n_removed} ({100*n_removed/n_genes_before:.1f}%)")

    print("\nNormalizing (TMM)...")
    norm_cpm, tmm_factors, ref_idx = tmm_normalize(filtered)
    print(f"  Reference sample: {filtered.columns[ref_idx]}")
    print("\n  TMM factors (should be close to 1.0 for a well-behaved dataset):")
    for s, f in tmm_factors.items():
        cond = metadata.loc[s, "condition"]
        batch = metadata.loc[s, "batch"]
        print(f"    {s}  {cond:8s}  batch {batch}  factor={f:.4f}")

    print("\nWriting results/normalized_counts.tsv ...")
    out = gene_info.loc[filtered.index].join(norm_cpm)
    out.to_csv(out_path, sep="\t")
    print(f"  Shape: {out.shape[0]} genes x {norm_cpm.shape[1]} samples")
    print("Done.")


if __name__ == "__main__":
    main()
