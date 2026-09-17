# Execution Mode: Recurring (Multi-run)
# QC script: library sizes, detected genes, PCA, sample-to-sample correlation.

import os
import sys


def load_data(counts_path, metadata_path):
    import pandas as pd
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
    gene_info = counts[["gene_symbol", "gene_length"]].copy()
    count_matrix = counts.drop(columns=["gene_symbol", "gene_length"])
    return count_matrix, gene_info, metadata


def compute_qc_metrics(count_matrix):
    import pandas as pd
    library_sizes = count_matrix.sum(axis=0)
    detected_genes = (count_matrix > 0).sum(axis=0)
    return pd.DataFrame({
        "library_size": library_sizes,
        "detected_genes": detected_genes,
    })


def normalize_log_cpm(count_matrix):
    import numpy as np
    lib_sizes = count_matrix.sum(axis=0)
    cpm = count_matrix.divide(lib_sizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)
    return log_cpm


def run_pca(log_cpm, n_components=2):
    import numpy as np
    # Center genes (subtract row mean), then run SVD on samples x genes
    X = log_cpm.T.values  # shape: samples x genes
    X_centered = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    coords = U[:, :n_components] * S[:n_components]
    variance_explained = (S ** 2) / (S ** 2).sum() * 100
    return coords, variance_explained


def plot_library_sizes(qc_df, metadata, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    condition_colors = {"control": "#4878CF", "treated": "#D65F5F"}
    batch_markers = {"A": "o", "B": "s"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    samples = qc_df.index.tolist()
    x = np.arange(len(samples))

    for ax, metric, label in [
        (axes[0], "library_size", "Library size (total counts)"),
        (axes[1], "detected_genes", "Detected genes (count > 0)"),
    ]:
        for i, sample in enumerate(samples):
            cond = metadata.loc[sample, "condition"]
            batch = metadata.loc[sample, "batch"]
            ax.bar(i, qc_df.loc[sample, metric],
                   color=condition_colors[cond], alpha=0.85, width=0.7)
            # batch marker overlay at top of bar
            ax.plot(i, qc_df.loc[sample, metric],
                    marker=batch_markers[batch], color="black",
                    markersize=6, zorder=5, linestyle="none")

        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(label)
        ax.set_xlabel("Sample")
        ax.spines[["top", "right"]].set_visible(False)

    # Legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor="#4878CF", label="Control"),
        Patch(facecolor="#D65F5F", label="Treated"),
        Line2D([0], [0], marker="o", color="black", linestyle="none",
               markersize=6, label="Batch A"),
        Line2D([0], [0], marker="s", color="black", linestyle="none",
               markersize=6, label="Batch B"),
    ]
    axes[1].legend(handles=legend_elements, loc="lower right", fontsize=8,
                   frameon=False)

    fig.suptitle("Per-sample QC metrics", fontsize=13, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_pca(coords, variance_explained, metadata, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    condition_colors = {"control": "#4878CF", "treated": "#D65F5F"}
    batch_markers = {"A": "o", "B": "s"}

    fig, ax = plt.subplots(figsize=(6, 5))

    samples = metadata.index.tolist()
    for i, sample in enumerate(samples):
        cond = metadata.loc[sample, "condition"]
        batch = metadata.loc[sample, "batch"]
        ax.scatter(coords[i, 0], coords[i, 1],
                   color=condition_colors[cond],
                   marker=batch_markers[batch],
                   s=80, edgecolors="white", linewidths=0.5, zorder=3)
        ax.annotate(sample, (coords[i, 0], coords[i, 1]),
                    textcoords="offset points", xytext=(5, 3), fontsize=7)

    ax.set_xlabel(f"PC1 ({variance_explained[0]:.1f}% variance)", fontsize=10)
    ax.set_ylabel(f"PC2 ({variance_explained[1]:.1f}% variance)", fontsize=10)
    ax.set_title("PCA of samples (log2 CPM)", fontsize=12)
    ax.axhline(0, color="lightgray", lw=0.7, zorder=1)
    ax.axvline(0, color="lightgray", lw=0.7, zorder=1)
    ax.spines[["top", "right"]].set_visible(False)

    legend_elements = [
        Patch(facecolor="#4878CF", label="Control"),
        Patch(facecolor="#D65F5F", label="Treated"),
        Line2D([0], [0], marker="o", color="gray", linestyle="none",
               markersize=7, label="Batch A"),
        Line2D([0], [0], marker="s", color="gray", linestyle="none",
               markersize=7, label="Batch B"),
    ]
    ax.legend(handles=legend_elements, fontsize=8, frameon=False,
              loc="best")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_correlation(log_cpm, metadata, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    corr = log_cpm.corr(method="pearson")
    samples = corr.index.tolist()
    n = len(samples)

    condition_colors = {"control": "#4878CF", "treated": "#D65F5F"}
    batch_colors = {"A": "#888888", "B": "#CCCCCC"}

    fig, ax = plt.subplots(figsize=(8, 7))

    im = ax.imshow(corr.values, vmin=0.85, vmax=1.0, cmap="Blues", aspect="auto")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(samples, fontsize=8)

    # Annotate cells with correlation value
    for i in range(n):
        for j in range(n):
            val = corr.values[i, j]
            text_color = "white" if val > 0.95 else "black"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                    fontsize=6, color=text_color)

    # Color-coded condition strip on top and left
    for k, sample in enumerate(samples):
        cond = metadata.loc[sample, "condition"]
        color = condition_colors[cond]
        ax.add_patch(plt.Rectangle((k - 0.5, -1.3), 1, 0.8,
                                   color=color, clip_on=False))
        ax.add_patch(plt.Rectangle((-1.3, k - 0.5), 0.8, 1,
                                   color=color, clip_on=False))

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Pearson r (log2 CPM)", fontsize=9)

    ax.set_title("Sample-to-sample correlation (log2 CPM)", fontsize=12, pad=20)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4878CF", label="Control"),
        Patch(facecolor="#D65F5F", label="Treated"),
    ]
    ax.legend(handles=legend_elements, fontsize=8, frameon=False,
              bbox_to_anchor=(1.18, 1.05), loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    counts_path = os.path.join(base, "data", "counts.tsv")
    metadata_path = os.path.join(base, "data", "sample_metadata.tsv")
    figures_dir = os.path.join(base, "results", "figures")
    os.makedirs(figures_dir, exist_ok=True)

    print("Loading data...")
    count_matrix, gene_info, metadata = load_data(counts_path, metadata_path)

    print("\n--- QC metrics ---")
    qc_df = compute_qc_metrics(count_matrix)
    qc_df["condition"] = metadata["condition"]
    qc_df["batch"] = metadata["batch"]
    print(qc_df.to_string())

    print("\n--- Summary by condition ---")
    import pandas as pd
    summary = qc_df.groupby("condition")[["library_size", "detected_genes"]].agg(
        ["mean", "min", "max"]
    )
    print(summary.to_string())

    print("\nNormalizing (log2 CPM)...")
    log_cpm = normalize_log_cpm(count_matrix)

    print("Running PCA...")
    coords, var_exp = run_pca(log_cpm)
    print(f"  PC1: {var_exp[0]:.1f}%  PC2: {var_exp[1]:.1f}%")

    print("Saving figures...")
    plot_library_sizes(
        qc_df, metadata,
        os.path.join(figures_dir, "qc_library_sizes.png")
    )
    plot_pca(
        coords, var_exp, metadata,
        os.path.join(figures_dir, "qc_pca.png")
    )
    plot_correlation(
        log_cpm, metadata,
        os.path.join(figures_dir, "qc_correlation.png")
    )

    print("\nFigures written to results/figures/:")
    print("  qc_library_sizes.png")
    print("  qc_pca.png")
    print("  qc_correlation.png")


if __name__ == "__main__":
    main()
