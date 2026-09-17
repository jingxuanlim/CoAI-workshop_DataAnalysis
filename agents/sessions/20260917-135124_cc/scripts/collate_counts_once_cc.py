# Execution Mode: One-off
# Imports per-sample Salmon quant.sf files via pytximport, summarizes to
# gene level using a tx2gene map, and writes a genes × samples integer count matrix.


def collate_salmon(quant_dirs, tx2gene_path: str, output_path: str) -> None:
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from pytximport import tximport

    # Normalize: Snakemake passes a str for single-input named fields
    if isinstance(quant_dirs, str):
        quant_dirs = [quant_dirs]
    quant_dirs = list(quant_dirs)

    paths = [Path(d) / "quant.sf" for d in quant_dirs]

    # Fix 2: explicit path→name mapping; never rely on positional ordering
    path_to_name = {str(p): Path(d).name for p, d in zip(paths, quant_dirs)}

    tx2gene = pd.read_csv(
        tx2gene_path,
        sep="\t",
        header=None,
        names=["transcript_id", "gene_id"],
    )

    result = tximport(
        paths,
        data_type="salmon",
        transcript_gene_map=tx2gene,
    )

    x = result.X.toarray() if hasattr(result.X, "toarray") else np.asarray(result.X)

    # Fix 1: Salmon EM produces fractional counts; round to integers for DE tools
    x = np.round(x).astype(int)

    # Map obs_names (full quant.sf paths) back to sample names via explicit dict
    col_names = [path_to_name.get(str(obs), str(obs)) for obs in result.obs_names]

    counts = pd.DataFrame(x, index=col_names, columns=result.var_names).T
    counts.index.name = "gene_id"
    counts.to_csv(output_path, sep="\t")


collate_salmon(
    snakemake.input.quant,
    snakemake.input.tx2gene,
    snakemake.output.tsv,
)
