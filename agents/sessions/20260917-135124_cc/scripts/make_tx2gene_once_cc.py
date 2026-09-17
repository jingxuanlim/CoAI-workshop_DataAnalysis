# Execution Mode: One-off
# Extracts a transcript-to-gene mapping from an Ensembl cDNA FASTA.
# Ensembl headers carry the gene ID in a "gene:<id>" field.
# Output: two-column TSV (transcript_id, gene_id), no header.


def extract_tx2gene(fa_path: str, output_path: str) -> None:
    import gzip
    import re
    import pandas as pd

    records = []
    opener = gzip.open if str(fa_path).endswith(".gz") else open
    with opener(fa_path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            tx_id = line[1:].split()[0]
            match = re.search(r"gene:(\S+)", line)
            gene_id = match.group(1) if match else tx_id
            records.append((tx_id, gene_id))
    pd.DataFrame(records).to_csv(output_path, sep="\t", index=False, header=False)


extract_tx2gene(snakemake.input.fa, snakemake.output.tsv)
