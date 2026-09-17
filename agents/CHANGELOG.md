# Changelog

All notable changes to this project will be documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com).

## [Unreleased]

### Added
- `environment.yml` — conda environment with FastQC, fastp, Salmon, MultiQC, pytximport, Snakemake
- `config.yaml` — workflow configuration (samples, reference URL, tool params)
- `Snakefile` — seven-rule pipeline: fastqc_raw → fastp → fastqc_trimmed → download_transcriptome → make_tx2gene → salmon_index → salmon_quant → collate → multiqc
- `scripts/make_tx2gene_once_cc.py` — extracts tx2gene map from Ensembl cDNA FASTA headers
- `scripts/collate_counts_once_cc.py` — imports Salmon quant.sf via pytximport, writes genes × samples TSV
- `data/test/SRR6357070_{1,2}.fastq.gz` — nf-core test paired-end sample (GSE110004, *S. cerevisiae*)
- `.gitignore` — added `resources/`, `.snakemake/`, `*.bam`, `data/test/*.fastq.gz`
- `README.md` — added RNA-seq preprocessing workflow section
