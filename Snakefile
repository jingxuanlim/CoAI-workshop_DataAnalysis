# No default configfile: here — always pass --configfile explicitly so configs
# are never silently merged (Snakemake merges configfile: with --configfile).
# Usage: snakemake --cores N --configfile config.yaml
#        snakemake --cores N --configfile config_zenodo.yaml

# ── Config ────────────────────────────────────────────────────────────────────
SAMPLES     = config["samples"]   # {name: {layout, fq1, fq2?}}
PE_SAMPLES  = [s for s, v in SAMPLES.items() if v["layout"] == "PE"]
SE_SAMPLES  = [s for s, v in SAMPLES.items() if v["layout"] == "SE"]
ALL_SAMPLES = list(SAMPLES.keys())

RESULTS  = config["results_dir"]
REF_DIR  = f"{config['resources_dir']}/{config['reference_version']}"

# Regex constraints keep PE/SE rules from ambiguously matching each other.
# "__nomatch__" is a sentinel that deliberately matches nothing real.
PE_RE = "|".join(PE_SAMPLES) if PE_SAMPLES else "__nomatch__"
SE_RE = "|".join(SE_SAMPLES) if SE_SAMPLES else "__nomatch__"

# ── Input helpers ─────────────────────────────────────────────────────────────
def fq1(wc): return SAMPLES[wc.sample]["fq1"]
def fq2(wc): return SAMPLES[wc.sample]["fq2"]

# ── Target lists (used by rule all and multiqc) ───────────────────────────────
FASTQC_RAW = (
    expand(f"{RESULTS}/fastqc/raw/{{s}}_1_fastqc.zip", s=PE_SAMPLES) +
    expand(f"{RESULTS}/fastqc/raw/{{s}}_2_fastqc.zip", s=PE_SAMPLES) +
    expand(f"{RESULTS}/fastqc/raw/{{s}}_fastqc.zip",   s=SE_SAMPLES)
)
FASTQC_TRIM = (
    expand(f"{RESULTS}/fastqc/trimmed/{{s}}_1_fastqc.zip", s=PE_SAMPLES) +
    expand(f"{RESULTS}/fastqc/trimmed/{{s}}_2_fastqc.zip", s=PE_SAMPLES) +
    expand(f"{RESULTS}/fastqc/trimmed/{{s}}_fastqc.zip",   s=SE_SAMPLES)
)
FASTP_JSONS  = expand(f"{RESULTS}/trimmed/{{s}}_fastp.json", s=ALL_SAMPLES)
SALMON_DIRS  = expand(f"{RESULTS}/salmon/{{s}}",             s=ALL_SAMPLES)


rule all:
    input:
        f"{RESULTS}/counts_matrix.tsv",
        f"{RESULTS}/multiqc/{config['multiqc']['report_name']}",


# ── Reference ─────────────────────────────────────────────────────────────────

rule download_transcriptome:
    output: fa=f"{REF_DIR}/transcriptome.fa.gz"
    params: url=config["reference"]["transcriptome_url"]
    shell: "curl -fsSL -o {output.fa} {params.url}"


# Fix 7: verify download integrity; populate reference.transcriptome_sha256 in
# config after the first successful download using:
#   shasum -a 256 resources/.../transcriptome.fa.gz
rule verify_transcriptome:
    input:  fa=f"{REF_DIR}/transcriptome.fa.gz"
    output: flag=touch(f"{REF_DIR}/.transcriptome_verified")
    params: expected=config["reference"].get("transcriptome_sha256", "")
    shell:
        """
        if [ -n "{params.expected}" ]; then
            echo "{params.expected}  {input.fa}" | shasum -a 256 -c -
        else
            echo "Warning: no transcriptome_sha256 in config; skipping verification"
        fi
        """


rule download_genome:
    output: fa=f"{REF_DIR}/genome.fa.gz"
    params: url=config["reference"]["genome_url"]
    shell: "curl -fsSL -o {output.fa} {params.url}"


# Fix 5: decoy-aware index — genome seqs act as decoys to reduce spurious mapping
rule make_decoys:
    input:  fa=f"{REF_DIR}/genome.fa.gz"
    output: txt=f"{REF_DIR}/decoys.txt"
    shell:  "gunzip -c {input.fa} | grep '^>' | cut -d ' ' -f 1 | sed 's/>//' > {output.txt}"


rule make_gentrome:
    input:
        tx=f"{REF_DIR}/transcriptome.fa.gz",
        genome=f"{REF_DIR}/genome.fa.gz",
    output: fa=f"{REF_DIR}/gentrome.fa.gz"
    shell:  "cat {input.tx} {input.genome} > {output.fa}"


rule make_tx2gene:
    input:  fa=f"{REF_DIR}/transcriptome.fa.gz"
    output: tsv=f"{REF_DIR}/tx2gene.tsv"
    script: "agents/sessions/20260917-135124_cc/scripts/make_tx2gene_once_cc.py"


rule salmon_index:
    input:
        gentrome=f"{REF_DIR}/gentrome.fa.gz",
        decoys=f"{REF_DIR}/decoys.txt",
        verified=f"{REF_DIR}/.transcriptome_verified",
    output: directory(f"{REF_DIR}/salmon_index")
    threads: config["salmon"]["threads"]
    benchmark: f"{RESULTS}/benchmarks/salmon_index.txt"
    shell:
        "salmon index -t {input.gentrome} -d {input.decoys} -i {output} --threads {threads}"


# ── QC (raw) ──────────────────────────────────────────────────────────────────
# FastQC names outputs from the input filename, so we run into a temp dir and
# rename. For PE the two output files are sorted alphabetically by filename,
# which matches R1/R2 order for all standard naming conventions (_1/_2, _R1/_R2).

rule fastqc_raw_pe:
    input: r1=fq1, r2=fq2
    output:
        html1=f"{RESULTS}/fastqc/raw/{{sample}}_1_fastqc.html",
        zip1=f"{RESULTS}/fastqc/raw/{{sample}}_1_fastqc.zip",
        html2=f"{RESULTS}/fastqc/raw/{{sample}}_2_fastqc.html",
        zip2=f"{RESULTS}/fastqc/raw/{{sample}}_2_fastqc.zip",
    wildcard_constraints: sample=PE_RE
    threads: config["fastqc"]["threads"]
    benchmark: f"{RESULTS}/benchmarks/{{sample}}_fastqc_raw.txt"
    shell:
        """
        tmpdir=$(mktemp -d)
        fastqc {input.r1} -o $tmpdir -t {threads}
        cp $tmpdir/*_fastqc.html {output.html1}
        cp $tmpdir/*_fastqc.zip  {output.zip1}
        rm $tmpdir/*
        fastqc {input.r2} -o $tmpdir -t {threads}
        cp $tmpdir/*_fastqc.html {output.html2}
        cp $tmpdir/*_fastqc.zip  {output.zip2}
        rm -rf $tmpdir
        """


rule fastqc_raw_se:
    input: r1=fq1
    output:
        html=f"{RESULTS}/fastqc/raw/{{sample}}_fastqc.html",
        zip=f"{RESULTS}/fastqc/raw/{{sample}}_fastqc.zip",
    wildcard_constraints: sample=SE_RE
    threads: config["fastqc"]["threads"]
    benchmark: f"{RESULTS}/benchmarks/{{sample}}_fastqc_raw.txt"
    shell:
        """
        tmpdir=$(mktemp -d)
        fastqc {input.r1} -o $tmpdir -t {threads}
        cp $tmpdir/*_fastqc.html {output.html}
        cp $tmpdir/*_fastqc.zip  {output.zip}
        rm -rf $tmpdir
        """


# ── Trimming ──────────────────────────────────────────────────────────────────

rule fastp_pe:
    input: r1=fq1, r2=fq2
    output:
        r1=f"{RESULTS}/trimmed/{{sample}}_1.fastq.gz",
        r2=f"{RESULTS}/trimmed/{{sample}}_2.fastq.gz",
        json=f"{RESULTS}/trimmed/{{sample}}_fastp.json",
        html=f"{RESULTS}/trimmed/{{sample}}_fastp.html",
    wildcard_constraints: sample=PE_RE
    threads: config["fastp"]["threads"]
    params: min_len=config["fastp"]["min_length"]
    benchmark: f"{RESULTS}/benchmarks/{{sample}}_fastp.txt"
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-o {output.r1} -O {output.r2} "
        "--json {output.json} --html {output.html} "
        "--length_required {params.min_len} "
        "--thread {threads}"


rule fastp_se:
    input: r1=fq1
    output:
        r1=f"{RESULTS}/trimmed/{{sample}}.fastq.gz",
        json=f"{RESULTS}/trimmed/{{sample}}_fastp.json",
        html=f"{RESULTS}/trimmed/{{sample}}_fastp.html",
    wildcard_constraints: sample=SE_RE
    threads: config["fastp"]["threads"]
    params: min_len=config["fastp"]["min_length"]
    benchmark: f"{RESULTS}/benchmarks/{{sample}}_fastp.txt"
    shell:
        "fastp -i {input.r1} "
        "-o {output.r1} "
        "--json {output.json} --html {output.html} "
        "--length_required {params.min_len} "
        "--thread {threads}"


# ── QC (trimmed) ──────────────────────────────────────────────────────────────

rule fastqc_trimmed_pe:
    input:
        r1=f"{RESULTS}/trimmed/{{sample}}_1.fastq.gz",
        r2=f"{RESULTS}/trimmed/{{sample}}_2.fastq.gz",
    output:
        html1=f"{RESULTS}/fastqc/trimmed/{{sample}}_1_fastqc.html",
        zip1=f"{RESULTS}/fastqc/trimmed/{{sample}}_1_fastqc.zip",
        html2=f"{RESULTS}/fastqc/trimmed/{{sample}}_2_fastqc.html",
        zip2=f"{RESULTS}/fastqc/trimmed/{{sample}}_2_fastqc.zip",
    wildcard_constraints: sample=PE_RE
    threads: config["fastqc"]["threads"]
    benchmark: f"{RESULTS}/benchmarks/{{sample}}_fastqc_trimmed.txt"
    shell:
        """
        tmpdir=$(mktemp -d)
        fastqc {input.r1} -o $tmpdir -t {threads}
        cp $tmpdir/*_fastqc.html {output.html1}
        cp $tmpdir/*_fastqc.zip  {output.zip1}
        rm $tmpdir/*
        fastqc {input.r2} -o $tmpdir -t {threads}
        cp $tmpdir/*_fastqc.html {output.html2}
        cp $tmpdir/*_fastqc.zip  {output.zip2}
        rm -rf $tmpdir
        """


rule fastqc_trimmed_se:
    input: r1=f"{RESULTS}/trimmed/{{sample}}.fastq.gz"
    output:
        html=f"{RESULTS}/fastqc/trimmed/{{sample}}_fastqc.html",
        zip=f"{RESULTS}/fastqc/trimmed/{{sample}}_fastqc.zip",
    wildcard_constraints: sample=SE_RE
    threads: config["fastqc"]["threads"]
    benchmark: f"{RESULTS}/benchmarks/{{sample}}_fastqc_trimmed.txt"
    shell:
        """
        tmpdir=$(mktemp -d)
        fastqc {input.r1} -o $tmpdir -t {threads}
        cp $tmpdir/*_fastqc.html {output.html}
        cp $tmpdir/*_fastqc.zip  {output.zip}
        rm -rf $tmpdir
        """


# ── Quantification ────────────────────────────────────────────────────────────

rule salmon_quant_pe:
    input:
        r1=f"{RESULTS}/trimmed/{{sample}}_1.fastq.gz",
        r2=f"{RESULTS}/trimmed/{{sample}}_2.fastq.gz",
        index=f"{REF_DIR}/salmon_index",
    output: directory(f"{RESULTS}/salmon/{{sample}}")
    wildcard_constraints: sample=PE_RE
    threads: config["salmon"]["threads"]
    benchmark: f"{RESULTS}/benchmarks/{{sample}}_salmon.txt"
    shell:
        "salmon quant -i {input.index} -l A "
        "-1 {input.r1} -2 {input.r2} "
        "--validateMappings -p {threads} -o {output}"


rule salmon_quant_se:
    input:
        r1=f"{RESULTS}/trimmed/{{sample}}.fastq.gz",
        index=f"{REF_DIR}/salmon_index",
    output: directory(f"{RESULTS}/salmon/{{sample}}")
    wildcard_constraints: sample=SE_RE
    threads: config["salmon"]["threads"]
    benchmark: f"{RESULTS}/benchmarks/{{sample}}_salmon.txt"
    shell:
        "salmon quant -i {input.index} -l A "
        "-r {input.r1} "
        "--validateMappings -p {threads} -o {output}"


# ── Count matrix ──────────────────────────────────────────────────────────────

rule collate:
    input:
        quant=SALMON_DIRS,
        tx2gene=f"{REF_DIR}/tx2gene.tsv",
    output: tsv=f"{RESULTS}/counts_matrix.tsv"
    script: "agents/sessions/20260917-135124_cc/scripts/collate_counts_once_cc.py"


# ── MultiQC ───────────────────────────────────────────────────────────────────

rule multiqc:
    input: FASTQC_RAW, FASTQC_TRIM, FASTP_JSONS, SALMON_DIRS
    output: f"{RESULTS}/multiqc/{config['multiqc']['report_name']}"
    params:
        outdir=f"{RESULTS}/multiqc",
        name=config["multiqc"]["report_name"],
    benchmark: f"{RESULTS}/benchmarks/multiqc.txt"
    shell:
        "multiqc {RESULTS}/fastqc/ {RESULTS}/trimmed/ {RESULTS}/salmon/ "
        "-o {params.outdir} -n {params.name} --force"
