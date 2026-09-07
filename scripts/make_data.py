"""Generate the synthetic RNA-seq dataset used in the workshop.

The script is deterministic. Running it again reproduces the files in data/
byte for byte. Attendees do not need to run it; the generated files are part of
the repository. It is included so that the composition of the dataset can be
inspected and, if needed, changed.

Design summary
--------------
16 samples: 8 control, 8 treated. Two processing batches, balanced across
conditions (4 control and 4 treated in each batch), so batch is not confounded
with condition.

4200 genes. 200 of them are differentially expressed with known effect sizes.
The differentially expressed genes are concentrated in four of the gene sets in
data/gene_sets.gmt, so the enrichment result also has a known answer.

Two deliberate features of the count matrix:

1. A moderate batch effect, visible in a PCA plot. An analysis that ignores it
   still recovers most of the true genes. An analysis that models it recovers
   more.
2. Twelve highly expressed genes, similar in behavior to ribosomal and
   mitochondrial transcripts, that together account for roughly a fifth of the
   counts. Their presence makes the choice of normalization method change the
   result.

Nothing else is added. The long tail of low-count genes comes from the count
model itself and is what a real experiment looks like.

Usage:
    python scripts/make_data.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260818
N_GENES = 4200
N_PER_GROUP = 8
N_DE = 200
N_NAMED_DE = 90        # differential genes drawn from the real-symbol list

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
SOLUTIONS = ROOT / "solutions"

# Real human gene symbols. They are used for the differentially expressed genes
# so that the biological-claim exercise in Block 3 has recognizable material to
# work with. Their expression pattern in this dataset is assigned at random and
# has no relationship to their real biology.
NAMED_GENES = """
TP53 MYC EGFR KRAS BRAF PTEN RB1 CDKN2A CCND1 CDK4 CDK6 MDM2 ATM CHEK2 BRCA1
BRCA2 PALB2 RAD51 XRCC1 ERCC1 MLH1 MSH2 MSH6 PMS2 APC CTNNB1 AXIN2 TCF7L2 WNT5A
FZD7 NOTCH1 NOTCH2 JAG1 DLL4 HES1 HEY1 GLI1 PTCH1 SMO SHH IL6 IL1B IL10 TNF
IFNG CXCL8 CXCL10 CCL2 CCL5 CCR7 CD4 CD8A CD19 CD68 PTPRC FOXP3 GATA3 TBX21
RORC STAT1 STAT3 STAT5A JAK1 JAK2 NFKB1 RELA IKBKB TLR4 MYD88 NLRP3 CASP1
CASP3 CASP8 CASP9 BAX BAK1 BCL2 BCL2L1 MCL1 BID PARP1 XIAP BIRC5 SOD1 SOD2
CAT GPX1 NQO1 HMOX1 NFE2L2 KEAP1 TXN PRDX1 GSR G6PD HK2 PKM LDHA LDHB PDK1
PDHA1 CS IDH1 IDH2 FH SDHA SDHB ACO2 ATP5F1A NDUFA1 NDUFS1 COX5A UQCRC1 SLC2A1
SLC7A11 SLC16A1 HIF1A VEGFA VEGFB ANGPT1 ANGPT2 KDR FLT1 PECAM1 CDH1 CDH2 VIM
SNAI1 SNAI2 TWIST1 ZEB1 ZEB2 MMP2 MMP9 TIMP1 TIMP2 COL1A1 COL3A1 FN1 SPARC
ACTA2 TAGLN TGFB1 TGFBR1 SMAD2 SMAD3 SMAD4 SMAD7 BMP2 BMP4 ID1 ID2 SERPINE1
THBS1 CTGF CYR61 YAP1 TAZ AKT1 PIK3CA MTOR RPS6KB1 EIF4E TSC1 TSC2 PRKAA1
FOXO1 FOXO3 SIRT1 PPARG PPARGC1A CEBPA SREBF1 FASN SCD ACACA CPT1A LPL ADIPOQ
LEP INS IGF1 IGF1R IRS1 GCK GLUT4 HNF4A HNF1A NR5A2 ESR1 PGR AR NR3C1 THRB
"""
NAMED_GENES = NAMED_GENES.split()

# Prefixes used to build the remaining gene symbols. They follow real human gene
# family naming so that the count matrix looks like a real export.
FAMILY_PREFIXES = [
    "ZNF", "OR", "SLC", "TMEM", "KRT", "COL", "ADAM", "MMP", "CYP", "UGT",
    "ABC", "KCN", "CACN", "SCN", "GRIN", "GABR", "HTR", "DRD", "ADRB", "CHRM",
    "RAB", "ARF", "RHO", "GNA", "PLC", "PRKC", "MAPK", "DUSP", "PPP1R", "PPP2R",
    "USP", "UBE2", "RNF", "TRIM", "FBXO", "CUL", "SKP", "ANAPC", "CDC", "CCNB",
    "HIST1H", "H2AF", "SMARC", "ARID", "KDM", "KMT2", "HDAC", "SIRT", "DNMT",
    "TET", "EZH", "SUZ", "BRD", "ATAD", "DDX", "DHX", "EIF", "RPL", "RPS",
    "MRPL", "MRPS", "TIMM", "TOMM", "SEC", "COPB", "VPS", "SNX", "RAB11FIP",
    "MYO", "KIF", "DYNC", "TUBB", "ACTN", "SPTB", "ANK", "TLN", "VCL", "PXN",
]

# The gene sets. Membership is assembled for this exercise and does not
# reproduce any real pathway database. This is stated in data/README.md and is
# the basis of the biological-claim verification exercise.
GENE_SET_NAMES = [
    "OXIDATIVE_PHOSPHORYLATION", "GLYCOLYSIS", "FATTY_ACID_METABOLISM",
    "CHOLESTEROL_HOMEOSTASIS", "DNA_REPAIR", "CELL_CYCLE_G2M_CHECKPOINT",
    "P53_SIGNALING", "APOPTOSIS", "AUTOPHAGY", "UNFOLDED_PROTEIN_RESPONSE",
    "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE",
    "INTERFERON_GAMMA_RESPONSE", "TNFA_SIGNALING_VIA_NFKB",
    "IL6_JAK_STAT3_SIGNALING", "COMPLEMENT_CASCADE", "COAGULATION",
    "ANGIOGENESIS", "HYPOXIA", "EPITHELIAL_MESENCHYMAL_TRANSITION",
    "EXTRACELLULAR_MATRIX_ORGANIZATION", "TGF_BETA_SIGNALING",
    "WNT_BETA_CATENIN_SIGNALING", "NOTCH_SIGNALING", "HEDGEHOG_SIGNALING",
    "PI3K_AKT_MTOR_SIGNALING", "MAPK_SIGNALING", "MYC_TARGETS",
    "E2F_TARGETS", "ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE",
    "ADIPOGENESIS", "MYOGENESIS", "SPERMATOGENESIS", "PEROXISOME",
    "PROTEIN_SECRETION", "MRNA_SPLICING", "RIBOSOME_BIOGENESIS",
    "REACTIVE_OXYGEN_SPECIES", "XENOBIOTIC_METABOLISM",
]

# The four gene sets that the differentially expressed genes are placed in.
TARGET_SETS = [
    "INFLAMMATORY_RESPONSE",
    "TNFA_SIGNALING_VIA_NFKB",
    "HYPOXIA",
    "DNA_REPAIR",
]


def build_gene_table(rng):
    """Return a DataFrame with gene_id, gene_symbol and gene_length."""
    symbols = list(NAMED_GENES)
    seen = set(symbols)
    while len(symbols) < N_GENES:
        prefix = FAMILY_PREFIXES[rng.integers(len(FAMILY_PREFIXES))]
        candidate = f"{prefix}{rng.integers(1, 900)}"
        if candidate not in seen:
            seen.add(candidate)
            symbols.append(candidate)

    order = rng.permutation(N_GENES)
    symbols = [symbols[i] for i in order]

    gene_ids = [f"ENSG{i:011d}" for i in range(1, N_GENES + 1)]
    lengths = np.round(rng.lognormal(mean=7.6, sigma=0.55, size=N_GENES)).astype(int)
    lengths = np.clip(lengths, 250, 60000)

    return pd.DataFrame(
        {"gene_id": gene_ids, "gene_symbol": symbols, "gene_length": lengths}
    )


def build_gene_sets(rng, genes, de_symbols):
    """Assign genes to gene sets.

    The four target sets receive a high proportion of the differentially
    expressed genes. Every other set receives genes drawn at random.
    """
    all_symbols = list(genes["gene_symbol"])
    de_symbols = list(de_symbols)
    rng.shuffle(de_symbols)

    sets = {}
    used_de = 0
    for name in GENE_SET_NAMES:
        size = int(rng.integers(40, 121))
        if name in TARGET_SETS:
            n_de = int(rng.integers(22, 34))
            members = de_symbols[used_de:used_de + n_de]
            used_de += n_de
            pool = [s for s in all_symbols if s not in members]
            filler = list(rng.choice(pool, size=size - len(members), replace=False))
            members = members + filler
        else:
            members = list(rng.choice(all_symbols, size=size, replace=False))
        rng.shuffle(members)
        sets[name] = members
    return sets


def main():
    rng = np.random.default_rng(SEED)
    DATA.mkdir(exist_ok=True)
    SOLUTIONS.mkdir(exist_ok=True)

    genes = build_gene_table(rng)

    # ---- sample table -----------------------------------------------------
    samples = []
    for group_index, condition in enumerate(["control", "treated"]):
        for replicate in range(1, N_PER_GROUP + 1):
            # Batch is balanced within each condition: replicates 1-4 are in
            # batch A, replicates 5-8 are in batch B.
            batch = "A" if replicate <= N_PER_GROUP // 2 else "B"
            samples.append(
                {
                    "sample_id": f"S{group_index * N_PER_GROUP + replicate:02d}",
                    "condition": condition,
                    "replicate": replicate,
                    "batch": batch,
                }
            )
    meta = pd.DataFrame(samples)
    n_samples = len(meta)

    # ---- baseline expression ---------------------------------------------
    baseline = rng.lognormal(mean=3.6, sigma=1.9, size=N_GENES)

    # Twelve genes that dominate the library, in the manner of ribosomal and
    # mitochondrial transcripts.
    high_idx = rng.choice(N_GENES, size=12, replace=False)
    genes.loc[high_idx, "gene_symbol"] = [
        "RPL13A", "RPL10", "RPLP0", "RPS18", "RPS27", "RPL41",
        "MT-CO1", "MT-CO2", "MT-CO3", "MT-ND4", "MT-ATP6", "MT-CYB",
    ]
    baseline[high_idx] = rng.uniform(9000, 26000, size=12)

    # ---- true differential expression ------------------------------------
    # Only genes with a reasonable baseline are eligible, because a gene with
    # almost no counts cannot show a measurable fold change.
    eligible = np.where(baseline > 25)[0]
    eligible = np.setdiff1d(eligible, high_idx)

    # Roughly half of the differential genes are taken from the list of real
    # human gene symbols, and they are given the larger effect sizes. This puts
    # recognizable symbols at the top of the results table, which is what the
    # biological-claim exercise in Block 3 needs. Their behavior here is
    # assigned by this script and has nothing to do with their real biology.
    symbol_array = genes["gene_symbol"].to_numpy()
    is_named = np.isin(symbol_array, NAMED_GENES)
    named_eligible = eligible[is_named[eligible]]
    other_eligible = eligible[~is_named[eligible]]

    n_named_de = min(N_NAMED_DE, len(named_eligible))
    named_de = rng.choice(named_eligible, size=n_named_de, replace=False)
    other_de = rng.choice(other_eligible, size=N_DE - n_named_de, replace=False)

    true_lfc = np.zeros(N_GENES)
    true_lfc[named_de] = rng.uniform(1.5, 3.0, size=n_named_de)
    true_lfc[other_de] = rng.uniform(0.8, 2.0, size=len(other_de))

    # Assign direction. Sixty per cent of the differential genes increase in the
    # treated condition.
    de_idx = np.concatenate([named_de, other_de])
    rng.shuffle(de_idx)
    n_up = int(round(0.6 * N_DE))
    up_idx = de_idx[:n_up]
    down_idx = de_idx[n_up:]
    true_lfc[down_idx] = -true_lfc[down_idx]

    # ---- batch effect -----------------------------------------------------
    # A per-gene shift applied to batch B, plus a small difference in overall
    # sequencing depth between the batches.
    batch_lfc = rng.normal(loc=0.0, scale=0.80, size=N_GENES)
    batch_depth = {"A": 1.00, "B": 0.88}

    # ---- expected counts --------------------------------------------------
    depth_factor = rng.uniform(0.85, 1.15, size=n_samples)

    mu = np.zeros((N_GENES, n_samples))
    for j, row in meta.iterrows():
        effect = true_lfc if row["condition"] == "treated" else np.zeros(N_GENES)
        batch_term = batch_lfc if row["batch"] == "B" else np.zeros(N_GENES)
        scale = depth_factor[j] * batch_depth[row["batch"]]
        mu[:, j] = baseline * np.power(2.0, effect + batch_term) * scale

    # ---- negative binomial sampling ---------------------------------------
    # Dispersion falls as expression rises, which is the relationship seen in
    # real RNA-seq data.
    dispersion = 0.08 + 12.0 / (baseline + 12.0)
    size = 1.0 / dispersion
    counts = np.empty_like(mu, dtype=np.int64)
    for j in range(n_samples):
        p = size / (size + mu[:, j])
        counts[:, j] = rng.negative_binomial(size, p)

    # ---- write files ------------------------------------------------------
    count_df = pd.DataFrame(counts, columns=list(meta["sample_id"]))
    count_df.insert(0, "gene_length", genes["gene_length"].values)
    count_df.insert(0, "gene_symbol", genes["gene_symbol"].values)
    count_df.insert(0, "gene_id", genes["gene_id"].values)
    count_df.to_csv(DATA / "counts.tsv", sep="\t", index=False)
    meta.to_csv(DATA / "sample_metadata.tsv", sep="\t", index=False)

    de_symbols = genes.loc[de_idx, "gene_symbol"]
    gene_sets = build_gene_sets(rng, genes, de_symbols)
    with open(DATA / "gene_sets.gmt", "w") as handle:
        for name, members in gene_sets.items():
            description = "assembled for the workshop; not a real pathway definition"
            handle.write("\t".join([name, description] + members) + "\n")

    # ---- answer key -------------------------------------------------------
    symbol_to_sets = {}
    for name, members in gene_sets.items():
        for member in members:
            symbol_to_sets.setdefault(member, []).append(name)

    key = pd.DataFrame(
        {
            "gene_id": genes["gene_id"],
            "gene_symbol": genes["gene_symbol"],
            "baseline_expression": np.round(baseline, 2),
            "true_log2_fold_change": np.round(true_lfc, 4),
            "is_differential": (true_lfc != 0).astype(int),
            "batch_log2_shift": np.round(batch_lfc, 4),
            "gene_sets": [
                ";".join(symbol_to_sets.get(s, [])) for s in genes["gene_symbol"]
            ],
        }
    )
    key.to_csv(SOLUTIONS / "answer_key.tsv", sep="\t", index=False)

    with open(SOLUTIONS / "answer_key_summary.txt", "w") as handle:
        handle.write("Composition of the synthetic dataset\n")
        handle.write("=" * 40 + "\n\n")
        handle.write(f"Seed: {SEED}\n")
        handle.write(f"Genes: {N_GENES}\n")
        handle.write(f"Samples: {n_samples} ({N_PER_GROUP} per condition)\n")
        handle.write("Batches: A and B, balanced within each condition\n\n")
        handle.write(f"Differentially expressed genes: {N_DE}\n")
        handle.write(f"  increased in treated: {len(up_idx)}\n")
        handle.write(f"  decreased in treated: {len(down_idx)}\n")
        handle.write("  absolute log2 fold change between 0.8 and 3.0\n")
        handle.write(f"  {n_named_de} of them use recognizable real gene symbols\n")
        handle.write("  and were given the larger effect sizes on purpose\n\n")
        handle.write("Gene sets containing an excess of differential genes:\n")
        for name in TARGET_SETS:
            members = set(gene_sets[name])
            n_hit = len(members & set(de_symbols))
            handle.write(f"  {name}: {n_hit} of {len(members)} members\n")
        handle.write("\nBatch effect: per-gene log2 shift, standard deviation 0.80,\n")
        handle.write("applied to batch B, plus a 12 per cent reduction in depth.\n")

    # ---- report -----------------------------------------------------------
    lib = counts.sum(axis=0)
    high_share = counts[high_idx].sum() / counts.sum()
    detected = (counts > 0).sum(axis=0)
    print("Wrote:")
    print(f"  data/counts.tsv            {N_GENES} genes x {n_samples} samples")
    print("  data/sample_metadata.tsv")
    print(f"  data/gene_sets.gmt         {len(gene_sets)} gene sets")
    print("  solutions/answer_key.tsv")
    print("  solutions/answer_key_summary.txt")
    print()
    print(f"Library size range: {lib.min():,} to {lib.max():,}")
    print(f"Detected genes per sample: {detected.min()} to {detected.max()}")
    print(f"Share of counts in the 12 highly expressed genes: {high_share:.1%}")


if __name__ == "__main__":
    main()
