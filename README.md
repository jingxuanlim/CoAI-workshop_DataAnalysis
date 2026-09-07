# Gene expression analysis with an AI coding assistant

A 45 minute hands-on session. You take a bulk RNA-seq count matrix from raw
counts to a set of figures using Claude Code, and then you find out how much of
what you produced is actually correct.

---

## Start here

**Before the session:** follow `docs/DataAnalysis_pre-setup.pdf`, then run:

```
python3 selfcheck.py
```

Every line must say PASS.

**During the session:** you need two files, in this order.

1. `docs/01-exercise-analysis.pdf` 
2. `docs/02-exercise-visualization.pdf`

Everything else supports those two.

---

## The data

Sixteen samples, eight control and eight treated, 4200 genes. The counts are
synthetic and were generated from a known model, which means the correct answer
exists and you can measure your result against it. That is the point of the
session and it is not possible with real data.

| File                       | Contents                               |
| -------------------------- | -------------------------------------- |
| `data/counts.tsv`          | Raw count matrix, genes by samples     |
| `data/sample_metadata.tsv` | Condition, replicate, processing batch |
| `data/gene_sets.gmt`       | Gene sets for the enrichment step      |

`data/README.md` describes how the data was made and what is deliberate about it.

---

## What is in this repository

```
data/          the count matrix, sample table and gene sets
docs/          setup, the two exercises, and the printed one-pagers
scripts/       tools you run: enrichment, scoring, data generation
checkpoints/   catch-up scripts, one per stage of the analysis
results/       where your own output goes. Empty to start with.
solutions/     the answer key. Do not open it before Exercise 02.
```

---

## Documents

Useful hand-out.

- `docs/cheatsheet.pdf` — the loop, the commands, three habits
- `docs/prompt-cards.pdf` — twelve prompts to read aloud
- `docs/figure-standards.pdf` — eight checks for any figure
- `docs/verification-checklist.pdf` — what to check before believing an analysis
- `docs/glossary.pdf` — the AI terms and the RNA-seq terms

---

## Running things yourself

```
python3 scripts/make_data.py                              # regenerate the data
python3 checkpoints/01_quality_assessment.py              # catch up to step 2
python3 checkpoints/02_filtered_and_normalized.py         # catch up to step 4
python3 checkpoints/03_differential_expression.py         # catch up to step 6
python3 checkpoints/04_result_figures.py                  # catch up to step 7
python3 scripts/enrich.py results/de_results.tsv          # gene set enrichment
python3 scripts/score_results.py results/de_results.tsv   # compare with the key
```

---

