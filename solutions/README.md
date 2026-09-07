# solutions

**Attendees: do not open this folder before Block 3.**

Knowing which genes really differ before you run the analysis changes what you
do, usually without your noticing. Block 3 has you compare against the key
through `scripts/score_results.py`, which is the intended route and reports the
comparison without showing you the gene list.

## Contents

| File                     | What it is                                                           |
| ------------------------ | -------------------------------------------------------------------- |
| `answer_key.tsv`         | Every gene: true log2 fold change, batch shift, gene set membership. |
| `answer_key_summary.txt` | The same thing in one screen.                                        |
| `reference_outputs/`     | Reference results, enrichment, benchmark and figures.                |

## What was built into the dataset

Do not read this before Block 3.

**A batch effect.** Samples processed in batch B differ from batch A for reasons
unrelated to the treatment. It is the largest single source of variance in the
data. In the PCA figure, PC1 separates the batches and PC2 separates the
conditions, which is the opposite of what most attendees expect. An analysis
that ignores batch still works. One that includes it in the model recovers about
20 more true genes.

**Twelve dominant genes.** Twelve genes carry about a fifth of all counts, as
ribosomal and mitochondrial transcripts do in a real experiment. Because of them,
normalizing by total count and normalizing by median of ratios give different
answers. This is what makes "why did you choose that normalization method?" a
question with a real answer.

Nothing else was added. The long tail of genes with very low counts comes from
the negative binomial count model itself.
