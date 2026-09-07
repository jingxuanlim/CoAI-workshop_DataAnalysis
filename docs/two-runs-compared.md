# The same prompt, run twice

Exercise 02, Check 1 asks you to compare your analysis with your neighbor's. Use
this instead when there is no time for that, when the room happens to converge,
or when you are reading on your own afterwards.

Two Claude Code sessions were given the identical starting prompt from
`docs/01-exercise-analysis.md`, in identical copies of this repository, one after
the other. Neither could see the other. Both were asked only for a plan.

---

## What both runs agreed on

- 4200 genes, 16 samples, 8 per condition, batch balanced across conditions.
- Filter first, then normalize, then test, then correct.
- Filtering rule: at least 10 counts in at least 8 samples.
- Normalization: median of ratios rather than total count.
- Multiple testing correction: Benjamini-Hochberg, not Bonferroni.
- Both declined to use `gene_length`, and both gave the same correct reason:
  length cancels when comparing one gene across samples.
- Both recommended including the batch term.
- Both stopped and asked a question before writing code.

That is substantial agreement, and it is worth saying plainly. The tool is not
erratic.

---

## Where they diverged

### The statistical test

**Run A** proposed a linear model on log2(normalized count + 0.5), with condition
and batch as predictors. It argued against the negative binomial model on the
grounds that its main advantage, variance shrinkage across genes, is the part
that is awkward to reproduce outside DESeq2, so an unshrunk negative binomial fit
is not obviously better than the simpler model.

**Run B** proposed a negative binomial regression per gene, with the same two
predictors. It argued that the negative binomial is the correct distribution for
count data because it allows the variance to exceed the mean. It then noted the
same shrinkage limitation and said to expect noisier results than DESeq2.

Both arguments are correct. They reach opposite conclusions from the same fact.
These two models do not produce the same gene list.

### The question each one asked

**Run A** asked whether to include the batch term.

**Run B** asked whether the `replicate` column indicates pairing — that is,
whether control replicate 1 and treated replicate 1 come from the same biological
source. If they did, the correct analysis would be a paired test, and every
number would change.

Run A did not raise this. It is a reasonable question about the experimental
design, and only one of the two runs thought to ask it.

### A number

Run A estimated that testing at an uncorrected threshold would give about 200
false positives. Run B said about 150, from a different assumption about how many
genes survive filtering. Neither had run any code. Both said so.

---

## What to take from this

The disagreement is not a defect and it is not a sign that one run was prompted
badly. It is a property of the tool.

Two consequences for your own work:

1. **A result you cannot regenerate is not a result.** Save the analysis as a
   script, put it in version control, and rerun the script. Do not recover a
   number by asking again.
2. **The plan you get is one plausible plan, not the plan.** Run A never
   considered pairing. If you had only run A, that question would not have come
   up, and nothing in its answer would have warned you.

Both runs volunteered that they had not been able to execute any code and that
their numbers were read from the files rather than computed. That is the
behavior you want, and it is worth noticing when you get it.
