# checkpoints

Four scripts. Each one takes you to the end of a stage of the analysis, so that
falling behind costs you a minute rather than the rest of the session.

| Script                          | Takes you to  | Writes                          |
| ------------------------------- | ------------- | ------------------------------- |
| `01_quality_assessment.py`      | end of step 2 | three figures                   |
| `02_filtered_and_normalized.py` | end of step 4 | `results/normalized_counts.tsv` |
| `03_differential_expression.py` | end of step 6 | `results/de_results.tsv`        |
| `04_result_figures.py`          | end of step 7 | four figures                    |

Run them from the top of the repository:

```
python3 checkpoints/03_differential_expression.py
```

Scripts 3 and 4 run the earlier ones automatically if their input is missing, so
you can jump straight to the one you need.

Read the script you run. They are short, they are commented, and they are one
version of what the exercise is asking you to produce. They are not the only
correct answer, and the code Claude writes for you will not look like this.
