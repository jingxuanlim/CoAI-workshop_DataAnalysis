# Figure standards

One page. Print it. Check every figure against it before you accept the figure,
whoever or whatever produced it.

---

### 1. Every axis has a label, and a unit where one exists

"log2 fold change (treated / control)" is a label. "log2FC" is a variable name.

### 2. Any transformation is stated on the axis

If the values are logged, the axis says so and gives the base. A reader should
never have to infer that from the numbers.

### 3. Diverging color scales are centerd at zero

Log2 fold change and z-scores run in both directions from a meaningful center.
The color at zero must be the neutral color. In matplotlib:
`TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)`. A default scale will center
itself on the middle of the data, which puts the neutral color somewhere
arbitrary.

### 4. Colors remain distinguishable to readers with color vision deficiency

About one man in twelve cannot separate red from green. Use blue and orange for
two categories. For continuous scales use viridis, or a blue to red diverging
scale. Never red to green.

### 5. The number of samples per group appears in the figure or its caption

`n = 8 per condition`. Without it, no reader can judge the result.

### 6. Axes are not truncated in a way that exaggerates a difference

A bar chart starts at zero. If an axis must be cut, the figure says so.

### 7. Heatmap rows are scaled, and the scaling is stated

Without scaling, a heatmap of expression shows which genes are highly expressed,
not which genes change. With row scaling it shows the pattern across samples.
These are different figures answering different questions. The color bar label
must say which one you made.

### 8. The figure is produced by a script that is saved and can be rerun

A figure you cannot regenerate is a figure you cannot correct. If it took ten
prompts and no file, you have a picture, not a result.

---

## Two failures that pass all eight checks

**The p-value histogram nobody drew.** A working test on data containing real
differences gives a flat histogram with a spike near zero. A hump in the middle,
a slope, or a spike at one means the test is wrong. This costs one line of code
and is left out of most published analyses.

**The threshold chosen after seeing the data.** Moving the fold change cut-off
from 1.0 to 0.8 because it gives a rounder number of genes is a decision that
never appears in the figure. Fix your thresholds before you look.
