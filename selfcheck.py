"""Check that this computer is ready for the workshop.

Run it from the top level of the repository:

    python selfcheck.py

It prints one line per check and a verdict at the end. Run it once when you
install everything, and again on the morning of the workshop.
"""

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

REQUIRED_PACKAGES = [
    ("numpy", "1.24"),
    ("pandas", "2.0"),
    ("scipy", "1.10"),
    ("matplotlib", "3.7"),
    ("statsmodels", "0.14"),
]

REQUIRED_FILES = [
    "data/counts.tsv",
    "data/sample_metadata.tsv",
    "data/gene_sets.gmt",
    "scripts/enrich.py",
    "scripts/score_results.py",
    "docs/01-exercise-analysis.pdf",
]

results = []


def check(name, passed, detail=""):
    mark = "PASS" if passed else "FAIL"
    print(f"[{mark}] {name}" + (f"  {detail}" if detail else ""))
    results.append((name, passed, detail))
    return passed


def version_at_least(found, wanted):
    def parts(v):
        out = []
        for piece in v.split("."):
            digits = "".join(c for c in piece if c.isdigit())
            out.append(int(digits) if digits else 0)
        return out
    a, b = parts(found), parts(wanted)
    a += [0] * (len(b) - len(a))
    b += [0] * (len(a) - len(b))
    return a >= b


print("Workshop self-check")
print("=" * 58)

# 1. Python version.
version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
check("Python 3.9 or later", sys.version_info >= (3, 9), f"found {version}")

# 2. Packages.
for package, minimum in REQUIRED_PACKAGES:
    try:
        module = importlib.import_module(package)
        found = getattr(module, "__version__", "unknown")
        ok = found == "unknown" or version_at_least(found, minimum)
        check(f"{package} {minimum} or later", ok, f"found {found}")
    except ImportError:
        check(f"{package} {minimum} or later", False, "not installed")

# 3. Matplotlib can draw and save a figure.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots()
    axes.plot([0, 1], [0, 1])
    target = ROOT / ".selfcheck_figure.png"
    figure.savefig(target)
    plt.close(figure)
    target.unlink()
    check("matplotlib can write a figure file", True)
except Exception as error:                                   # noqa: BLE001
    check("matplotlib can write a figure file", False, str(error)[:60])

# 4. Claude Code.
claude = shutil.which("claude")
if claude:
    try:
        proc = subprocess.run([claude, "--version"], capture_output=True,
                              text=True, timeout=30)
        check("Claude Code is installed", proc.returncode == 0,
              proc.stdout.strip()[:40])
    except Exception:                                        # noqa: BLE001
        check("Claude Code is installed", False, "found but would not start")
else:
    check("Claude Code is installed", False, "'claude' not found on PATH")

# 5. Repository contents.
missing = [f for f in REQUIRED_FILES if not (ROOT / f).exists()]
check("Workshop files are present", not missing,
      "missing: " + ", ".join(missing) if missing else "")

# 6. The data loads and has the expected shape.
try:
    import pandas as pd
    counts = pd.read_csv(ROOT / "data" / "counts.tsv", sep="\t", nrows=5)
    meta = pd.read_csv(ROOT / "data" / "sample_metadata.tsv", sep="\t")
    shape_ok = len(meta) == 16 and len(counts.columns) == 19
    check("The count matrix loads and has 16 samples", shape_ok,
          f"{len(meta)} samples in the metadata table")
except Exception as error:                                   # noqa: BLE001
    check("The count matrix loads and has 16 samples", False, str(error)[:60])

failed = [name for name, ok, _ in results if not ok]

print("=" * 58)
if not failed:
    print("Ready. Nothing further to do before the workshop.")
    print()
    print("Authentication is not checked here. Start Claude Code once with the")
    print("command 'claude' inside this folder and confirm that it responds.")
    sys.exit(0)

print(f"{len(failed)} check(s) failed:")
for name in failed:
    print(f"  - {name}")
print()
print("Fixes are in docs/00-setup.md, section 'If a check fails'.")
print("If you cannot resolve it, come to the session anyway. You can pair with")
print("another attendee and every block has material that works without Claude.")
sys.exit(1)
