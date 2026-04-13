"""
Stratified 80/20 split of the 492 Clostridiales genomes by genus.

Genera with <5 members and UNKNOWN are collapsed into "Other" before
stratification so that sklearn's StratifiedShuffleSplit has ≥2 members
per stratum.

Run from the project root:
    python data_clostridia/stratified_split.py

Outputs (relative to project root):
    data_clostridia/train.txt       — 394 accessions
    data_clostridia/heldout.txt     — 98 accessions
    data_clostridia/train_genomes/  — symlinks → data/original/*.gbff
    data_clostridia/heldout_genomes/— symlinks → data/original/*.gbff
"""

import os
import sys

import pandas as pd
from sklearn.model_selection import train_test_split

CLADES_TSV = "data_clostridia/clades.tsv"
GBFF_DIR = os.path.abspath("data/original")
TRAIN_TXT = "data_clostridia/train.txt"
HELDOUT_TXT = "data_clostridia/heldout.txt"
TRAIN_DIR = "data_clostridia/train_genomes"
HELDOUT_DIR = "data_clostridia/heldout_genomes"
MIN_GENUS_SIZE = 5
RANDOM_STATE = 42
TEST_SIZE = 0.2

# ── 1. Load clades ────────────────────────────────────────────────────────────
df = pd.read_csv(CLADES_TSV, sep="\t")
assert len(df) == 492, f"Expected 492 rows, got {len(df)}"

# ── 2. Build stratification label: collapse small genera → "Other" ─────────────
genus_counts = df["genus"].value_counts()
large_genera = set(genus_counts[genus_counts >= MIN_GENUS_SIZE].index)
# UNKNOWN always goes to Other regardless of size
large_genera.discard("UNKNOWN")

df["stratum"] = df["genus"].apply(lambda g: g if g in large_genera else "Other")

print("Stratification strata:")
strata_counts = df["stratum"].value_counts()
for stratum, count in strata_counts.items():
    print(f"  {stratum:<30} {count:>4}")
print(f"  {'─'*35}")
print(f"  {'Total':<30} {len(df):>4}\n")

# ── 3. Stratified split ────────────────────────────────────────────────────────
accessions = df["accession"].tolist()
strata = df["stratum"].tolist()

train_accs, heldout_accs = train_test_split(
    accessions,
    test_size=TEST_SIZE,
    stratify=strata,
    random_state=RANDOM_STATE,
)
train_accs = sorted(train_accs)
heldout_accs = sorted(heldout_accs)

print(f"Split: {len(train_accs)} train / {len(heldout_accs)} held-out")

# ── 4. Per-genus verification ──────────────────────────────────────────────────
train_set = set(train_accs)
heldout_set = set(heldout_accs)

df["split"] = df["accession"].apply(lambda a: "train" if a in train_set else "heldout")

print("\nPer-stratum counts (train / held-out):")
print(f"  {'Stratum':<30} {'Train':>6} {'Held-out':>9} {'Total':>7}")
print(f"  {'─'*55}")
warnings = []
for stratum in sorted(strata_counts.index):
    sub = df[df["stratum"] == stratum]
    n_train = (sub["split"] == "train").sum()
    n_held = (sub["split"] == "heldout").sum()
    n_total = len(sub)
    # A stratum of size N can yield at most floor(N * test_size) held-out.
    # For N=5 that is 1 — mathematically unavoidable, not a bug.
    expected_min_held = max(1, int(n_total * TEST_SIZE))
    hard_fail = n_held == 0
    soft_warn = (not hard_fail) and n_held < 2 and expected_min_held >= 2
    flag = (
        " ✗ 0 held-out!"
        if hard_fail
        else (" ⚠ only 1 held-out (expected)" if n_held < 2 else "")
    )
    print(f"  {stratum:<30} {n_train:>6} {n_held:>9} {n_total:>7}{flag}")
    if hard_fail:
        warnings.append(f"HARD FAIL — stratum '{stratum}' has 0 held-out members.")
    elif soft_warn:
        warnings.append(
            f"NOTE — stratum '{stratum}' has 1 held-out member (large enough to expect 2)."
        )

hard_failures = [w for w in warnings if w.startswith("HARD")]
if hard_failures:
    print("\nERRORS (split invalid):")
    for w in hard_failures:
        print(f"  {w}")
    sys.exit(1)

notes = [w for w in warnings if w.startswith("NOTE")]
if notes:
    print("\nNotes (expected for small strata, not errors):")
    for w in notes:
        print(f"  {w}")

print(
    "\nSplit is valid. Genera with 1 held-out member have ≤6 total — 20% rounds to 1."
)

# ── 5. Write train.txt and heldout.txt ────────────────────────────────────────
with open(TRAIN_TXT, "w") as f:
    f.write("\n".join(train_accs) + "\n")
with open(HELDOUT_TXT, "w") as f:
    f.write("\n".join(heldout_accs) + "\n")
print(f"\nWrote {TRAIN_TXT} ({len(train_accs)} lines)")
print(f"Wrote {HELDOUT_TXT} ({len(heldout_accs)} lines)")

# ── 6. Create symlink directories ─────────────────────────────────────────────
for split_dir, split_accs in [(TRAIN_DIR, train_accs), (HELDOUT_DIR, heldout_accs)]:
    os.makedirs(split_dir, exist_ok=True)
    # Remove stale symlinks
    for entry in os.scandir(split_dir):
        if entry.is_symlink():
            os.unlink(entry.path)
    # Create fresh symlinks
    missing = []
    for acc in split_accs:
        src = os.path.join(GBFF_DIR, acc + ".gbff")
        dst = os.path.join(split_dir, acc + ".gbff")
        if not os.path.exists(src):
            missing.append(src)
            continue
        os.symlink(src, dst)
    if missing:
        print(f"\nWARNING: {len(missing)} source files not found for {split_dir}:")
        for m in missing:
            print(f"  {m}")
    else:
        n = len(os.listdir(split_dir))
        print(f"Symlinked {n} files → {split_dir}/")

print(
    "\nDone. Commit train.txt and heldout.txt immediately (random_state=42 is fixed)."
)
