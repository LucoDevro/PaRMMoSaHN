"""
Build clades.tsv — master taxonomy database for the 492 Clostridiales genomes.

Columns: accession, filename, genus, species, strain, matched_in_excel

Run from the project root:
    python data_clostridia/build_clades.py
"""

import os
import re

import pandas as pd

GBFF_DIR = "data/original"
EXCEL = "additional_file1_genomes_metadata.xlsx"
OUT = "data_clostridia/clades.tsv"

# ── 1. Collect all accessions from .gbff filenames ───────────────────────────
gbff_files = sorted(f for f in os.listdir(GBFF_DIR) if f.endswith(".gbff"))
accessions = [f.replace(".gbff", "") for f in gbff_files]
print(f"GBFF files found: {len(gbff_files)}")

# ── 2. Load metadata Excel ────────────────────────────────────────────────────
df_meta = pd.read_excel(
    EXCEL,
    usecols=[
        "Assembly Accession",
        "Organism Name",
        "Organism Taxonomic ID",
        "Organism Infraspecific Names Strain",
    ],
)
df_meta.columns = ["accession", "organism_name", "tax_id", "strain"]
df_meta = df_meta.dropna(subset=["accession"])

# ── 3. Match our 492 against the Excel ───────────────────────────────────────
our_set = set(accessions)
matched = df_meta[df_meta["accession"].isin(our_set)].copy()
unmatched_accs = sorted(our_set - set(df_meta["accession"]))
print(f"Matched in Excel:   {len(matched)}")
print(f"Unmatched (no metadata): {len(unmatched_accs)}")
if unmatched_accs:
    print("  Unmatched accessions:")
    for a in unmatched_accs:
        print(f"    {a}")


# ── 4. Parse genus and species from Organism Name ─────────────────────────────
def parse_organism(name):
    """Return (genus, species) from an NCBI Organism Name string."""
    if not isinstance(name, str):
        return "UNKNOWN", "UNKNOWN"
    name = name.strip()
    # Strip leading/trailing brackets (reclassified genera, e.g. [Clostridium])
    name_clean = re.sub(r"^\[|\]", "", name)
    parts = name_clean.split()
    genus = parts[0] if parts else "UNKNOWN"
    species = (
        " ".join(parts[:2]) if len(parts) >= 2 else parts[0] if parts else "UNKNOWN"
    )
    return genus, species


matched[["genus", "species"]] = matched["organism_name"].apply(
    lambda x: pd.Series(parse_organism(x))
)

# ── 5. Build rows for unmatched accessions ────────────────────────────────────
unmatched_rows = pd.DataFrame(
    {
        "accession": unmatched_accs,
        "organism_name": "UNKNOWN",
        "tax_id": None,
        "strain": None,
        "genus": "UNKNOWN",
        "species": "UNKNOWN",
    }
)

# ── 6. Combine and add filename / matched flag ────────────────────────────────
df = pd.concat([matched, unmatched_rows], ignore_index=True)
df["filename"] = df["accession"] + ".gbff"
df["matched_in_excel"] = df["organism_name"] != "UNKNOWN"

# Reorder columns
df = df[
    [
        "accession",
        "filename",
        "genus",
        "species",
        "strain",
        "tax_id",
        "organism_name",
        "matched_in_excel",
    ]
]
df = df.sort_values("accession").reset_index(drop=True)

# ── 7. Write clades.tsv ───────────────────────────────────────────────────────
df.to_csv(OUT, sep="\t", index=False)
print(f"\nWrote {len(df)} rows → {OUT}")

# ── 8. Report genus distribution ──────────────────────────────────────────────
genus_counts = df["genus"].value_counts()
print("\nGenus distribution (all 492 genomes):")
print(f"  {'Genus':<30} {'Count':>5}")
print(f"  {'-'*30} {'-----':>5}")
for genus, count in genus_counts.items():
    print(f"  {genus:<30} {count:>5}")

# Warn about stratification feasibility
small_genera = genus_counts[genus_counts < 5]
print(
    f"\nGenera with <5 members (will be collapsed to 'Other' for stratification): {len(small_genera)}"
)
print(
    f"  {', '.join(small_genera.index.tolist()[:20])}"
    + (" ..." if len(small_genera) > 20 else "")
)
print(f"\nGenera with ≥5 members (stratifiable): {(genus_counts >= 5).sum()}")
print(f"  {', '.join(genus_counts[genus_counts >= 5].index.tolist())}")
