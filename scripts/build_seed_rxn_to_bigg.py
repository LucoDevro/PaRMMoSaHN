#!/usr/bin/env python3
"""Build src/parmmosahn/data/seed_rxn_to_bigg.tsv from MetaNetX reac_xref.tsv.

Chains SEED reaction ID -> MNXref MNXR ID -> BiGG reaction ID using the
official MetaNetX cross-reference table for reactions. Sister of
build_seed_to_bigg.py (which does the same for metabolites).

Usage
-----
    python scripts/build_seed_rxn_to_bigg.py \\
        --reac-xref ~/thesis/data/external/reac_xref.tsv \\
        --output    src/parmmosahn/data/seed_rxn_to_bigg.tsv

reac_xref.tsv is downloaded from
https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv (current is
MNXref 4.5, CC-BY 4.0, dated 2025-08-13). The file is ~78 MB and is NOT
bundled in the wheel.

The output TSV is two-column (seed_rxn_id, bigg_rxn_id) with a header.
Idempotent: rerunning on the same reac_xref.tsv produces a byte-identical
output.

Notes on chain semantics
------------------------
Reactions are more contested across namespaces than metabolites because
the same MNX reaction can correspond to BiGG reactions that differ only in
direction or compartment annotation (e.g. GLCpts vs GLCpts_periplasm vs
GLCpts_external). We follow the same first-wins convention as the
metabolite builder: the first bigg.reaction line in source order is taken
as canonical. This matches the convention used elsewhere in the package
and produces reproducible output.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

DEFAULT_REAC_XREF = Path.home() / "thesis/data/external/reac_xref.tsv"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parent.parent / "src/parmmosahn/data/seed_rxn_to_bigg.tsv"
)

# Manual reaction-side overrides. Empty for the first cut — metabolite
# stereoisomers don't have a clean reaction analogue, and direction /
# compartment ambiguity is better handled by the first-wins MNX chain
# than by hand-curated overrides. Reserved here so the script's shape
# matches the metabolite builder and future supervisor-verified entries
# have a documented home.
BIGG_RXN_OVERRIDES: dict[str, str] = {}


def parse_reac_xref(path: Path) -> tuple[dict[str, str], dict[str, str], str | None]:
    """Stream reac_xref.tsv once and build the SEED->MNXR and MNXR->BiGG dicts.

    Returns
    -------
    seed_to_mnxr : {rxn00001: MNXR1}
    mnxr_to_bigg : {MNXR1: GLCpts}
    version      : MNXref version string from the file header, or None
    """
    seed_to_mnxr: dict[str, str] = {}
    mnxr_to_bigg: dict[str, str] = {}
    version: str | None = None

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#VERSION:") and version is None:
                version = line.split(":", 1)[1].strip()
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            source, mnxr_id = parts[0], parts[1]
            if source.startswith("seed.reaction:"):
                seed_to_mnxr[source.split(":", 1)[1]] = mnxr_id
            elif source.startswith("bigg.reaction:"):
                if mnxr_id not in mnxr_to_bigg:
                    bigg = source.split(":", 1)[1]
                    # Strip the SBML SId-safe "R_" prefix retained by some
                    # MetaNetX-imported BiGG entries (~11 % of the table).
                    # Verified empirically: BiGG REST returns 404 for "R_TREH"
                    # but resolves "TREH" to the Trehalase reaction.
                    if bigg.startswith("R_"):
                        bigg = bigg[2:]
                    mnxr_to_bigg[mnxr_id] = bigg
    return seed_to_mnxr, mnxr_to_bigg, version


def build_mapping(
    seed_to_mnxr: dict[str, str],
    mnxr_to_bigg: dict[str, str],
    overrides: dict[str, str],
) -> dict[str, str]:
    """Chain SEED -> MNXR -> BiGG, applying overrides last (overrides win)."""
    result: dict[str, str] = {}
    for rxn, mnxr in seed_to_mnxr.items():
        bigg = mnxr_to_bigg.get(mnxr)
        if bigg:
            result[rxn] = bigg
    for rxn, bigg in overrides.items():
        result[rxn] = bigg
    return result


def write_tsv(mapping: dict[str, str], output: Path) -> None:
    """Write the two-column TSV in deterministic (sorted) order."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("seed_rxn_id\tbigg_rxn_id\n")
        for rxn in sorted(mapping):
            fh.write(f"{rxn}\t{mapping[rxn]}\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--reac-xref",
        type=Path,
        default=DEFAULT_REAC_XREF,
        help=f"Path to MetaNetX reac_xref.tsv (default: {DEFAULT_REAC_XREF})",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output TSV path (default: {DEFAULT_OUTPUT})",
    )
    args = p.parse_args(argv)

    if not args.reac_xref.exists():
        print(f"ERROR: {args.reac_xref} not found.", file=sys.stderr)
        print(
            "Download reac_xref.tsv from "
            "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",
            file=sys.stderr,
        )
        return 2

    print(
        f"Reading {args.reac_xref} ({args.reac_xref.stat().st_size // (1024**2)} MB)..."
    )
    seed_to_mnxr, mnxr_to_bigg, version = parse_reac_xref(args.reac_xref)
    print(f"  MNXref version     : {version}")
    print(f"  SEED -> MNXR rows  : {len(seed_to_mnxr):,}")
    print(f"  MNXR -> BiGG rows  : {len(mnxr_to_bigg):,}")

    mapping = build_mapping(seed_to_mnxr, mnxr_to_bigg, BIGG_RXN_OVERRIDES)
    chain_only = sum(1 for r in mapping if r not in BIGG_RXN_OVERRIDES)
    print(f"  SEED -> BiGG chain hits   : {chain_only:,}")
    print(f"  Manual overrides applied  : {len(BIGG_RXN_OVERRIDES)}")
    print(f"  TOTAL output rows         : {len(mapping):,}")

    write_tsv(mapping, args.output)
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
