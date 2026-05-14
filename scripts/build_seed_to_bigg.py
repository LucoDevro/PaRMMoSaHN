#!/usr/bin/env python3
"""Build src/parmmosahn/data/seed_to_bigg.tsv from MetaNetX chem_xref.tsv.

Chains SEED compound ID -> MNXref MNXM ID -> BiGG metabolite ID using the
official MetaNetX cross-reference table. Applies a small set of manual
overrides for stereoisomer/anomer pairs where SEED and BiGG map to different
MNXref IDs for the same named compound (these would otherwise drop silently).

Usage
-----
    python scripts/build_seed_to_bigg.py \\
        --chem-xref ~/thesis/data/external/chem_xref.tsv \\
        --output    src/parmmosahn/data/seed_to_bigg.tsv

The chem_xref.tsv file is 648 MB and is NOT bundled in the wheel; download
it from https://www.metanetx.org/mnxdoc/mnxref.html (current is MNXref 4.5,
CC-BY 4.0, dated 2025-08-13). The thesis Methods section (3.3.5) should
cite this script as the source of the bundled TSV.

The output TSV is two-column (seed_cpd_id, bigg_met_id) with a header. The
script is idempotent: rerunning on the same chem_xref.tsv produces a
byte-identical output.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

DEFAULT_CHEM_XREF = Path.home() / "thesis/data/external/chem_xref.tsv"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parent.parent / "src/parmmosahn/data/seed_to_bigg.tsv"
)

# Stereoisomer / anomer / tautomer overrides. SEED and BiGG sometimes pin
# their canonical entry to different MNXM IDs even though the named compound
# is the same; the SEED->MNX->BiGG chain drops silently in those cases.
# Each entry was supervisor-verified or manually checked against the BiGG
# universal/metabolites REST endpoint. Sourced from
# parmmosahn-thesis archive build_substrate_mapping.py BIGG_OVERRIDES.
BIGG_OVERRIDES: dict[str, str] = {
    "cpd00082": "fru",  # D-Fructose: SEED MNXM1363862 vs BiGG MNXM1364098
    "cpd00108": "gal",  # D-Galactose: anomer split
    "cpd00224": "arab__L",  # L-Arabinose
    "cpd00208": "lcts",  # D-Lactose
    "cpd00396": "rmn",  # L-Rhamnose (supervisor-corrected SEED ID)
    "cpd01061": "raffin",  # D-Raffinose
    "cpd00181": "starch",  # Starch
    "cpd00794": "tre",  # D-Trehalose
    "cpd03198": "melib",  # D-Melibiose (alpha/beta anomer)
    "cpd00089": "g1p",  # D-Glucose-1-Phosphate
    "cpd00280": "galur",  # D-Galacturonate
    "cpd00453": "12ppd__R",  # 1,2-Propanediol (R-enantiomer in BiGG)
    "cpd00320": "asp__D",  # D-Aspartic Acid
    "cpd00611": "thr__D",  # D-Threonine
}


def parse_chem_xref(path: Path) -> tuple[dict[str, str], dict[str, str], str | None]:
    """Stream chem_xref.tsv once and build the SEED->MNX and MNX->BiGG dicts.

    For MNX->BiGG, the FIRST entry encountered for a given MNX is kept;
    MetaNetX lists canonical BiGG forms first in source order, so this
    matches the convention used by build_substrate_mapping.py.

    Returns
    -------
    seed_to_mnx : {cpd00027: MNXM41}
    mnx_to_bigg : {MNXM41: glc__D}
    version     : MNXref version string from the file header, or None
    """
    seed_to_mnx: dict[str, str] = {}
    mnx_to_bigg: dict[str, str] = {}
    version: str | None = None

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#VERSION:") and version is None:
                # First #VERSION: line is MetaNetX itself; later #VERSION: lines
                # are per-source-resource versions (KEGG, BiGG, SEED, etc.)
                version = line.split(":", 1)[1].strip()
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            source, mnx_id = parts[0], parts[1]
            if source.startswith("seed.compound:"):
                seed_to_mnx[source.split(":", 1)[1]] = mnx_id
            elif source.startswith("bigg.metabolite:"):
                if mnx_id not in mnx_to_bigg:
                    mnx_to_bigg[mnx_id] = source.split(":", 1)[1]
    return seed_to_mnx, mnx_to_bigg, version


def build_mapping(
    seed_to_mnx: dict[str, str],
    mnx_to_bigg: dict[str, str],
    overrides: dict[str, str],
) -> dict[str, str]:
    """Chain SEED -> MNX -> BiGG, applying overrides last (overrides win)."""
    result: dict[str, str] = {}
    for cpd, mnx in seed_to_mnx.items():
        bigg = mnx_to_bigg.get(mnx)
        if bigg:
            result[cpd] = bigg
    for cpd, bigg in overrides.items():
        result[cpd] = bigg
    return result


def write_tsv(mapping: dict[str, str], output: Path) -> None:
    """Write the two-column TSV in deterministic (sorted) order."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("seed_cpd_id\tbigg_met_id\n")
        for cpd in sorted(mapping):
            fh.write(f"{cpd}\t{mapping[cpd]}\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--chem-xref",
        type=Path,
        default=DEFAULT_CHEM_XREF,
        help=f"Path to MetaNetX chem_xref.tsv (default: {DEFAULT_CHEM_XREF})",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output TSV path (default: {DEFAULT_OUTPUT})",
    )
    args = p.parse_args(argv)

    if not args.chem_xref.exists():
        print(f"ERROR: {args.chem_xref} not found.", file=sys.stderr)
        print(
            "Download chem_xref.tsv from https://www.metanetx.org/mnxdoc/mnxref.html",
            file=sys.stderr,
        )
        return 2

    print(
        f"Reading {args.chem_xref} ({args.chem_xref.stat().st_size // (1024**2)} MB)..."
    )
    seed_to_mnx, mnx_to_bigg, version = parse_chem_xref(args.chem_xref)
    print(f"  MNXref version    : {version}")
    print(f"  SEED -> MNX rows  : {len(seed_to_mnx):,}")
    print(f"  MNX  -> BiGG rows : {len(mnx_to_bigg):,}")

    mapping = build_mapping(seed_to_mnx, mnx_to_bigg, BIGG_OVERRIDES)
    chain_only = sum(1 for c in mapping if c not in BIGG_OVERRIDES)
    print(f"  SEED -> BiGG chain hits   : {chain_only:,}")
    print(f"  Manual overrides applied  : {len(BIGG_OVERRIDES)}")
    print(f"  TOTAL output rows         : {len(mapping):,}")

    write_tsv(mapping, args.output)
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
