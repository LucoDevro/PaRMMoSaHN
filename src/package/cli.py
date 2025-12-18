"""CLI entrypoints (skeleton).

The CLI is intentionally minimal until pipeline steps are implemented. Commands
currently emit TODO messages instead of executing logic.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from package.config import PipelineConfig
from package.pipeline import PipelineRunner


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="PaRMMoSaHN CLI (WIP)")
    sub = parser.add_subparsers(dest="command", required=True)

    run_all = sub.add_parser("run-all", help="Run full pipeline (WIP)")
    run_all.add_argument("genomes_dir", type=Path)
    run_all.add_argument("output_dir", type=Path)

    sub.add_parser("build-pangenome", help="Build pangenome only (WIP)")
    sub.add_parser("build-pan-model", help="Build pan-model only (WIP)")
    sub.add_parser("derive-strain", help="Derive strain models (WIP)")
    sub.add_parser("curate", help="Curate models (WIP)")
    sub.add_parser("aggregate", help="Aggregate results (WIP)")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "run-all":
        config = PipelineConfig(
            genomes_dir=args.genomes_dir,
            pangenome_dir=args.output_dir / "01_pangenome",
            models_dir=args.output_dir / "02_models",
            results_dir=args.output_dir / "03_results",
        )
        runner = PipelineRunner(config)
        # TODO: implement PipelineRunner.run_all when logic is migrated.
        print("TODO: implement run-all. Config:", config)
        return 0

    print(f"TODO: command '{args.command}' is not implemented yet.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
