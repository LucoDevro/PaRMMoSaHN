"""Pangenome domain models and builder skeletons."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, List


@dataclass
class GenomeRecord:
    accession: str
    path: Path
    metadata: dict


@dataclass
class Pangenome:
    genomes: List[GenomeRecord]
    graph: Any  # TODO: replace Any with concrete graph structure
    annotations: Any = None  # TODO: detail annotation schema


class PangenomeBuilder:
    """Constructs a pangenome object from genome inputs."""

    def __init__(self, genomes_dir: Path):
        self.genomes_dir = genomes_dir

    def discover_genomes(self) -> List[GenomeRecord]:
        """Find genome files and return records.

        TODO: implement file discovery for FASTA/GBFF and capture metadata.
        """
        raise NotImplementedError

    def build(self) -> Pangenome:
        """Build the pangenome graph/structure from genome inputs.

        TODO: port clustering/graph construction logic from notebooks.
        """
        raise NotImplementedError
