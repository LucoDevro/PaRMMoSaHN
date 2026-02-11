"""I/O helpers for loading and saving inputs/outputs.

These are intentionally minimal skeletons; wire them to the actual formats from
the notebooks (FASTA/GBFF/HDF5/CSV) when migrating logic.
"""

from pathlib import Path
from typing import Any


class FileRepository:
    """Handles reads/writes to the filesystem.

    TODO: Implement concrete loaders/savers based on notebook logic.
    """

    def __init__(self, base_dir: Path):
        self.base_dir = base_dir

    def load_fasta(self, path: Path) -> Any:
        raise NotImplementedError("TODO: implement FASTA parsing")

    def load_gbff(self, path: Path) -> Any:
        raise NotImplementedError("TODO: implement GenBank parsing")

    def load_h5(self, path: Path) -> Any:
        raise NotImplementedError("TODO: implement HDF5 reading")

    def save_h5(self, obj: Any, path: Path) -> None:
        raise NotImplementedError("TODO: implement HDF5 writing")

    def load_csv(self, path: Path) -> Any:
        raise NotImplementedError("TODO: implement CSV loading (e.g., media)")

    def save_csv(self, obj: Any, path: Path) -> None:
        raise NotImplementedError("TODO: implement CSV writing")
