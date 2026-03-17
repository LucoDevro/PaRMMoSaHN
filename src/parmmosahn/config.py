"""Configuration models for pangenome reconstruction workflows.

This module defines simple dataclasses to pass around filesystem paths and
parameter settings instead of relying on hardcoded values inside notebooks.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class PipelineConfig:
    """Holds high-level configuration for the full pipeline.

    TODO: refine fields with more specific types or enums where appropriate, and add validation logic if needed.
    """

    genomes_dir: Path
    pangenome_dir: Path
    models_dir: Path
    results_dir: Path
    media_dir: Optional[Path] = None
    random_seed: int = 42
    tmp_dir: Optional[Path] = None
    extra: dict = field(default_factory=dict)

    def ensure_dirs(self) -> None:
        """Create required directories if they do not exist."""
        for path in [self.pangenome_dir, self.models_dir, self.results_dir]:
            path.mkdir(parents=True, exist_ok=True)
        if self.tmp_dir:
            self.tmp_dir.mkdir(parents=True, exist_ok=True)
