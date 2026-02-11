"""Model-building skeletons for pan and strain-specific models."""

from dataclasses import dataclass
from typing import Any


@dataclass
class PanModel:
    graph: Any  # TODO: replace with actual model type (e.g., cobra.Model)
    metadata: dict


@dataclass
class StrainModel:
    id: str
    model: Any  # TODO: actual model type
    metadata: dict


class PanModelBuilder:
    """Builds a pan-model from a pangenome."""

    def __init__(self, media: Any = None):
        self.media = media

    def build(self, pangenome) -> PanModel:  # type: ignore[arg-type]
        """Construct the pan-model from the pangenome structure.

        TODO: port stoichiometric/constraint-building logic from notebooks.
        """
        raise NotImplementedError


class StrainModelDeriver:
    """Derives strain-specific models from a pan-model."""

    def __init__(self, strategy: str = "default"):
        self.strategy = strategy

    def derive(self, pan_model: PanModel, strain_id: str) -> StrainModel:
        """Slice or specialize the pan-model into a strain-specific model.

        TODO: implement derivation logic from notebooks.
        """
        raise NotImplementedError
