"""Curation strategy skeletons."""

from typing import Any


class Curator:
    """Applies manual or automatic curation steps to models."""

    def __init__(self, mode: str = "auto"):
        self.mode = mode

    def curate(self, model: Any) -> Any:
        """Curate a model according to the chosen mode.

        TODO: implement automatic rules and hook in manual steps.
        """
        raise NotImplementedError
