"""Results aggregation skeleton."""

from dataclasses import dataclass
from typing import Any, List


@dataclass
class ResultsBundle:
    strain_models: List[Any]
    summaries: Any  # TODO: type for summary tables/figures


class ResultsAggregator:
    """Collects and summarizes outputs across strains."""

    def aggregate(self, models: List[Any]) -> ResultsBundle:
        """Aggregate per-strain outputs into a bundle.

        TODO: implement collection of stats, tables, and plots.
        """
        raise NotImplementedError
