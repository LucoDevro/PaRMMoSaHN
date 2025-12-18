"""Pipeline orchestrator skeleton."""

from pathlib import Path

from package.config import PipelineConfig
from package.curation import Curator
from package.io import FileRepository
from package.models import PanModelBuilder, StrainModelDeriver
from package.pangenome import PangenomeBuilder
from package.results import ResultsAggregator


class PipelineRunner:
    """Coordinates the end-to-end pangenome → models → curation → results flow."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.repo = FileRepository(config.genomes_dir)
        self.pangenome_builder = PangenomeBuilder(config.genomes_dir)
        self.pan_model_builder = PanModelBuilder()
        self.strain_deriver = StrainModelDeriver()
        self.curator = Curator(mode="auto")
        self.aggregator = ResultsAggregator()

    def run_all(self) -> Path:
        """Execute the full pipeline.

        TODO: implement orchestration steps and return path to results bundle.
        """
        raise NotImplementedError

    def build_pangenome(self):  # type: ignore[override]
        raise NotImplementedError

    def build_pan_model(self):  # type: ignore[override]
        raise NotImplementedError

    def derive_strain_models(self):  # type: ignore[override]
        raise NotImplementedError

    def curate_models(self):  # type: ignore[override]
        raise NotImplementedError

    def aggregate_results(self):  # type: ignore[override]
        raise NotImplementedError
