from pathlib import Path
from typing import Dict


class SampleContext:
    """
    Represents a single biological sample in an SFtool execution.
    """

    def __init__(self, sample_data: dict):
        self.sample_id: str = sample_data["sample_id"]
        self.vcf: Path = Path(sample_data["vcf"]).resolve()
        self.role: str = sample_data.get("relation", "proband")

        # Populated during execution
        self.vcf_outputs: Dict[str, Path] = {}
        self.results: Dict[str, dict] = {}

    def __repr__(self) -> str:
        return f"SampleContext(sample_id={self.sample_id}, role={self.role})"


class ExecutionContext:
    """
    Run-level context built from samples_data returned by validate_all().
    """

    def __init__(
            self,
            samples_data: dict,
            output_dir: str,
            tmp_dir: str | None = None,
    ):
        self.execution = samples_data["execution"]
        self.samples_data = samples_data

        # Paths
        self.output_dir = Path(output_dir).resolve()
        self.tmp_dir = (
            Path(tmp_dir).resolve()
            if tmp_dir
            else self.output_dir / "tmp"
        )

        # Samples (1 or 2)
        self.samples: Dict[str, SampleContext] = {
            s["sample_id"]: SampleContext(s)
            for s in samples_data["samples"]
        }

        # Runtime artefacts
        self.bed_files: Dict[str, Path] = {}
        self.json_files: Dict[str, Path] = {}
        self.versions: Dict[str, str] = {}

        self._init_directories()

    # ---------- directories ----------

    def _init_directories(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.tmp_dir.mkdir(parents=True, exist_ok=True)

        for sample_id in self.samples:
            (self.output_dir / sample_id).mkdir(parents=True, exist_ok=True)

    # ---------- execution-level properties ----------

    @property
    def assembly(self) -> str:
        return self.execution["reference_genome"]

    def iter_samples(self):
        return self.samples.values()

    def __repr__(self) -> str:
        return f"ExecutionContext(samples={len(self.samples)}, assembly={self.assembly})"
