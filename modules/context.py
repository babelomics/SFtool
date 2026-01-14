# context.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import uuid

from modules.config import Config


class SampleContext:
    """
    Represents a single biological sample in an SFtool execution.

    SampleContext contains:
        • sample-specific biological inputs
        • sample-specific tool / reference paths
        • a reference to ExecutionContext for run-level configuration

    Resolution rule:
        sample-level value > execution-level default
    """

    def __init__(self, sample_data: dict, exec_ctx: "ExecutionContext"):
        # ------------------------------------------------------------------
        # Required sample-level fields
        # ------------------------------------------------------------------
        self.sample_id: str = sample_data["sample_id"]
        self.vcf: Path = Path(sample_data["vcf_path"]).resolve()

        self.role: str = sample_data.get("relation", "proband")
        self.categories: List[str] = sample_data.get("categories", [])

        self.stripy_path = sample_data.get("stripy_path")
        self.smaca_path = sample_data.get("smaca_path")
        self.hpo_path = sample_data.get("hpo_path")

        # ------------------------------------------------------------------
        # Outputs populated during execution
        # ------------------------------------------------------------------
        self.vcf_outputs: Dict[str, Path] = {}
        self.results: Dict[str, dict] = {}

    @staticmethod
    def _resolve_path(
            sample_value: Optional[str],
            execution_default: Optional[Path],
    ) -> Optional[Path]:
        if sample_value:
            return Path(sample_value).resolve()
        return execution_default

    def __repr__(self) -> str:
        return (
            f"SampleContext(sample_id={self.sample_id}, "
            f"role={self.role}, categories={self.categories})"
        )


class ExecutionContext:
    """
    Run-level execution context.

    SINGLE SOURCE OF TRUTH for:
        • run-level parameters
        • directory layout
        • configuration
    """

    def __init__(
            self,
            execution_meta: dict,
            config: Config,
            output_dir: str | Path,
            tmp_dir: str | Path | None = None,
    ):
        # ------------------------------------------------------------------
        # Execution identifier
        # ------------------------------------------------------------------
        self.run_id: str = (
                datetime.utcnow().strftime("%Y%m%d_%H%M%S")
                + "_"
                + uuid.uuid4().hex[:8]
        )

        # ------------------------------------------------------------------
        # Configuration (run-level)
        # ------------------------------------------------------------------
        self.config: Config = config

        # ------------------------------------------------------------------
        # Run-level parameters
        # ------------------------------------------------------------------
        self.assembly: Optional[str] = execution_meta.get("reference_genome")
        self.mode: Optional[str] = execution_meta.get("mode")
        self.clinvar_evidence: Optional[int] = execution_meta.get("clinvar_evidence")
        self.profile: Optional[str] = execution_meta.get("profile")
        self.RR_mode: Optional[str] = execution_meta.get("RR_mode")
        self.variant_confirmation: Optional[str] = execution_meta.get("variant_confirmation")

        # ------------------------------------------------------------------
        # Output directory layout
        # ------------------------------------------------------------------
        self.base_output_dir: Path = Path(output_dir).resolve()
        self.base_output_dir.mkdir(parents=True, exist_ok=True)

        self.run_dir: Path = self.base_output_dir / self.run_id
        self.run_dir.mkdir(parents=True, exist_ok=True)

        self.tmp_dir: Path = (
            Path(tmp_dir).resolve()
            if tmp_dir
            else (self.run_dir / "tmp")
        )
        self.tmp_dir.mkdir(parents=True, exist_ok=True)

        # ------------------------------------------------------------------
        # Outputs populated during execution
        # ------------------------------------------------------------------
        self.outputs: Dict[str, Dict[str, Dict[str, Path]]] = {
            "catalogs": {
                "bed_files": {},
                "json_files": {},
            }
        }

        # ------------------------------------------------------------------
        # Sample contexts (populated by validation)
        # ------------------------------------------------------------------
        self.samples: List[SampleContext] = []

    def add_sample(self, sample: SampleContext):
        self.samples.append(sample)

    def get_sample(self, sample_id: str) -> Optional[SampleContext]:
        for sample in self.samples:
            if sample.sample_id == sample_id:
                return sample
        return None

    def __repr__(self) -> str:
        return (
            f"ExecutionContext(run_id={self.run_id}, "
            f"assembly={self.assembly}, "
            f"mode={self.mode}, "
            f"samples={len(self.samples)})"
        )
