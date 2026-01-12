# context.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from datetime import datetime
import uuid


class SampleContext:
    """
    Represents a single biological sample in an SFtool execution.

    SampleContext contains:
        • sample-specific biological inputs
        • sample-specific tool / reference paths
        • a reference to ExecutionContext for run-level configuration

    IMPORTANT:
        Some fields (e.g. stripy_path, smaca_path, hpo_path) may also appear
        at run-level as defaults. Sample-level values always take precedence.
    """

    def __init__(self, sample_data: dict, exec_ctx: "ExecutionContext"):
        # ------------------------------------------------------------------
        # Required sample-level fields
        # ------------------------------------------------------------------
        self.sample_id: str = sample_data["sample_id"]
        self.vcf: Path = Path(sample_data["vcf_path"]).resolve()

        self.role: str = sample_data.get("relation", "proband")
        self.categories: List[str] = sample_data.get("categories", [])

        # Reference to run-level context
        self.exec: ExecutionContext = exec_ctx

        # ------------------------------------------------------------------
        # Sample-specific tool / reference paths
        #
        # Resolution rule:
        #   sample-level value > execution-level default
        # ------------------------------------------------------------------
        self.stripy_path: Optional[Path] = self._resolve_path(
            sample_data.get("stripy_path"),
            exec_ctx.stripy_path,
        )

        self.smaca_path: Optional[Path] = self._resolve_path(
            sample_data.get("smaca_path"),
            exec_ctx.smaca_path,
        )

        self.hpo_path: Optional[Path] = self._resolve_path(
            sample_data.get("hpo_path"),
            exec_ctx.hpo_path,
        )

        # ------------------------------------------------------------------
        # Outputs populated during execution
        # ------------------------------------------------------------------
        self.vcf_outputs: Dict[str, Path] = {}
        self.results: Dict[str, dict] = {}

        # ------------------------------------------------------------------
        # TODO (refactor):
        #   • Ensure validate_all() places sample-specific paths ONLY in samples
        #   • Clearly document which tools allow per-sample overrides
        #   • Add validation to ensure required paths are resolved
        # ------------------------------------------------------------------

    @staticmethod
    def _resolve_path(
            sample_value: Optional[str],
            execution_default: Optional[Path],
    ) -> Optional[Path]:
        """
        Resolve a path using sample-level override with execution-level fallback.
        """
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
    Run-level context built from samples_data returned by validate_all().

    ExecutionContext is the SINGLE SOURCE OF TRUTH for run-level invariants.

    It defines a run-specific directory layout under the user-provided
    output directory, using the generated run_id as directory name.
    """

    def __init__(
            self,
            samples_data: dict,
            output_dir: str,
    ):
        # ----------------------------------------------------------
        # Execution identifier (generated per run)
        # ----------------------------------------------------------
        self.run_id: str = (
                datetime.utcnow().strftime("%Y%m%d_%H%M%S")
                + "_"
                + uuid.uuid4().hex[:8]
        )

        # ------------------------------------------------------------------
        # Raw validated input (temporary, for backward compatibility)
        # ------------------------------------------------------------------
        self._samples_data = samples_data
        execution = samples_data.get("execution", {})

        # ------------------------------------------------------------------
        # Run-level identifiers and parameters
        # ------------------------------------------------------------------
        self.assembly: Optional[str] = execution.get("reference_genome")
        self.mode: Optional[str] = execution.get("mode")
        self.clinvar_evidence: Optional[int] = execution.get("clinvar_evidence")
        self.profile: Optional[str] = execution.get("profile")

        # ------------------------------------------------------------------
        # Output directory layout
        #
        # base_output_dir : user-provided (--outdir)
        # run_dir         : base_output_dir / <run_id>
        # tmp_dir         : run_dir / tmp
        # ------------------------------------------------------------------
        self.base_output_dir: Path = Path(output_dir).resolve()
        self.base_output_dir.mkdir(parents=True, exist_ok=True)

        self.run_dir: Path = self.base_output_dir / self.run_id
        self.run_dir.mkdir(parents=True, exist_ok=True)

        self.tmp_dir: Path = self.run_dir / "tmp"
        self.tmp_dir.mkdir(parents=True, exist_ok=True)

        # ------------------------------------------------------------------
        # Run-level resources
        # ------------------------------------------------------------------
        self.bed_files: Dict[str, Path] = {
            k: Path(v).resolve()
            for k, v in execution.get("bed_files", {}).items()
        }

        self.json_files: Dict[str, Path] = {
            k: Path(v).resolve()
            for k, v in execution.get("json_files", {}).items()
        }

        # ------------------------------------------------------------------
        # DEFAULT tool / reference paths (optional)
        #
        # These are defaults and may be overridden per sample.
        # ------------------------------------------------------------------
        self.stripy_path: Optional[Path] = (
            Path(execution.get("stripy_path")).resolve()
            if execution.get("stripy_path") else None
        )

        self.smaca_path: Optional[Path] = (
            Path(execution.get("smaca_path")).resolve()
            if execution.get("smaca_path") else None
        )

        self.hpo_path: Optional[Path] = (
            Path(execution.get("hpo_path")).resolve()
            if execution.get("hpo_path") else None
        )

        # ------------------------------------------------------------------
        # Sample contexts
        # ------------------------------------------------------------------
        self.samples: List[SampleContext] = [
            SampleContext(sample, self)
            for sample in samples_data.get("samples", [])
        ]

        # ------------------------------------------------------------------
        # TODO (refactor):
        #   • Remove self._samples_data once all consumers use contexts
        #   • Add explicit validation that required paths are resolved
        #   • Consider immutability for ExecutionContext
        # ------------------------------------------------------------------

    def get_sample(self, sample_id: str) -> Optional[SampleContext]:
        for sample in self.samples:
            if sample.sample_id == sample_id:
                return sample
        return None

    def __repr__(self) -> str:
        return (
            f"ExecutionContext(run_id={self.run_id}, "
            f"reference_genome={self.reference_genome}, "
            f"mode={self.mode}, "
            f"samples={len(self.samples)})"
        )
