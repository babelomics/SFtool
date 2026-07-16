# context.py
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import uuid

from sftool.core.config import Config
from sftool.variant_confirmation.models import (
    VariantConfirmationRequest
)


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
        self.sex: str = sample_data["sex"]
        self.vcf: Path = Path(sample_data["vcf_path"]).resolve()
        pgx_vcf_path = sample_data.get("pgx_vcf_path")
        self.pgx_vcf = Path(pgx_vcf_path).resolve() if pgx_vcf_path else None

        self.role: str = sample_data.get("relation", "proband")
        self.categories: List[str] = sample_data.get("categories", [])

        self.stripy_path = sample_data.get("stripy_path")
        self.smaca_path = sample_data.get("smaca_path")
        self.hpo_terms = sample_data.get("hpo_terms", [])

        variant_confirmation_data = sample_data.get(
            "variant_confirmation"
        )

        self.variant_confirmation_request: Optional[
            VariantConfirmationRequest
        ] = (
            VariantConfirmationRequest.from_dict(
                variant_confirmation_data
            )
            if variant_confirmation_data
            else None
        )
        # ------------------------------------------------------------------
        # Outputs populated during execution
        # ------------------------------------------------------------------
        self.vcf_outputs: dict = {
            "normalized": None,
            "intersected": {},
            "genebe_annotated": {},
            "PGx_preprocessed": None,
            "variant_confirmation": {
                "raw": None,
                "normalized": None,
                "genebe_annotated": None,
                "matches": None,
            },
        }
        self.reports: dict = {}
        self.variant_collections: dict = {
            "PR": {
                "snv_indels_genebe_clinvar": {},
            },
            "RR": {
                "snv_indels_genebe_clinvar": {},
                "STRs": {},
                "SMN1_copy": {},
            },
            "PGx": {
                "pharmCAT_variants": {}
            }
        }
        self.variant_selection: dict = {
            "PR": {
                "snv_indels_genebe_clinvar": {},
            },
            "RR": {
                "snv_indels_genebe_clinvar": {},
                "STRs": {},
                "SMN1_copy": {},
            },
            "PGx": {
                "pharmCAT_variants": {}
            }
        }
        self.results: Dict[str, object] = {
            "variant_confirmation": None
        }

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
        self.modes: List[str] = execution_meta.get("modes",[])
        self.clinvar_evidence: Optional[int] = execution_meta.get("clinvar_evidence")
        self.variant_classification_sources: list[str] = execution_meta.get("variant_classification_sources")
        self.RR_mode: Optional[str] = execution_meta.get("RR_mode")



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
        self.outputs: Dict[str, Dict[str, Dict[str, str | Path]]] = {
            "catalogs": {
                "bed_files": {},
                "json_files": {},
            },
            "clinvar": {
                "clinvar_db": "",
                "clinvar_summary_db": "",
                "clinvar_db_version": "",
                "clinvar_db_assembly": "",
                "PR_json": "",
                "RR_json": ""
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
            f"modes={self.modes}, "
            f"samples={len(self.samples)})"
        )
