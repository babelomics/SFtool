"""

Execution bootstrap:
- validate samples_info and config files
- build Config and ExecutionContext object
- check runtime dependencies

"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, Any, List
from modules.context import ExecutionContext, SampleContext
from modules.config import Config
from modules.misc.runtime import check_runtime_dependencies
from modules.misc.errors import ValidationError
from modules.misc.vcf_utils import (
    validate_chr_prefix,
    check_vcf_positions_present
)

# =====================================================================
# Public API
# =====================================================================

def bootstrap_execution(
        samples_json: str | Path,
        config_json: str | Path,
        output_dir: str | Path,
        tmp_dir: str | Path | None = None,
) -> ExecutionContext:
    """
    Load, validate and normalize all inputs, returning an ExecutionContext.
    """

    samples_json = Path(samples_json)
    config_json = Path(config_json)
    output_dir = Path(output_dir)
    tmp_dir = Path(tmp_dir) if tmp_dir else None

    _validate_file_exists(samples_json, "samples_info JSON")
    _validate_file_exists(config_json, "config JSON")

    samples_info = _load_json(samples_json)
    config_data = _load_json(config_json)

    # -----------------------------------------------------------------
    # Legacy-equivalent validations (dict-level)
    # -----------------------------------------------------------------
    validate_samples_info(samples_info)
    validate_config(config_data, samples_info)

    # -----------------------------------------------------------------
    # Instantiate Config (run-level)
    # -----------------------------------------------------------------
    config = Config(config_data)
    check_runtime_dependencies(config.paths)

    # -----------------------------------------------------------------
    # Build ExecutionContext
    # -----------------------------------------------------------------
    ctx = ExecutionContext(
        execution_meta=samples_info["execution"],
        config=config,
        output_dir=output_dir,
        tmp_dir=tmp_dir,
    )

    # -----------------------------------------------------------------
    # Build SampleContexts
    # -----------------------------------------------------------------
    for sample_data in samples_info["samples"]:
        sample_ctx = SampleContext(sample_data=sample_data, exec_ctx=ctx)
        ctx.add_sample(sample_ctx)

    # -----------------------------------------------------------------
    # Cross-object semantic validation
    # -----------------------------------------------------------------
    validate_execution_context(ctx)

    return ctx


# =====================================================================
# JSON utilities
# =====================================================================

def _load_json(path: Path) -> Dict[str, Any]:
    try:
        with path.open() as fh:
            return json.load(fh)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in {path}: {e}") from e


def _validate_file_exists(path: Path, label: str):
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if not path.is_file():
        raise ValueError(f"{label} is not a file: {path}")


# =====================================================================
# Legacy-equivalent validation functions
# =====================================================================

def validate_samples_info(samples_info: Dict[str, Any]):
    """
    Equivalent to legacy validate_samples_info()
    """
    if "execution" not in samples_info:
        raise ValueError("samples_info must contain an 'execution' block")

    if "samples" not in samples_info:
        raise ValueError("samples_info must contain a 'samples' block")

    validate_execution_block(samples_info["execution"], num_samples=len(samples_info["samples"]))

    if not isinstance(samples_info["samples"], list) or not samples_info["samples"]:
        raise ValueError("'samples' must be a non-empty list")


    validate_sample_block(samples_info["samples"], samples_info["execution"]["mode"])


def validate_execution_block(exec_data: Dict[str, Any], num_samples):
    # -------- Required fields --------
    if "mode" not in exec_data:
        raise ValidationError("Missing required field: execution.mode")

    mode = exec_data["mode"]

    if mode not in ["secondary_findings_discovery", "variant_confirmation"]:
        raise ValidationError(f"Invalid execution.mode: {mode}")

    # -------- Set defaults --------
    exec_data.setdefault("reference_genome", "GRCh37")
    exec_data.setdefault("clinvar_evidence", 1)
    exec_data.setdefault("RR_mode", "screening")

    # -------- Validate reference genome --------
    if exec_data["reference_genome"] not in ["GRCh37", "GRCh38"]:
        raise ValidationError("execution.reference_genome must be GRCh37 or GRCh38")

    # -------- Validate ClinVar evidence --------
    ce = exec_data["clinvar_evidence"]
    if not isinstance(ce, int) or not (1 <= ce <= 5):
        raise ValidationError("execution.clinvar_evidence must be an integer between 1 and 5")

    # -------- MODE: VARIANT CONFIRMATION --------
    if mode == "variant_confirmation":
        if num_samples != 1:
            raise ValidationError("variant_confirmation mode requires exactly ONE sample.")

        vc = exec_data.get("variant_confirmation", {})
        if not vc.get("enabled", False):
            raise ValidationError("variant_confirmation.enabled must be true in variant_confirmation mode.")

        for field in ["chrom", "pos", "ref", "alt"]:
            if field not in vc or vc[field] in ["", None]:
                raise ValidationError(f"Missing field: variant_confirmation.{field}")

        return  # Nothing else applies to this mode

    # =====================================================
    # MODE: SECONDARY FINDINGS DISCOVERY
    # =====================================================
    # Set default profile
    exec_data.setdefault("profile", "advanced")

    if exec_data["profile"] not in ["basic", "advanced"]:
        raise ValidationError("execution.profile must be 'basic' or 'advanced'")


    # -------- RR_mode rules --------
    rr_mode = exec_data["RR_mode"]

    if num_samples == 1:
        # One-sample → only screening allowed
        if rr_mode != "screening":
            raise ValidationError(
                "RR_mode must be 'screening' when only one sample is provided."
            )
    elif num_samples == 2:
        # Two samples → screening OR advanced
        if rr_mode not in ["screening", "advanced"]:
            raise ValidationError(
                "RR_mode for two samples must be 'screening' or 'advanced'."
            )

def validate_sample_block(samples: Dict[str, Any], mode):
    if not isinstance(samples, list):
        raise ValidationError("samples must be a list")

    if not (1 <= len(samples) <= 2):
        raise ValidationError("samples must contain 1 or 2 entries")

    relations = [s["relation"] for s in samples]

    # =====================================================
    # MODE: VARIANT CONFIRMATION
    # =====================================================
    if mode == "variant_confirmation":
        if len(samples) != 1:
            raise ValidationError("variant_confirmation mode requires exactly one sample.")

        if relations[0] != "proband":
            raise ValidationError("In variant_confirmation mode, the only allowed relation is 'proband'.")

        # Validate required fields and optional paths
        s = samples[0]

        if "vcf_path" not in s:
            raise ValidationError("Sample missing field: vcf_path")

        if not os.path.exists(s["vcf_path"]):
            raise ValidationError(f"VCF file not found for sample {s['sample_id']}: {s['vcf_path']}")


        s.setdefault("hpo_terms", [])
        s.setdefault("stripy_path", "")
        s.setdefault("smaca_path", "")

        return  # Done after validating the single sample

    # =====================================================
    # MODE: SECONDARY FINDINGS
    # =====================================================
    if len(samples) == 1:
        # Only allowed relation = proband
        if relations[0] != "proband":
            raise ValidationError("For one-sample secondary findings mode, relation must be 'proband'.")

    elif len(samples) == 2:
        # Must be exactly parent1 + parent2
        if sorted(relations) != ["parent1", "parent2"]:
            raise ValidationError(
                "For two-sample secondary findings mode, relations must be exactly ['parent1', 'parent2']."
            )

    # Validate VCF paths and optional fields for all samples
    for s in samples:
        if "sample_id" not in s:
            raise ValidationError("Each sample must contain sample_id")

        if "relation" not in s:
            raise ValidationError("Each sample must contain relation")

        if "vcf_path" not in s:
            raise ValidationError(f"Sample {s.get('sample_id', '?')} missing vcf_path")

        if not os.path.exists(s["vcf_path"]):
            raise ValidationError(
                f"VCF file not found for sample {s['sample_id']}: {s['vcf_path']}"
            )
        if "sex" not in s:
            raise ValidationError("Missing required field 'sex' for sample " + s['sample_id'])

        sex = s["sex"]

        if not isinstance(sex, str):
            raise ValidationError(
                f"Invalid type for 'sex' in sample {s['sample_id']}: "
                f"expected string, got {type(sex).__name__}"
            )

        if sex not in {"male", "female", "unknown"}:
            raise ValueError(
                f"Invalid value for 'sex' in sample {s['sample_id']}: '{sex}'. "
                "Allowed values are: male, female, unknown"
            )

        # Optional fields
        s.setdefault("hpo_terms", [])
        s.setdefault("stripy_path", "")
        s.setdefault("smaca_path", "")



def validate_config(config: Dict[str, Any], samples_info: dict):
    required = [
        "paths", "references", "catalogs",
        "clinvar", "genebe_credentials", "smaca_thresholds"
    ]

    for key in required:
        if key not in config:
            raise ValidationError(f"Missing '{key}' block in config.json")

    # Shorthand variables
    references = config["references"]
    catalogs = config["catalogs"]
    clinvar = config["clinvar"]

    # ----------------------------------------------------
    # Validate reference genomes exist
    # ----------------------------------------------------
    if "genomes" not in references:
        raise ValidationError("Missing 'references.genomes' block in config.json")

    for name, path in references["genomes"].items():
        if not isinstance(path, str) or path == "":
            raise ValidationError(f"Invalid reference genome path for '{name}'")

        if not os.path.exists(path):
            raise ValidationError(
                f"Reference genome '{name}' does not exist at: {path}"
            )

    # ----------------------------------------------------
    # Validate gene_to_phenotype file existence
    # ----------------------------------------------------
    if "gene_to_phenotype_file" not in references:
        raise ValidationError("Missing 'references.gene_to_phenotype_file' in config.json")

    g2p = references["gene_to_phenotype_file"]

    if not os.path.exists(g2p):
        raise ValidationError(
            f"gene_to_phenotype_file does not exist: {g2p}"
        )

    # ----------------------------------------------------
    # Validate pharmCAT_positions_vcf file existence when PGx category exists
    # ----------------------------------------------------
    requested_categories = {
        category
        for sample in samples_info["samples"]
        for category in sample.get("categories", [])
    }

    if "PGx" in requested_categories:
        if "pharmCAT_positions_vcf" not in references:
            raise ValidationError("Missing 'references.pharmCAT_positions_vcf' in config.json")

        pharmCAT_positions = references["pharmCAT_positions_vcf"]

        if not os.path.exists(pharmCAT_positions):
            raise ValidationError(
                f"pharmCAT_positions_vcf does not exist: {pharmCAT_positions}"
            )




    # ----------------------------------------------------
    # Validate catalogs (PR, RR, STR, PGx) file existence
    # ----------------------------------------------------
    for label, path in catalogs.items():
        if not os.path.exists(path):
            raise ValidationError(
                f"Catalog file for '{label}' does not exist: {path}"
            )

    # ----------------------------------------------------
    # Validate clinvar.version follows YYYYMMDD
    # ----------------------------------------------------
    if "version" not in clinvar:
        raise ValidationError("Missing 'clinvar.version' field in config.json")

    version = clinvar["version"]

    # must be string of 8 digits
    if not (isinstance(version, str) and (version == "latest" or (len(version) == 8 and version.isdigit()))):
            raise ValidationError(
                f"Invalid clinvar.version '{version}'. Expected YYYYMMDD (8 digits) or 'latest' string for the downloading of the latest Clinvar version."
            )

    if version != "latest":
        # split into components
        year = int(version[0:4])
        month = int(version[4:6])
        day = int(version[6:8])

        # year range
        if not (2015 <= year <= 2030):
            raise ValidationError(
                f"Invalid clinvar.version year '{year}'. Must be 2000–2030."
            )

        # month range
        if not (1 <= month <= 12):
            raise ValidationError(
                f"Invalid clinvar.version month '{month:02d}'. Must be 01–12."
            )

        # days per month (default February 28, updated later if leap year)
        days_in_month = {
            1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30,
            7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31
        }

        # leap year adjustment
        if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
            days_in_month[2] = 29

        # day range
        if not (1 <= day <= days_in_month[month]):
            raise ValidationError(
                f"Invalid clinvar.version day '{day:02d}' for month {month:02d}."
            )

    # ----------------------------------------------------
    # Validate clinvar.db_path existence
    # ----------------------------------------------------
    if "db_path" not in clinvar:
        raise ValidationError("Missing 'clinvar.db_path' in config.json")

    if not os.path.exists(clinvar["db_path"]):
        raise ValidationError(
            f"ClinVar database path does not exist: {clinvar['db_path']}"
        )



# =====================================================================
# Context-level semantic validation
# =====================================================================

def validate_execution_context(ctx: ExecutionContext):
    """
    Validations that require fully instantiated context objects.
    """

    _validate_unique_sample_ids(ctx)
    _validate_vcf_paths(ctx)
    _validate_categories(ctx)
    _validate_pgx_inputs(ctx)



def _validate_unique_sample_ids(ctx: ExecutionContext):
    seen = set()
    for sample in ctx.samples:
        if sample.sample_id in seen:
            raise ValueError(f"Duplicated sample_id detected: {sample.sample_id}")
        seen.add(sample.sample_id)


def _validate_vcf_paths(ctx: ExecutionContext):
    for sample in ctx.samples:
        if not sample.vcf.exists():
            raise FileNotFoundError(
                f"VCF not found for sample '{sample.sample_id}': {sample.vcf}"
            )


def _validate_categories(ctx: ExecutionContext):
    valid_categories = {"PR", "RR", "PGx"}
    for sample in ctx.samples:
        for cat in sample.categories:
            if cat not in valid_categories:
                raise ValueError(
                    f"Invalid category '{cat}' for sample '{sample.sample_id}'. "
                    f"Valid categories: {sorted(valid_categories)}"
                )


def _validate_pgx_inputs(ctx: ExecutionContext):
    """
    Validate PGx-specific VCF requirements.

    Rules:
    - All VCFs must match execution.reference_genome.
    - If PGx is requested:
        - reference genome must be GRCh38.
        - the VCF used for PGx must have chr-prefixed chromosomes.
        - the VCF used for PGx must contain all PharmCAT positions.
    - If PGx is requested and two VCFs are provided:
        - both vcf_path and pgx_vcf_path must have chr-prefixed chromosomes.
    """
    for sample in ctx.samples:
        if "PGx" in sample.categories:

            if ctx.assembly != "GRCh38":
                raise ValidationError(
                    f"Sample {sample.sample_id}: PGx requires GRCh38. "
                    f"Current reference_genome is {ctx.assembly}."
                )

            pgx_vcf = sample.pgx_vcf or sample.vcf

            validate_chr_prefix(pgx_vcf)

            if sample.pgx_vcf:
                validate_chr_prefix(
                    vcf_path=sample.vcf
                )

            check_vcf_positions_present(
                input_vcf=pgx_vcf,
                required_vcf=ctx.config.references.pharmCAT_positions_vcf,
                output_file=ctx.tmp_dir / "PGx" / f"{sample.sample_id}.missing_pharmcat_positions.tsv"
            )