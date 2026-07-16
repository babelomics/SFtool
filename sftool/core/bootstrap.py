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
from typing import Dict, Any
from sftool.core.context import ExecutionContext, SampleContext
from sftool.core.config import Config
from sftool.utils.runtime import check_runtime_dependencies
from sftool.utils.errors import ValidationError
from sftool.utils.vcf_utils import (
    validate_chr_prefix,
    check_vcf_positions_present
)


SUPPORTED_EXECUTION_MODES = {
    "secondary_findings_discovery",
    "variant_confirmation",
}

SUPPORTED_SAMPLE_GENDER = {
    "male",
    "female",
    "unknown"
}

SUPPORTED_SAMPLE_ROLES = {
    "proband",
    "parent1",
    "parent2"
}

SUPPORTED_VARIANT_CLASSIFICATION_SOURCES = [
    ["genebe"],
    ["genebe", "clinvar"]

]

SUPPORTED_REFERENCE_GENOMES = {
    "GRCh37",
    "GRCh38"
}

SUPPORTED_RR_MODES = {
    "screening",
    "advanced"
}
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

    validate_file_exists(samples_json, "samples_info JSON")
    validate_file_exists(config_json, "config JSON")

    samples_info = load_json(samples_json)
    config_data = load_json(config_json)

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

def load_json(path: Path) -> Dict[str, Any]:
    try:
        with path.open() as fh:
            return json.load(fh)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in {path}: {e}") from e


def validate_file_exists(path: Path, label: str):
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
    if not isinstance(samples_info, dict):
        raise ValidationError(
            "samples_info must be a JSON object"
        )

    if "execution" not in samples_info:
        raise ValueError("samples_info must contain an 'execution' block")

    if "samples" not in samples_info:
        raise ValueError("samples_info must contain a 'samples' block")

    execution = samples_info["execution"]
    samples = samples_info["samples"]

    if not isinstance(execution, dict):
        raise ValidationError(
            "'execution' must be an object"
        )

    if not isinstance(samples, list):
        raise ValidationError(
            "'samples' must be a list"
        )

    if not 1 <= len(samples) <= 2:
        raise ValidationError(
            "'samples' must contain 1 or 2 entries"
        )

    validate_execution_block(execution, num_samples=len(samples))

    validate_sample_block(samples, modes=execution["modes"])


def validate_execution_block(exec_data: Dict[str, Any], num_samples):
    # -------- Required fields --------
    if "modes" not in exec_data:
        raise ValidationError("Missing required field: execution.modes")

    modes = exec_data["modes"]

    if not isinstance(modes, list):
        raise ValidationError(
            "execution.modes must be a list"
        )

    if not modes:
        raise ValidationError(
            "execution.modes must contain at least one workflow"
        )

    if any(not isinstance(mode, str) for mode in modes):
        raise ValidationError(
            "Every execution mode must be a string"
        )

    if len(modes) != len(set(modes)):
        raise ValidationError(
            "execution.modes must not contain duplicated workflows"
        )

    unsupported_modes = set(modes) - SUPPORTED_EXECUTION_MODES

    if unsupported_modes:
        raise ValidationError(
            "Unsupported execution mode(s): "
            + ", ".join(sorted(unsupported_modes))
        )

    # With two samples, secondary_findings_discovery is mandatory
    if num_samples == 2 and "secondary_findings_discovery" not in modes:
        raise ValidationError(
            "Two-sample executions require 'secondary_findings_discovery'."
        )

    # -------- Set defaults --------
    exec_data.setdefault("reference_genome", "GRCh37")
    exec_data.setdefault("clinvar_evidence", 1)

    # -------- Validate reference genome --------
    if exec_data["reference_genome"] not in SUPPORTED_REFERENCE_GENOMES:
        raise ValidationError("execution.reference_genome must be GRCh37 or GRCh38")

    # -------- Validate ClinVar evidence --------
    ce = exec_data["clinvar_evidence"]
    if not isinstance(ce, int) or not (1 <= ce <= 5) or isinstance(ce, bool):
        raise ValidationError("execution.clinvar_evidence must be an integer between 1 and 5")

    # =====================================================
    # MODE: SECONDARY FINDINGS DISCOVERY
    # =====================================================

    if "secondary_findings_discovery" in modes:
        # Set default variant_classification_sources
        exec_data.setdefault("variant_classification_sources", ["genebe", "clinvar"])
        exec_data.setdefault("RR_mode", "screening")

        if exec_data["variant_classification_sources"] not in SUPPORTED_VARIANT_CLASSIFICATION_SOURCES:
            raise ValidationError("execution.variant_classification_sources must be ['genebe'] or ['genebe', 'clinvar']")

        # -------- RR_mode rules --------
        rr_mode = exec_data["RR_mode"]

        # One-sample → only screening allowed
        if num_samples == 1 and rr_mode != "screening":
                raise ValidationError("RR_mode must be 'screening' when only one sample is provided.")
        # Two samples → screening OR advanced
        if num_samples == 2 and rr_mode not in SUPPORTED_RR_MODES:
                raise ValidationError("RR_mode for two samples must be 'screening' or 'advanced'.")

def validate_sample_block(samples: Dict[str, Any], modes: list[str]):
    if not isinstance(samples, list):
        raise ValidationError("samples must be a list")

    if not (1 <= len(samples) <= 2):
        raise ValidationError("samples must contain 1 or 2 entries")

    for sample in samples:
        validate_common_sample_fields(sample)

    if "secondary_findings_discovery" in modes:
        validate_secondary_findings_samples(samples)

    validate_variant_confirmation_requests(samples, modes)


def validate_common_sample_fields(sample: Dict[str, Any]):

    ### Sample validation: relation, sex, vcf_path

    if not isinstance(sample, dict):
        raise ValidationError(
            "Each sample must be an object"
        )

    sample_id = sample.get("sample_id")

    if not isinstance(sample_id, str) or not sample_id.strip():
        raise ValidationError(
            "Each sample must contain "
            "a non-empty sample_id"
        )

    relation = sample.get("relation")

    if relation not in SUPPORTED_SAMPLE_ROLES:
        raise ValidationError(
            f"Invalid relation for sample "
            f"{sample_id}: {relation}"
        )

    if "vcf_path" not in sample:
        raise ValidationError(
            f"Sample {sample_id} missing vcf_path"
        )

    vcf_path = sample["vcf_path"]

    if not isinstance(vcf_path, str) or not vcf_path.strip():
        raise ValidationError(
            f"Sample {sample_id}: "
            "vcf_path must be a non-empty string"
        )

    if not os.path.exists(vcf_path):
        raise ValidationError(
            f"VCF file not found for sample "
            f"{sample_id}: {vcf_path}"
        )

    if "sex" not in sample:
        raise ValidationError(
            f"Missing required field 'sex' "
            f"for sample {sample_id}"
        )

    sex = sample["sex"]

    if not isinstance(sex, str):
        raise ValidationError(
            f"Invalid type for 'sex' in sample "
            f"{sample_id}: expected string, "
            f"got {type(sex).__name__}"
        )

    if sex not in SUPPORTED_SAMPLE_GENDER:
        raise ValidationError(
            f"Invalid value for 'sex' in sample "
            f"{sample_id}: '{sex}'. "
            "Allowed values are: "
            "male, female, unknown"
        )

    sample.setdefault("hpo_terms", [])
    sample.setdefault("stripy_path", "")
    sample.setdefault("smaca_path", "")
    sample.setdefault("pgx_vcf_path", "")
    sample.setdefault("categories", [])


def validate_secondary_findings_samples(samples: list[Dict[str, Any]]):

    ### Validation relationship rules for secondary findings mode

    relations = [
        sample["relation"]
        for sample in samples
    ]

    if len(samples) == 1:
        if relations[0] != "proband":
            raise ValidationError(
                "For one-sample secondary findings "
                "execution, relation must be 'proband'."
            )

    elif len(samples) == 2:
        if sorted(relations) != [
            "parent1",
            "parent2",
        ]:
            raise ValidationError(
                "For two-sample secondary findings "
                "execution, relations must be exactly "
                "['parent1', 'parent2']."
            )


def validate_variant_confirmation_requests(samples: list[Dict[str, Any]], modes: list[str]):

    #### Validate variant_confirmation mode. The following rules are implemented:
    #       1. variant_confirmation is optional per sample
    #       2. if variant_confirmation mode is present, at least one sample must contain a variant
    #       3. variant must be a non-empty string (no spaces)
    #       4. if variant confirmation is not present, a variant cannot be included for any sample

    confirmation_enabled = ("variant_confirmation" in modes)

    samples_with_request = [
        sample
        for sample in samples
        if "variant_confirmation" in sample
    ]

    if confirmation_enabled and not samples_with_request:
        raise ValidationError(
            "At least one sample must define "
            "'variant_confirmation' when "
            "'variant_confirmation' is enabled "
            "in execution.modes."
        )

    if not confirmation_enabled and samples_with_request:
        sample_ids = [
            sample.get("sample_id", "?")
            for sample in samples_with_request
        ]

        raise ValidationError(
            "variant_confirmation was provided "
            "for sample(s) "
            f"{', '.join(sample_ids)}, "
            "but the workflow is not enabled "
            "in execution.modes."
        )

    for sample in samples_with_request:
        sample_id = sample["sample_id"]
        request = sample["variant_confirmation"]

        if not isinstance(request, dict):
            raise ValidationError(
                f"Sample {sample_id}: "
                "variant_confirmation must be "
                "an object."
            )

        if set(request.keys()) != {"variant"}:
            raise ValidationError(
                f"Sample {sample_id}: "
                "variant_confirmation must contain "
                "exactly one field: 'variant'."
            )

        variant = request["variant"]

        if not isinstance(variant, str) or not variant.strip():
            raise ValidationError(
                f"Sample {sample_id}: "
                "variant_confirmation.variant must "
                "be a non-empty string."
            )

        request["variant"] = variant.strip()

def validate_config(config: Dict[str, Any], samples_info: dict | None = None):
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

    requested_categories = set()

    if samples_info is not None:
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