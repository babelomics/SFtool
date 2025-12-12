import json
import os
import sys

from modules.misc.errors import ValidationError

# =====================================================
# LOAD JSON
# =====================================================
def load_json(path):
    if not os.path.exists(path):
        raise ValidationError(f"JSON file not found: {path}")

    try:
        with open(path, "r") as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        raise ValidationError(f"Invalid JSON format in {path}: {e}")


# =====================================================
# VALIDATE EXECUTION BLOCK
# =====================================================
def validate_execution_block(exec_data, num_samples):
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


# =====================================================
# VALIDATE SAMPLES
# =====================================================
def validate_sample_block(samples, mode):
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

        s.setdefault("hpo_path", "")
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

        # Optional fields
        s.setdefault("hpo_path", "")
        s.setdefault("stripy_path", "")
        s.setdefault("smaca_path", "")


# =====================================================
# VALIDATE BOTH FILES
# =====================================================
def validate_samples_info(samples_json_path):
    data = load_json(samples_json_path)

    if "execution" not in data:
        raise ValidationError("Missing 'execution' block in samples_info.json")

    if "samples" not in data:
        raise ValidationError("Missing 'samples' block in samples_info.json")

    execution = data["execution"]
    samples = data["samples"]

    validate_execution_block(execution, num_samples=len(samples))
    validate_sample_block(samples, execution["mode"])

    return data


def validate_config(config_json_path):
    data = load_json(config_json_path)

    required = [
        "paths", "references", "catalogs",
        "clinvar", "genebe_credentials", "smaca_thresholds"
    ]

    for key in required:
        if key not in data:
            raise ValidationError(f"Missing '{key}' block in config.json")

    # Shorthand variables
    references = data["references"]
    catalogs = data["catalogs"]
    clinvar = data["clinvar"]

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
    if not isinstance(version, str) or len(version) != 8 or not version.isdigit():
        raise ValidationError(
            f"Invalid clinvar.version '{version}'. Expected YYYYMMDD (8 digits)."
        )

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

    return data


def validate_all(samples_path, config_path):
    samples_data = validate_samples_info(samples_path)
    config_data = validate_config(config_path)
    return samples_data, config_data

