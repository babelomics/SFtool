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

    # RR_mode rules
    if num_samples == 1:
        exec_data.setdefault("RR_mode", "auto")
        if exec_data["RR_mode"] != "auto":
            raise ValidationError("RR_mode must be 'auto' when only one sample is provided.")
    else:
        # 2 samples
        if exec_data.get("RR_mode") not in ["screening", "advanced"]:
            raise ValidationError("RR_mode must be 'screening' or 'advanced' for two samples.")


# =====================================================
# VALIDATE SAMPLES
# =====================================================
def validate_sample_block(samples, mode):
    if not isinstance(samples, list):
        raise ValidationError("samples must be a list")

    if len(samples) == 0 or len(samples) > 2:
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

    return data


def validate_all(samples_path, config_path):
    samples_data = validate_samples_info(samples_path)
    config_data = validate_config(config_path)
    return samples_data, config_data

