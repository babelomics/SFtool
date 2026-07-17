from __future__ import annotations

import json
import os
import tempfile

from pathlib import Path
from typing import Optional, Sequence

from sftool.core.context import ExecutionContext, SampleContext
from sftool.utils.geneBe_utils import run_genebe
from sftool.variant_confirmation.matcher import CandidateMatchingOutput
from sftool.variant_confirmation.models import (
    VariantCandidate,
    VariantConfirmationRequest,
    VariantConfirmationResult,
)


VARIANT_CONFIRMATION_DIRECTORY = "variant_confirmation"
CONVERSION_FILENAME = "conversion.json"
GENEBE_FILENAME = "diagnostic_candidates.genebe.vcf.gz"


def get_variant_confirmation_output_dir(
        ctx: ExecutionContext,
        sample: SampleContext,
) -> Path:
    """
    Return and create the Variant Confirmation output directory for a sample.
    """
    output_dir = (
            Path(ctx.run_dir)
            / sample.sample_id
            / VARIANT_CONFIRMATION_DIRECTORY
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    return output_dir


def require_normalized_patient_vcf(
        sample: SampleContext,
) -> Path:
    """
    Return the normalized patient VCF required for exact candidate matching.

    Raises
    ------
    RuntimeError
        If sample preprocessing has not generated the normalized VCF.
    """
    normalized_vcf = sample.vcf_outputs.get(
        "normalized"
    )

    if normalized_vcf is None:
        raise RuntimeError(
            "Normalized patient VCF is missing for sample "
            f"{sample.sample_id!r}. Sample preprocessing must run "
            "before Variant Confirmation."
        )

    normalized_vcf = Path(
        normalized_vcf
    )

    if not normalized_vcf.is_file():
        raise RuntimeError(
            "Normalized patient VCF not found for sample "
            f"{sample.sample_id!r}: {normalized_vcf}"
        )

    return normalized_vcf


def write_conversion_json(
        request: VariantConfirmationRequest,
        candidates: Sequence[VariantCandidate],
        output_dir: str | Path,
) -> Path:
    """
    Serialize the parsed request and converted genomic candidates.

    The file is written atomically as ``conversion.json``.
    """
    output_dir = Path(
        output_dir
    )
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path = (
            output_dir
            / CONVERSION_FILENAME
    )

    payload = {
        "request": request.to_dict(),
        "candidates": [
            candidate.to_dict()
            for candidate in candidates
        ],
    }

    temporary_path = None

    try:
        with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=output_dir,
                prefix=f".{CONVERSION_FILENAME}.",
                suffix=".tmp",
                delete=False,
        ) as temporary_file:
            json.dump(
                payload,
                temporary_file,
                indent=2,
                ensure_ascii=False,
            )
            temporary_file.write("\n")

            temporary_path = Path(
                temporary_file.name
            )

        os.replace(
            temporary_path,
            output_path,
        )

    except (OSError, TypeError, ValueError) as exc:
        if (
                temporary_path is not None
                and temporary_path.exists()
        ):
            try:
                temporary_path.unlink()
            except OSError:
                pass

        raise RuntimeError(
            "Could not write Variant Confirmation conversion JSON "
            f"{output_path}: {exc}"
        ) from exc

    return output_path


def has_detected_candidates(
        matching_output: CandidateMatchingOutput,
) -> bool:
    """
    Return True when at least one genomic candidate was found.
    """
    return any(
        variant_match.found
        for variant_match in matching_output.matches
    )


def run_variant_confirmation_genebe(
        ctx: ExecutionContext,
        input_vcf: str | Path,
        output_dir: str | Path,
) -> Path:
    """
    Annotate detected diagnostic candidates using GeneBe.

    ``run_genebe`` is called with an explicit output path so the standard
    Variant Confirmation filename does not depend on PR/RR category naming.
    """
    output_path = (
            Path(output_dir)
            / GENEBE_FILENAME
    )

    return run_genebe(
        norm_vcf=input_vcf,
        category=None,
        assembly=ctx.assembly,
        genebe_path=ctx.config.paths.genebe,
        java_path=ctx.config.paths.java,
        api_key=ctx.config.genebe_credentials.api_key,
        username=ctx.config.genebe_credentials.username,
        tmp_dir=ctx.tmp_dir,
        output_file=output_path,
    )


def store_variant_confirmation_outputs(
        sample: SampleContext,
        conversion_json: Path,
        raw_candidate_vcf: Path,
        normalized_candidate_vcf: Path,
        matching_vcf: Path,
        genebe_annotated_vcf: Optional[Path],
        result_json: Path,
        result: VariantConfirmationResult,
) -> None:
    """
    Store generated paths and the structured result in SampleContext.
    """
    outputs = sample.vcf_outputs.setdefault(
        VARIANT_CONFIRMATION_DIRECTORY,
        {},
    )

    outputs.update(
        {
            "conversion": conversion_json,
            "raw": raw_candidate_vcf,
            "normalized": normalized_candidate_vcf,
            "matches": matching_vcf,
            "genebe_annotated": genebe_annotated_vcf,
            "result_json": result_json,
        }
    )

    sample.results[
        VARIANT_CONFIRMATION_DIRECTORY
    ] = result
