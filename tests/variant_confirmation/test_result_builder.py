import json

from pathlib import Path

import pytest

from sftool.variant_confirmation.matcher import CandidateMatchingOutput
from sftool.variant_confirmation.models import (
    VariantCandidate,
    VariantConfirmationRequest,
    VariantMatch,
)
from sftool.variant_confirmation.result_builder import (
    InvalidResultBuilderInputError,
    VariantConfirmationResultBuilder,
)


def build_request(
        variant: str = "NM_000251.3:c.2030C>A",
        representation_type: str = "hgvsc",
) -> VariantConfirmationRequest:
    """
    Build a parsed variant confirmation request for result-builder tests.
    """
    return VariantConfirmationRequest(
        variant=variant,
        representation_type=representation_type,
    )


def build_normalized_candidate(
        candidate_id: str = "candidate_1",
        chromosome: str = "2",
        position: int = 47476391,
        reference: str = "C",
        alternate: str = "A",
        assembly: str = "GRCh38",
) -> VariantCandidate:
    """
    Build a normalized diagnostic candidate for result-builder tests.
    """
    candidate = VariantCandidate(
        candidate_id=candidate_id,
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
        assembly=assembly,
    )

    candidate.set_normalized_coordinates(
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
    )

    return candidate


def build_variant_match(
        candidate: VariantCandidate,
        found: bool,
        genotype: str | None = "0/1",
        quality: float | None = 99.0,
        filters: list[str] | None = None,
        sample_format: dict | None = None,
) -> VariantMatch:
    """
    Build a match aligned with the normalized coordinates of a candidate.
    """
    variant_match = VariantMatch(
        candidate_id=candidate.candidate_id,
        chromosome=candidate.normalized_chromosome,
        position=candidate.normalized_position,
        reference=candidate.normalized_reference,
        alternate=candidate.normalized_alternate,
    )

    if found:
        variant_match.set_match_data(
            genotype=genotype,
            quality=quality,
            filters=filters if filters is not None else ["PASS"],
            sample_format=(
                sample_format
                if sample_format is not None
                else {
                    "GT": genotype,
                    "DP": 42,
                }
            ),
        )

    return variant_match


def build_matching_output(
        matches: list[VariantMatch],
        tmp_path: Path,
) -> CandidateMatchingOutput:
    """
    Build the matcher output wrapper without invoking bcftools.
    """
    return CandidateMatchingOutput(
        matches=matches,
        vcf_path=(
                tmp_path
                / "diagnostic_candidates.matches.vcf.gz"
        ),
    )


def write_annotated_vcf(
        output_path: Path,
) -> Path:
    """
    Write a minimal GeneBe-annotated VCF with two transcript annotations.
    """
    first_annotation = "|".join(
        [
            "MSH2",
            ".",
            "NM_000251.3",
            "missense_variant",
            ".",
            ".",
            ".",
            "PM2&PP3",
            "Likely_pathogenic",
            "NM_000251.3:c.2030C>A",
            "NP_000242.1:p.Pro677His",
        ]
    )

    second_annotation = "|".join(
        [
            "MSH2",
            ".",
            "NM_000251.4",
            "missense_variant",
            ".",
            ".",
            ".",
            "PM2",
            "Uncertain_significance",
            "NM_000251.4:c.2030C>A",
            "NP_000242.2:p.Pro677His",
        ]
    )

    output_path.write_text(
        (
            "##fileformat=VCFv4.2\n"
            '##INFO=<ID=gene_symbol_base,Number=.,Type=String,'
            'Description="Gene symbols">\n'
            '##INFO=<ID=acmg_by_gene_base,Number=.,Type=String,'
            'Description="GeneBe ACMG annotations">\n'
            '##INFO=<ID=dbsnp_base,Number=.,Type=String,'
            'Description="dbSNP identifiers">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "2\t47476391\t.\tC\tA\t.\tPASS\t"
            f"gene_symbol_base=MSH2;"
            f"acmg_by_gene_base={first_annotation},{second_annotation};"
            "dbsnp_base=rs123456\n"
        ),
        encoding="utf-8",
    )

    return output_path


# ---------------------------------------------------------------------------
# Builds a confirmed structured result from one detected candidate.
# ---------------------------------------------------------------------------
def test_build_confirmed_result(
        tmp_path,
):
    """
    Verify the status and patient VCF fields for one detected candidate.
    """
    request = build_request()
    candidate = build_normalized_candidate()

    variant_match = build_variant_match(
        candidate=candidate,
        found=True,
    )

    builder = VariantConfirmationResultBuilder()

    result = builder.build(
        request=request,
        candidates=[
            candidate
        ],
        matching_output=build_matching_output(
            matches=[
                variant_match
            ],
            tmp_path=tmp_path,
        ),
    )

    assert result.get_status() == "confirmed"
    assert result.is_found() is True
    assert result.is_ambiguous() is False

    result_data = result.to_dict()

    assert result_data["candidate_count"] == 1
    assert result_data["match_count"] == 1
    assert result_data["candidates"][0]["candidate_id"] == "candidate_1"
    assert result_data["matches"][0]["candidate_id"] == "candidate_1"
    assert result_data["matches"][0]["genotype"] == "0/1"
    assert result_data["matches"][0]["quality"] == 99.0
    assert result_data["matches"][0]["filters"] == [
        "PASS"
    ]
    assert result_data["matches"][0]["sample_format"] == {
        "GT": "0/1",
        "DP": 42,
    }


# ---------------------------------------------------------------------------
# Retains a candidate that is absent from the patient VCF.
# ---------------------------------------------------------------------------
def test_build_not_found_result(
        tmp_path,
):
    """
    Verify that a non-detected candidate remains fully traceable.
    """
    request = build_request()
    candidate = build_normalized_candidate()

    variant_match = build_variant_match(
        candidate=candidate,
        found=False,
    )

    builder = VariantConfirmationResultBuilder()

    result = builder.build(
        request=request,
        candidates=[
            candidate
        ],
        matching_output=build_matching_output(
            matches=[
                variant_match
            ],
            tmp_path=tmp_path,
        ),
    )

    result_data = result.to_dict()

    assert result_data["status"] == "not_found"
    assert result_data["found"] is False
    assert result_data["candidate_count"] == 1
    assert result_data["match_count"] == 0

    assert result_data["matches"][0]["found"] is False
    assert result_data["matches"][0]["genotype"] is None
    assert result_data["matches"][0]["quality"] is None
    assert result_data["matches"][0]["annotations"] == []

    assert result_data["warnings"] == [
        builder.NOT_FOUND_WARNING
    ]


# ---------------------------------------------------------------------------
# Represents one input resolving into multiple genomic candidates.
# ---------------------------------------------------------------------------
def test_build_ambiguous_result(
        tmp_path,
):
    """
    Verify candidate order, ambiguity status and ambiguity warning.
    """
    request = build_request(
        variant="NP_000242.1:p.Pro677His",
        representation_type="hgvsp",
    )

    first_candidate = build_normalized_candidate(
        candidate_id="candidate_1",
        position=47476391,
        alternate="A",
    )

    second_candidate = build_normalized_candidate(
        candidate_id="candidate_2",
        position=47476391,
        alternate="G",
    )

    first_match = build_variant_match(
        candidate=first_candidate,
        found=True,
    )

    second_match = build_variant_match(
        candidate=second_candidate,
        found=False,
    )

    builder = VariantConfirmationResultBuilder()

    result = builder.build(
        request=request,
        candidates=[
            first_candidate,
            second_candidate,
        ],
        matching_output=build_matching_output(
            matches=[
                second_match,
                first_match,
            ],
            tmp_path=tmp_path,
        ),
    )

    result_data = result.to_dict()

    assert result_data["status"] == "confirmed_ambiguous"
    assert result_data["ambiguous"] is True

    assert [
               candidate["candidate_id"]
               for candidate in result_data["candidates"]
           ] == [
               "candidate_1",
               "candidate_2",
           ]

    assert [
               variant_match["candidate_id"]
               for variant_match in result_data["matches"]
           ] == [
               "candidate_1",
               "candidate_2",
           ]

    assert result_data["warnings"].count(
        builder.AMBIGUOUS_CONVERSION_WARNING
    ) == 1


# ---------------------------------------------------------------------------
# Attaches GeneBe annotations to the matching detected variant.
# ---------------------------------------------------------------------------
def test_attach_genebe_annotations(
        tmp_path,
):
    """
    Verify parsing and association of multiple GeneBe transcript annotations.
    """
    request = build_request()
    candidate = build_normalized_candidate()

    variant_match = build_variant_match(
        candidate=candidate,
        found=True,
    )

    annotated_vcf_path = write_annotated_vcf(
        tmp_path
        / "diagnostic_candidates.genebe.vcf"
    )

    builder = VariantConfirmationResultBuilder()

    result = builder.build(
        request=request,
        candidates=[
            candidate
        ],
        matching_output=build_matching_output(
            matches=[
                variant_match
            ],
            tmp_path=tmp_path,
        ),
        annotated_vcf_path=annotated_vcf_path,
    )

    annotations = result.to_dict()["matches"][0]["annotations"]

    assert len(annotations) == 2

    first_annotation = annotations[0]

    assert first_annotation["gene"] == "MSH2"
    assert first_annotation["transcript"] == "NM_000251.3"
    assert first_annotation["consequence"] == "missense_variant"
    assert first_annotation["hgvsc"] == "NM_000251.3:c.2030C>A"
    assert first_annotation["hgvsp"] == "NP_000242.1:p.Pro677His"
    assert first_annotation["dbsnp"] == "rs123456"
    assert first_annotation["acmg_classification"] == (
        "Likely_pathogenic"
    )
    assert first_annotation["acmg_criteria"] == [
        "PM2",
        "PP3",
    ]

    assert annotations[1]["transcript"] == "NM_000251.4"
    assert builder.MISSING_ANNOTATION_WARNING not in result.warnings


# ---------------------------------------------------------------------------
# Rejects results when a candidate has no corresponding VariantMatch.
# ---------------------------------------------------------------------------
def test_reject_candidate_without_match(
        tmp_path,
):
    """
    Verify that every candidate must have exactly one matching result.
    """
    first_candidate = build_normalized_candidate(
        candidate_id="candidate_1",
        position=100,
    )

    second_candidate = build_normalized_candidate(
        candidate_id="candidate_2",
        position=200,
    )

    first_match = build_variant_match(
        candidate=first_candidate,
        found=True,
    )

    builder = VariantConfirmationResultBuilder()

    with pytest.raises(
            InvalidResultBuilderInputError,
            match="candidate_2",
    ):
        builder.build(
            request=build_request(),
            candidates=[
                first_candidate,
                second_candidate,
            ],
            matching_output=build_matching_output(
                matches=[
                    first_match
                ],
                tmp_path=tmp_path,
            ),
        )


# ---------------------------------------------------------------------------
# Writes the complete structured result to the standard JSON output.
# ---------------------------------------------------------------------------
def test_write_structured_result_json(
        tmp_path,
):
    """
    Verify atomic JSON creation and consistency with the in-memory result.
    """
    request = build_request()
    candidate = build_normalized_candidate()

    variant_match = build_variant_match(
        candidate=candidate,
        found=True,
    )

    builder = VariantConfirmationResultBuilder()

    output_directory = (
            tmp_path
            / "variant_confirmation"
    )

    output = builder.build_to_directory(
        request=request,
        candidates=[
            candidate
        ],
        matching_output=build_matching_output(
            matches=[
                variant_match
            ],
            tmp_path=tmp_path,
        ),
        output_directory=output_directory,
    )

    expected_path = (
            output_directory
            / "diagnostic_variant_results.json"
    )

    assert output.json_path == expected_path
    assert expected_path.is_file()

    with expected_path.open(
            "r",
            encoding="utf-8",
    ) as input_handle:
        json_data = json.load(
            input_handle
        )

    assert json_data == output.result.to_dict()

    remaining_files = [
        path.name
        for path in output_directory.iterdir()
    ]

    assert remaining_files == [
        "diagnostic_variant_results.json"
    ]