import gzip

from sftool.variant_confirmation.clinvar_lookup import (
    VariantConfirmationClinVarLookup,
)
from sftool.variant_confirmation.models import VariantMatch


def _write_clinvar_database(path):
    path.write_text(
        "\t".join([
            "Type",
            "Name",
            "GeneSymbol",
            "ClinicalSignificance",
            "ClinSigSimple",
            "RS# (dbSNP)",
            "VariationID",
            "PhenotypeIDS",
            "PhenotypeList",
            "Assembly",
            "Chromosome",
            "Start",
            "Stop",
            "ReviewStatus",
            "SubmitterCategories",
            "PositionVCF",
            "ReferenceAlleleVCF",
            "AlternateAlleleVCF",
        ])
        + "\n"
        + "\t".join([
            "single nucleotide variant",
            "NM_000546.6(TP53):c.215C>G",
            "TP53",
            "Pathogenic",
            "1",
            "1042522",
            "12345",
            "Orphanet:524;OMIM:151623",
            "Li-Fraumeni syndrome",
            "GRCh38",
            "17",
            "7676154",
            "7676154",
            "reviewed by expert panel",
            "3",
            "7676154",
            "G",
            "C",
        ])
        + "\n",
        encoding="utf-8",
        )


def _write_submission_database(path):
    with gzip.open(
            path,
            "wt",
            encoding="utf-8",
    ) as handle:
        handle.write(
            "#VariationID\tClinicalSignificance\t"
            "ContributesToAggregateClassification\n"
        )
        handle.write(
            "12345\tPathogenic\tyes\n"
        )
        handle.write(
            "12345\tLikely pathogenic\tyes\n"
        )
        handle.write(
            "12345\tBenign\tno\n"
        )


def _build_match(
        *,
        candidate_id="candidate_1",
        found=True,
        chromosome="chr17",
        position=7676154,
        reference="G",
        alternate="C",
):
    variant_match = VariantMatch(
        candidate_id=candidate_id,
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
        found=False,
    )

    if found:
        variant_match.set_match_data(
            genotype="0/1",
            quality=99.0,
            filters=["PASS"],
            sample_format={"GT": "0/1"},
        )

    return variant_match


def test_lookup_returns_complete_clinvar_annotation_without_evidence_filter(
        tmp_path,
):
    """
    Verify that an exact detected candidate retrieves the complete ClinVar
    annotation and submission summary without applying clinvar_evidence.
    """
    clinvar_db = tmp_path / "clinvar_database_GRCh38.txt"
    submission_db = tmp_path / "clinvar_submission.txt.gz"

    _write_clinvar_database(clinvar_db)
    _write_submission_database(submission_db)

    annotations = VariantConfirmationClinVarLookup().lookup(
        matches=[_build_match()],
        clinvar_db=clinvar_db,
        clinvar_submission_db=submission_db,
    )

    annotation = annotations["candidate_1"]

    assert annotation["clinical_significance"] == "Pathogenic"
    assert annotation["review_status"] == (
        "(3) reviewed by expert panel"
    )
    assert annotation["stars"] == 3
    assert annotation["clinvar_id"] == "12345"
    assert annotation["rs"] == "rs1042522"
    assert annotation["orpha"] == "524"
    assert annotation["omim"] == "151623"
    assert annotation["clinical_significance_summary"] == (
        "Pathogenic (1); Likely pathogenic (1)"
    )


def test_lookup_matches_coordinates_ignoring_chr_prefix(
        tmp_path,
):
    """
    Verify that candidate and ClinVar chromosomes match when only the candidate
    contains the conventional ``chr`` prefix.
    """
    clinvar_db = tmp_path / "clinvar_database_GRCh38.txt"
    submission_db = tmp_path / "clinvar_submission.txt.gz"

    _write_clinvar_database(clinvar_db)
    _write_submission_database(submission_db)

    annotations = VariantConfirmationClinVarLookup().lookup(
        matches=[
            _build_match(
                chromosome="chr17",
            )
        ],
        clinvar_db=clinvar_db,
        clinvar_submission_db=submission_db,
    )

    assert "candidate_1" in annotations


def test_lookup_ignores_candidates_not_detected_in_patient_vcf(
        tmp_path,
):
    """
    Verify that candidates not found in the patient VCF are not queried or
    annotated with ClinVar information.
    """
    clinvar_db = tmp_path / "clinvar_database_GRCh38.txt"
    submission_db = tmp_path / "clinvar_submission.txt.gz"

    _write_clinvar_database(clinvar_db)
    _write_submission_database(submission_db)

    annotations = VariantConfirmationClinVarLookup().lookup(
        matches=[
            _build_match(
                found=False,
            )
        ],
        clinvar_db=clinvar_db,
        clinvar_submission_db=submission_db,
    )

    assert annotations == {}


def test_lookup_returns_empty_result_when_variant_is_absent_from_clinvar(
        tmp_path,
):
    """
    Verify that a detected candidate absent from the local ClinVar database
    produces no annotation and does not fail the confirmation workflow.
    """
    clinvar_db = tmp_path / "clinvar_database_GRCh38.txt"
    submission_db = tmp_path / "clinvar_submission.txt.gz"

    _write_clinvar_database(clinvar_db)
    _write_submission_database(submission_db)

    annotations = VariantConfirmationClinVarLookup().lookup(
        matches=[
            _build_match(
                position=123456,
            )
        ],
        clinvar_db=clinvar_db,
        clinvar_submission_db=submission_db,
    )

    assert annotations == {}
