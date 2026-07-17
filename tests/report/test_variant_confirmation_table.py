from types import SimpleNamespace

import pytest

from sftool.report.variant_confirmation_table import (
    build_variant_confirmation_table,
)
from sftool.variant_confirmation.models import (
    VariantCandidate,
    VariantConfirmationRequest,
    VariantConfirmationResult,
    VariantMatch,
)


def _build_ctx(gene_to_phenotype_file):
    return SimpleNamespace(
        config=SimpleNamespace(
            references=SimpleNamespace(
                gene_to_phenotype_file=str(gene_to_phenotype_file)
            )
        )
    )


def _build_sample(result, hpo_terms=None):
    return SimpleNamespace(
        results={"variant_confirmation": result},
        hpo_terms=hpo_terms or [],
    )


def _build_result(
        *,
        found=True,
        annotations=None,
        genotype="0/1",
        sample_format=None,
):
    request = VariantConfirmationRequest(
        variant="NM_000546.6:c.215C>G",
        representation_type="hgvsc",
    )

    candidate = VariantCandidate(
        candidate_id="candidate_1",
        chromosome="17",
        position=7676154,
        reference="G",
        alternate="C",
        assembly="GRCh38",
    )
    candidate.set_normalized_coordinates(
        chromosome="17",
        position=7676154,
        reference="G",
        alternate="C",
    )

    variant_match = VariantMatch(
        candidate_id="candidate_1",
        chromosome="17",
        position=7676154,
        reference="G",
        alternate="C",
        found=False,
    )

    if found:
        variant_match.set_match_data(
            genotype=genotype,
            quality=99.0,
            filters=["PASS"],
            sample_format=sample_format or {
                "GT": genotype,
                "AD": [46, 38],
                "DP": 84,
                "GQ": 99,
            },
        )

    for annotation in annotations or []:
        variant_match.add_annotation(annotation)

    result = VariantConfirmationResult(request)
    result.add_candidate(candidate)
    result.add_match(variant_match)

    return result


def test_build_variant_confirmation_table_includes_annotation_format_and_hpo(
        tmp_path,
):
    """
    Verify that a detected diagnostic variant with GeneBe annotations is
    correctly converted into a single report row including:

    - input representation
    - GeneBe annotation fields
    - ClinVar information
    - genotype and zygosity
    - FORMAT serialization
    - related sample HPO terms
    """

    gene_hpo_file = tmp_path / "genes_to_phenotype.txt"
    gene_hpo_file.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\n"
        "7157\tTP53\tHP:0002664\n"
        "7157\tTP53\tHP:0003002\n",
        encoding="utf-8",
    )

    result = _build_result(
        annotations=[
            {
                "gene": "TP53",
                "transcript": "NM_000546.6",
                "consequence": "missense_variant",
                "hgvsc": "NM_000546.6:c.215C>G",
                "hgvsp": "NP_000537.3:p.Pro72Arg",
                "dbsnp": "rs1042522",
                "acmg_classification": "Pathogenic",
                "acmg_criteria": ["PS3", "PM2"],
                "clinvar": {
                    "clinical_significance": "Pathogenic",
                    "review_status": "reviewed_by_expert_panel",
                    "variation_id": "12345",
                },
            }
        ]
    )

    table = build_variant_confirmation_table(
        _build_ctx(gene_hpo_file),
        _build_sample(
            result,
            hpo_terms=["HP:0002664", "HP:9999999"],
        ),
    )

    assert table.tab_name == "Variant confirmation"
    assert len(table.rows) == 1

    row = table.rows[0]

    assert row["Input_variant"] == "NM_000546.6:c.215C>G"
    assert row["Input_format"] == "hgvsc"
    assert row["Found_in_sample"] == "Yes"
    assert row["Gene"] == "TP53"
    assert row["Genotype"] == "0/1"
    assert row["Zygosity"] == "HET"
    assert row["related_HPOs_for_sample"] == "HP:0002664"
    assert row["VCF_Sample_FORMAT"] == (
        "GT=0/1; AD=46,38; DP=84; GQ=99"
    )


def test_build_variant_confirmation_table_keeps_not_found_candidate(
        tmp_path,
):
    """
    Verify that diagnostic candidates not detected in the patient VCF are still
    included in the report with sample-dependent fields left empty.
    """

    gene_hpo_file = tmp_path / "genes_to_phenotype.txt"
    gene_hpo_file.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\n",
        encoding="utf-8",
    )

    table = build_variant_confirmation_table(
        _build_ctx(gene_hpo_file),
        _build_sample(
            _build_result(found=False)
        ),
    )

    assert len(table.rows) == 1

    row = table.rows[0]

    assert row["Found_in_sample"] == "No"
    assert row["Gene"] == ""
    assert row["Genotype"] == ""
    assert row["Zygosity"] == ""
    assert row["VCF_Sample_FORMAT"] == ""
    assert row["related_HPOs_for_sample"] == "NA"


def test_build_variant_confirmation_table_creates_one_row_per_annotation(
        tmp_path,
):

    """
    Verify that one report row is generated for every GeneBe annotation
    associated with the same diagnostic candidate.
    """
    gene_hpo_file = tmp_path / "genes_to_phenotype.txt"
    gene_hpo_file.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\n"
        "1\tGENE1\tHP:0000001\n"
        "2\tGENE2\tHP:0000002\n",
        encoding="utf-8",
    )

    result = _build_result(
        annotations=[
            {"gene": "GENE1", "clinvar": {}},
            {"gene": "GENE2", "clinvar": {}},
        ]
    )

    table = build_variant_confirmation_table(
        _build_ctx(gene_hpo_file),
        _build_sample(
            result,
            hpo_terms=["HP:0000001", "HP:0000002"],
        ),
    )

    assert [row["Gene"] for row in table.rows] == [
        "GENE1",
        "GENE2",
    ]
    assert [
               row["related_HPOs_for_sample"]
               for row in table.rows
           ] == [
               "HP:0000001",
               "HP:0000002",
           ]


def test_build_variant_confirmation_table_returns_none_without_result(
        tmp_path,
):
    gene_hpo_file = tmp_path / "genes_to_phenotype.txt"
    gene_hpo_file.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\n",
        encoding="utf-8",
    )

    sample = SimpleNamespace(
        results={},
        hpo_terms=[],
    )

    assert build_variant_confirmation_table(
        _build_ctx(gene_hpo_file),
        sample,
    ) is None
