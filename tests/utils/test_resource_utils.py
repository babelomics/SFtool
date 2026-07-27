import pytest

from sftool.utils.resource_utils import (
    ResourceSpecificationError,
    load_bundled_resources,
    validate_bundled_resources,
)


def test_load_bundled_resources() -> None:
    specification = load_bundled_resources()

    assert specification["schema_version"] == 1

    assert "clinvar" in specification
    assert "hpo" in specification
    assert "pharmcat" in specification
    assert "reference_genomes" in specification


def test_bundled_clinvar_definition() -> None:
    specification = load_bundled_resources()
    clinvar = specification["clinvar"]

    assert clinvar["version"]
    assert clinvar["archive_base_url"]
    assert "{version}" in clinvar["variant_summary_filename"]
    assert "{version}" in clinvar["submission_summary_filename"]

def test_bundled_reference_genome_assemblies() -> None:
    specification = load_bundled_resources()
    genomes = specification["reference_genomes"]

    assert "GRCh37" in genomes
    assert "GRCh38" in genomes

def valid_resource_specification() -> dict:
    """
    Return a valid resource specification for validation tests.
    """

    return {
        "schema_version": 1,
        "clinvar": {
            "version": "2026-07",
            "archive_base_url": (
                "https://ftp.ncbi.nlm.nih.gov/"
                "pub/clinvar/tab_delimited/archive/"
            ),
            "variant_summary_filename": (
                "variant_summary_{version}.txt.gz"
            ),
            "submission_summary_filename": (
                "submission_summary_{version}.txt.gz"
            ),
        },
        "hpo": {
            "version": "test-version",
            "url": "https://example.org/genes_to_phenotype.txt",
            "filename": "genes_to_phenotype.txt",
        },
        "pharmcat": {
            "version": "test-version",
            "url": "https://example.org/pharmcat_positions.vcf",
            "filename": "pharmcat_positions.vcf",
        },
        "reference_genomes": {
            "GRCh37": {
                "url": None,
                "filename": None,
            },
            "GRCh38": {
                "url": None,
                "filename": None,
            },
        },
    }

def test_validate_bundled_resources_accepts_valid_specification() -> None:
    specification = valid_resource_specification()

    validate_bundled_resources(specification)

def test_validate_bundled_resources_rejects_missing_section() -> None:
    specification = valid_resource_specification()
    del specification["hpo"]

    with pytest.raises(
            ResourceSpecificationError,
            match="hpo",
    ):
        validate_bundled_resources(specification)