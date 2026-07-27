from pathlib import Path
from unittest.mock import Mock

import pytest

from sftool.utils import catalog_utils
from sftool.utils.catalog_utils import (
    CatalogGenerationError,
    GeneCoordinate,
    build_catalog_resources,
    collect_gene_coordinates,
    copy_rr_str_catalog,
    write_bed_file,
)


def test_collect_gene_coordinates_uses_grch37_alias(monkeypatch):
    lookup = Mock(
        return_value={
            "Chromosome": "6",
            "Start": 100,
            "End": 200,
        }
    )
    monkeypatch.setattr(
        catalog_utils,
        "get_gene_location_ensembl",
        lookup,
    )

    coordinates = collect_gene_coordinates(
        ["MMUT"],
        "GRCh37",
    )

    lookup.assert_called_once_with("MUT", "GRCh37")
    assert coordinates[0].gene_symbol == "MMUT"


def test_write_bed_file_generates_both_chromosome_styles(tmp_path):
    coordinates = [
        GeneCoordinate("1", 10, 20, "GENE1"),
        GeneCoordinate("MT", 30, 40, "GENE2"),
    ]

    plain = tmp_path / "PR.bed"
    prefixed = tmp_path / "PR.chr.bed"

    write_bed_file(coordinates, plain, chr_prefix=False)
    write_bed_file(coordinates, prefixed, chr_prefix=True)

    assert plain.read_text() == (
        "1\t10\t20\tGENE1\n"
        "MT\t30\t40\tGENE2\n"
    )
    assert prefixed.read_text() == (
        "chr1\t10\t20\tGENE1\n"
        "chrM\t30\t40\tGENE2\n"
    )


def test_build_catalog_resources_creates_expected_files(
        tmp_path,
        monkeypatch,
):
    source_csv = tmp_path / "PR.csv"
    source_csv.write_text(
        "Gene,Phenotype\nGENE1,Condition\n",
        encoding="latin1",
    )

    monkeypatch.setattr(
        catalog_utils,
        "collect_gene_coordinates",
        Mock(
            return_value=[
                GeneCoordinate("1", 10, 20, "GENE1")
            ]
        ),
    )

    outputs = build_catalog_resources(
        category="PR",
        assembly="GRCh38",
        source_csv=source_csv,
        output_dir=tmp_path / "catalogs",
    )

    assert outputs["bed"].exists()
    assert outputs["chr_bed"].exists()
    assert outputs["json"].exists()


def test_build_catalog_resources_rejects_rr_str(tmp_path):
    with pytest.raises(CatalogGenerationError):
        build_catalog_resources(
            category="RR_STR",
            assembly="GRCh38",
            source_csv=tmp_path / "RR_STR.csv",
            output_dir=tmp_path / "catalogs",
        )


def test_copy_rr_str_catalog_preserves_file(tmp_path, monkeypatch):
    source = tmp_path / "RR_STR.csv"
    source.write_text("Gene\nFMR1\n", encoding="latin1")

    monkeypatch.setattr(
        catalog_utils,
        "get_bundled_catalog_resource",
        lambda category: source,
    )

    destination = copy_rr_str_catalog(
        tmp_path / "resources"
    )

    assert destination.read_bytes() == source.read_bytes()