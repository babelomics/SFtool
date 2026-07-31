from pathlib import Path

import pytest

from sftool.core.resources import RuntimeResources
from sftool.utils.errors import ResourceManifestError


@pytest.fixture
def runtime_resources(tmp_path: Path) -> RuntimeResources:
    manifest_path = tmp_path / "resources.json"
    resources_root = tmp_path / "resources"

    execution = {
        "assembly": "GRCh38",
        "reference_genome": {
            "fasta": resources_root / "references/GRCh38/genome.fa",
            "fai": resources_root / "references/GRCh38/genome.fa.fai",
            "source": "manifest",
            "fasta_sha256": "reference-fasta-checksum",
            "fai_sha256": "reference-index-checksum",
        },
        "catalogs": {
            "PR": {
                "version": "3.1",
                "bed": {
                    "path": resources_root / "catalogs/GRCh38/PR.bed",
                    "sha256": "pr-bed-checksum",
                },
                "chr_bed": {
                    "path": resources_root / "catalogs/GRCh38/PR.chr.bed",
                    "sha256": "pr-chr-bed-checksum",
                },
                "json": {
                    "path": resources_root / "catalogs/GRCh38/PR.json",
                    "sha256": "pr-json-checksum",
                },
            },
        },
        "clinvar": {
            "PR": {
                "path": (
                        resources_root
                        / "clinvar/GRCh38/PR/evidence_1.json"
                ),
                "sha256": "clinvar-pr-checksum",
            },
        },
        "clinvar_database": {
            "path": (
                    resources_root
                    / "clinvar/GRCh38/clinvar_database.json"
            ),
            "sha256": "clinvar-database-checksum",
        },
        "rr_str": None,
        "hpo": {
            "version": "v2026-06-23",
            "file": {
                "path": (
                        resources_root
                        / "hpo/genes_to_phenotype.txt"
                ),
                "sha256": "hpo-checksum",
            },
        },
        "pharmcat": {
            "version": "3.0.1",
            "positions_vcf": {
                "path": (
                        resources_root
                        / "pharmcat/pharmcat_positions.vcf"
                ),
                "sha256": "pharmcat-checksum",
            },
        },
    }

    installed = {
        "root": resources_root,
        "clinvar": {
            "version": "2026-06",
        },
    }

    return RuntimeResources(
        manifest_path=manifest_path,
        manifest={"schema_version": 2},
        installed=installed,
        execution=execution,
    )


def test_exposes_resources_root(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.root == tmp_path / "resources"


def test_exposes_assembly(
        runtime_resources: RuntimeResources,
) -> None:
    assert runtime_resources.assembly == "GRCh38"


def test_exposes_reference_genome(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.reference_genome == (
            tmp_path
            / "resources/references/GRCh38/genome.fa"
    )


def test_exposes_reference_genome_index(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.reference_genome_index == (
            tmp_path
            / "resources/references/GRCh38/genome.fa.fai"
    )


def test_exposes_reference_genome_source(
        runtime_resources: RuntimeResources,
) -> None:
    assert runtime_resources.reference_genome_source == "manifest"


def test_exposes_reference_checksums(
        runtime_resources: RuntimeResources,
) -> None:
    assert (
            runtime_resources.reference_genome_sha256
            == "reference-fasta-checksum"
    )
    assert (
            runtime_resources.reference_genome_index_sha256
            == "reference-index-checksum"
    )


def test_exposes_catalog_json(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.catalog_json("PR") == (
            tmp_path
            / "resources/catalogs/GRCh38/PR.json"
    )


def test_exposes_catalog_bed(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.catalog_bed("PR") == (
            tmp_path
            / "resources/catalogs/GRCh38/PR.bed"
    )


def test_exposes_catalog_chr_bed(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.catalog_chr_bed("PR") == (
            tmp_path
            / "resources/catalogs/GRCh38/PR.chr.bed"
    )


def test_exposes_filtered_clinvar_file(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.clinvar_file("PR") == (
            tmp_path
            / "resources/clinvar/GRCh38/PR/evidence_1.json"
    )


def test_exposes_assembly_clinvar_database(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.clinvar_database == (
            tmp_path
            / "resources/clinvar/GRCh38/clinvar_database.json"
    )


def test_exposes_clinvar_version(
        runtime_resources: RuntimeResources,
) -> None:
    assert runtime_resources.clinvar_version == "2026-06"


def test_exposes_hpo_resource(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.hpo_file == (
            tmp_path
            / "resources/hpo/genes_to_phenotype.txt"
    )
    assert runtime_resources.hpo_version == "v2026-06-23"


def test_exposes_pharmcat_resource(
        runtime_resources: RuntimeResources,
        tmp_path: Path,
) -> None:
    assert runtime_resources.pharmcat_positions_vcf == (
            tmp_path
            / "resources/pharmcat/pharmcat_positions.vcf"
    )
    assert runtime_resources.pharmcat_version == "3.0.1"


def test_rr_str_is_none_when_not_selected(
        runtime_resources: RuntimeResources,
) -> None:
    assert runtime_resources.rr_str_catalog is None
    assert runtime_resources.rr_str_version is None


def test_rejects_unselected_catalog(
        runtime_resources: RuntimeResources,
) -> None:
    with pytest.raises(
            ResourceManifestError,
            match=(
                    "No catalog resource was selected "
                    "for category RR"
            ),
    ):
        runtime_resources.catalog_json("RR")


def test_rejects_unselected_clinvar_resource(
        runtime_resources: RuntimeResources,
) -> None:
    with pytest.raises(
            ResourceManifestError,
            match=(
                    "No ClinVar resource was selected "
                    "for category RR"
            ),
    ):
        runtime_resources.clinvar_file("RR")