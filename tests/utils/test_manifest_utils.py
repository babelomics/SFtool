from pathlib import Path

import pytest

from sftool.utils.manifest_utils import (
    build_catalog_manifest,
    build_installed_manifest,
    describe_installed_file,
)

from sftool.utils.resource_utils import ResourceOperationError
from sftool.utils.manifest_utils import describe_installed_file

def test_describe_installed_file(tmp_path: Path) -> None:
    output_root = tmp_path / "resources"
    output_root.mkdir()

    resource = output_root / "hpo" / "hp.obo"
    resource.parent.mkdir()
    resource.write_text("test")

    descriptor = describe_installed_file(
        resource,
        output_root=output_root,
    )

    assert descriptor["path"] == "hpo/hp.obo"
    assert "sha256" in descriptor

def test_describe_installed_file_outside_root(tmp_path: Path) -> None:
    output_root = tmp_path / "resources"
    output_root.mkdir()

    resource = tmp_path / "outside.txt"
    resource.write_text("test", encoding="utf-8")

    with pytest.raises(
            ResourceOperationError,
            match="outside the installation root",
    ):
        describe_installed_file(
            resource,
            output_root=output_root,
        )

def test_build_catalog_manifest(tmp_path: Path) -> None:
    output_root = tmp_path / "resources"
    output_root.mkdir()

    bed = output_root / "catalogs" / "PR" / "GRCh38" / "pr.bed"
    json = output_root / "catalogs" / "PR" / "GRCh38" / "pr.json"

    bed.parent.mkdir(parents=True)
    bed.write_text("bed")
    json.write_text("{}")

    manifest = build_catalog_manifest(
        {
            "PR": {
                "GRCh38": {
                    "bed": bed,
                    "json": json,
                }
            }
        },
        output_root=output_root,
    )

    assert manifest["PR"]["GRCh38"]["bed"]["path"] == "catalogs/PR/GRCh38/pr.bed"
    assert manifest["PR"]["GRCh38"]["json"]["path"] == "catalogs/PR/GRCh38/pr.json"


def test_build_installed_manifest(tmp_path: Path) -> None:
    output_root = tmp_path / "resources"
    output_root.mkdir()

    hpo = output_root / "hpo" / "hp.obo"
    pharmcat = output_root / "pharmcat" / "positions.vcf.gz"

    hpo.parent.mkdir()
    pharmcat.parent.mkdir()

    hpo.write_text("ontology")
    pharmcat.write_text("vcf")

    manifest = build_installed_manifest(
        output_root=output_root,
        resource_version="stable",
        catalog_resources={},
        clinvar_resources={"version": "2026-07", "files": {}},
        hpo_resource={
            "version": "2026-06",
            "path": hpo,
        },
        pharmcat_resource={
            "version": "3.0.0",
            "positions_vcf": pharmcat,
        },
    )

    assert manifest["resource_version"] == "stable"
    assert manifest["datasets"]["hpo"]["file"]["path"] == "hpo/hp.obo"
    assert (
            manifest["datasets"]["pharmcat"]["positions_vcf"]["path"]
            == "pharmcat/positions.vcf.gz"
    )