import pytest
import hashlib
from io import BytesIO
from pathlib import Path
from unittest.mock import patch
from urllib.error import URLError


from sftool.utils.resource_utils import (
    ResourceOperationError,
    ResourceSpecificationError,
    calculate_sha256,
    download_file,
    ensure_directory,
    read_json,
    render_resource_filename,
    resolve_resource_url,
    write_json,
    load_bundled_resources,
    validate_bundled_resources,
    download_versioned_resource,
    validate_hpo_gene_to_phenotype,
    validate_vcf_resource,
)

HPO_HEADER = (
    "ncbi_gene_id\t"
    "gene_symbol\t"
    "hpo_id\t"
    "hpo_name\t"
    "frequency\t"
    "disease_id"
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
            "installed_filename": "genes_to_phenotype.txt",
        },
        "pharmcat": {
            "version": "test-version",
            "url": "https://example.org/pharmcat_positions.vcf",
            "installed_filename": "pharmcat_positions.vcf",
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

def test_ensure_directory_creates_nested_directory(tmp_path: Path) -> None:
    directory = tmp_path / "resources" / "downloads"

    result = ensure_directory(directory)

    assert result == directory
    assert directory.is_dir()


def test_ensure_directory_rejects_existing_file(tmp_path: Path) -> None:
    path = tmp_path / "resources"
    path.write_text("content", encoding="utf-8")

    with pytest.raises(ResourceOperationError, match="not a directory"):
        ensure_directory(path)


def test_render_resource_filename() -> None:
    result = render_resource_filename(
        "variant_summary_{version}.txt.gz",
        version="2026-07",
    )

    assert result == "variant_summary_2026-07.txt.gz"


def test_render_resource_filename_rejects_unknown_placeholder() -> None:
    with pytest.raises(ResourceSpecificationError):
        render_resource_filename(
            "variant_summary_{release}.txt.gz",
            version="2026-07",
        )

@pytest.mark.parametrize(
    "base_url",
    [
        "https://example.org/resources",
        "https://example.org/resources/",
    ],
)
def test_resolve_resource_url(base_url: str) -> None:
    result = resolve_resource_url(base_url, "resource.txt.gz")

    assert result == "https://example.org/resources/resource.txt.gz"

def test_calculate_sha256(tmp_path: Path) -> None:
    content = b"SFtool"
    path = tmp_path / "resource.txt"
    path.write_bytes(content)

    result = calculate_sha256(path)

    assert result == hashlib.sha256(content).hexdigest()

def test_write_and_read_json(tmp_path: Path) -> None:
    path = tmp_path / "resources" / "resources.json"
    data = {
        "schema_version": 1,
        "resources": {},
    }

    result = write_json(data, path)

    assert result == path
    assert read_json(path) == data


def test_read_json_rejects_invalid_json(tmp_path: Path) -> None:
    path = tmp_path / "resources.json"
    path.write_text("{invalid", encoding="utf-8")

    with pytest.raises(ResourceOperationError, match="Could not read JSON"):
        read_json(path)

def test_download_file(tmp_path: Path) -> None:
    destination = tmp_path / "downloads" / "resource.txt"

    with patch(
            "sftool.utils.resource_utils.urlopen",
            return_value=BytesIO(b"downloaded content"),
    ):
        result = download_file(
            "https://example.org/resource.txt",
            destination,
        )

    assert result == destination
    assert destination.read_bytes() == b"downloaded content"
    assert not destination.with_name(
        f"{destination.name}.part"
    ).exists()


def test_download_file_reuses_existing_file(tmp_path: Path) -> None:
    destination = tmp_path / "resource.txt"
    destination.write_bytes(b"existing content")

    with patch(
            "sftool.utils.resource_utils.urlopen",
    ) as mocked_urlopen:
        result = download_file(
            "https://example.org/resource.txt",
            destination,
        )

    assert result == destination
    assert destination.read_bytes() == b"existing content"
    mocked_urlopen.assert_not_called()


def test_download_file_cleans_up_after_failure(tmp_path: Path) -> None:
    destination = tmp_path / "resource.txt"
    partial_path = destination.with_name(
        f"{destination.name}.part"
    )

    with patch(
            "sftool.utils.resource_utils.urlopen",
            side_effect=URLError("network error"),
    ):
        with pytest.raises(
                ResourceOperationError,
                match="Could not download resource",
        ):
            download_file(
                "https://example.org/resource.txt",
                destination,
            )

    assert not destination.exists()
    assert not partial_path.exists()

def test_download_versioned_resource_creates_expected_path(
        tmp_path,
        monkeypatch,
):
    specification = {
        "version": "v2026-06-23",
        "url": "https://example.org/genes_to_phenotype.txt",
        "installed_filename": "genes_to_phenotype_{version}.txt",
    }

    def fake_download_file(url, destination, **kwargs):
        assert url == specification["url"]

        destination.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        destination.write_text(
            f"{HPO_HEADER}\n"
            "1\tA1BG\tHP:0000001\tAll\t-\tOMIM:123456\n",
            encoding="utf-8",
        )

        return destination

    monkeypatch.setattr(
        "sftool.utils.resource_utils.download_file",
        fake_download_file,
    )

    result = download_versioned_resource(
        output_root=tmp_path,
        resource_name="hpo",
        specification=specification,
        validator=validate_hpo_gene_to_phenotype,
    )

    expected_path = (
            tmp_path
            / "hpo"
            / "genes_to_phenotype_v2026-06-23.txt"
    )

    assert expected_path.is_file()
    assert result["version"] == "v2026-06-23"
    assert result["source_url"] == specification["url"]
    assert result["path"] == (
        "hpo/genes_to_phenotype_v2026-06-23.txt"
    )
    assert result["sha256"]


def test_validate_hpo_gene_to_phenotype_accepts_valid_file(
        tmp_path,
):
    hpo_path = tmp_path / "genes_to_phenotype.txt"

    hpo_path.write_text(
        "ncbi_gene_id\t"
        "gene_symbol\t"
        "hpo_id\t"
        "hpo_name\t"
        "frequency\t"
        "disease_id\n"
        "1\tA1BG\tHP:0000001\tAll\t-\tOMIM:123456\n",
        encoding="utf-8",
    )

    validate_hpo_gene_to_phenotype(hpo_path)


def test_validate_hpo_gene_to_phenotype_rejects_invalid_header(
        tmp_path,
):
    hpo_path = tmp_path / "genes_to_phenotype.txt"

    hpo_path.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\n"
        "1\tA1BG\tHP:0000001\n",
        encoding="utf-8",
    )

    with pytest.raises(
            ResourceOperationError,
            match="unexpected header",
    ):
        validate_hpo_gene_to_phenotype(hpo_path)


def test_validate_vcf_resource_accepts_valid_vcf(
        tmp_path,
):
    vcf_path = tmp_path / "pharmcat_positions.vcf"

    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t.\tA\tG\t.\tPASS\t.\n",
        encoding="utf-8",
    )

    validate_vcf_resource(vcf_path)