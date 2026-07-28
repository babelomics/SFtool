from pathlib import Path

import pytest
import csv
import gzip
from click.testing import CliRunner

from sftool.utils import clinvar_utils
from sftool.utils.clinvar_utils import (
    download_bundled_clinvar_snapshot,
    build_clinvar_databases,
    CLINVAR_DATABASE_COLUMNS,
    build_clinvar_assembly_database
)
from sftool.utils.resource_utils import (
    ResourceSpecificationError,
)


def test_download_bundled_clinvar_snapshot(
        tmp_path,
        monkeypatch,
):
    downloaded: list[tuple[str, Path, bool]] = []

    def fake_download_file(
            url,
            destination,
            *,
            overwrite=False,
            **kwargs,
    ):
        destination = Path(destination)
        destination.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        destination.write_bytes(b"clinvar-test-data")

        downloaded.append(
            (url, destination, overwrite)
        )

        return destination

    monkeypatch.setattr(
        clinvar_utils,
        "download_file",
        fake_download_file,
    )

    specification = {
        "version": "2026-06",
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
    }

    result = download_bundled_clinvar_snapshot(
        output_root=tmp_path,
        clinvar_specification=specification,
    )

    source_directory = (
            tmp_path / "clinvar" / "source"
    )

    expected_variant_summary = (
            source_directory
            / "variant_summary_2026-06.txt.gz"
    )
    expected_submission_summary = (
            source_directory
            / "submission_summary_2026-06.txt.gz"
    )

    assert expected_variant_summary.exists()
    assert expected_submission_summary.exists()

    assert result["version"] == "2026-06"
    assert result["variant_summary"] == str(
        expected_variant_summary
    )
    assert result["submission_summary"] == str(
        expected_submission_summary
    )

    assert len(downloaded) == 2

    assert downloaded[0][0].endswith(
        "/variant_summary_2026-06.txt.gz"
    )
    assert downloaded[1][0].endswith(
        "/submission_summary_2026-06.txt.gz"
    )


def test_download_bundled_clinvar_snapshot_reuses_files(
        tmp_path,
        monkeypatch,
):
    calls = 0

    def fake_download_file(
            url,
            destination,
            *,
            overwrite=False,
            **kwargs,
    ):
        nonlocal calls
        calls += 1

        destination = Path(destination)
        destination.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        if not destination.exists():
            destination.write_bytes(b"clinvar-test-data")

        return destination

    monkeypatch.setattr(
        clinvar_utils,
        "download_file",
        fake_download_file,
    )

    specification = {
        "version": "2026-06",
        "archive_base_url": "https://example.org/archive/",
        "variant_summary_filename": (
            "variant_summary_{version}.txt.gz"
        ),
        "submission_summary_filename": (
            "submission_summary_{version}.txt.gz"
        ),
    }

    first = download_bundled_clinvar_snapshot(
        output_root=tmp_path,
        clinvar_specification=specification,
    )
    second = download_bundled_clinvar_snapshot(
        output_root=tmp_path,
        clinvar_specification=specification,
    )

    assert first == second
    assert calls == 4


def test_download_bundled_clinvar_snapshot_rejects_missing_field(
        tmp_path,
):
    specification = {
        "version": "2026-06",
        "archive_base_url": "https://example.org/archive/",
        "variant_summary_filename": (
            "variant_summary_{version}.txt.gz"
        ),
    }

    with pytest.raises(
            ResourceSpecificationError,
            match="submission_summary_filename",
    ):
        download_bundled_clinvar_snapshot(
            output_root=tmp_path,
            clinvar_specification=specification,
        )


def test_build_clinvar_assembly_database_filters_assembly(
        tmp_path,
):
    source_path = tmp_path / "variant_summary.txt.gz"

    header = list(CLINVAR_DATABASE_COLUMNS)

    grch37_row = [
        "single nucleotide variant",
        "variant 37",
        "GENE1",
        "Pathogenic",
        "1",
        "123",
        "1001",
        "MONDO:1",
        "Disease 1",
        "GRCh37",
        "1",
        "100",
        "100",
        "criteria provided, single submitter",
        "1",
        "100",
        "A",
        "G",
    ]

    grch38_row = [
        "single nucleotide variant",
        "variant 38",
        "GENE2",
        "Pathogenic",
        "1",
        "456",
        "1002",
        "MONDO:2",
        "Disease 2",
        "GRCh38",
        "2",
        "200",
        "200",
        "reviewed by expert panel",
        "1",
        "200",
        "C",
        "T",
    ]

    with gzip.open(
            source_path,
            "wt",
            encoding="utf-8",
            newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(header)
        writer.writerow(grch37_row)
        writer.writerow(grch38_row)

    output_path = build_clinvar_assembly_database(
        variant_summary_path=source_path,
        output_root=tmp_path / "resources",
        assembly="GRCh37",
        version="2026-06",
    )

    rows = list(
        csv.reader(
            output_path.open(encoding="utf-8"),
            delimiter="\t",
        )
    )

    assert rows == [
        header,
        grch37_row,
    ]

def test_build_clinvar_databases_generates_both_assemblies(
        tmp_path,
        monkeypatch,
):
    calls = []

    def fake_build_clinvar_assembly_database(
            *,
            variant_summary_path,
            output_root,
            assembly,
            version,
            overwrite,
    ):
        calls.append(assembly)

        output_path = (
                output_root
                / "clinvar"
                / assembly
                / f"clinvar_database_{assembly}_{version}.txt"
        )
        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        output_path.write_text(
            "header\n",
            encoding="utf-8",
        )

        return output_path

    monkeypatch.setattr(
        clinvar_utils,
        "build_clinvar_assembly_database",
        fake_build_clinvar_assembly_database,
    )

    resources = build_clinvar_databases(
        variant_summary_path=tmp_path / "variant_summary.txt.gz",
        output_root=tmp_path / "resources",
        version="2026-06",
    )

    assert calls == [
        "GRCh37",
        "GRCh38",
    ]

    assert set(resources) == {
        "GRCh37",
        "GRCh38",
    }