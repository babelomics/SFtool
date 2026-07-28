from pathlib import Path

import pytest
from click.testing import CliRunner

from sftool.utils import clinvar_utils
from sftool.utils.clinvar_utils import (
    download_bundled_clinvar_snapshot,
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

