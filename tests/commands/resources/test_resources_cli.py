from pathlib import Path

from click.testing import CliRunner

from sftool.cli import main


def test_resources_group_is_registered() -> None:
    runner = CliRunner()

    result = runner.invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "resources" in result.output
    assert "Manage SFtool datasets" in result.output


def test_resources_help() -> None:
    runner = CliRunner()

    result = runner.invoke(
        main,
        ["resources", "--help"],
    )

    assert result.exit_code == 0
    assert "Manage SFtool datasets" in result.output
    assert "setup" in result.output


def test_resources_setup_help() -> None:
    runner = CliRunner()

    result = runner.invoke(
        main,
        ["resources", "setup", "--help"],
    )

    assert result.exit_code == 0
    assert "--output-dir" in result.output
    assert "--clinvar-evidence" in result.output
    assert "--resource-version" in result.output
    assert "--download-reference-genomes" in result.output


def test_setup_requires_output_directory() -> None:
    runner = CliRunner()

    result = runner.invoke(
        main,
        ["resources", "setup"],
    )

    assert result.exit_code != 0
    assert "Missing option '--output-dir'" in result.output


def test_setup_uses_default_options() -> None:
    runner = CliRunner()

    with runner.isolated_filesystem():
        result = runner.invoke(
            main,
            [
                "resources",
                "setup",
                "--output-dir",
                "resources",
            ],
        )

    assert result.exit_code == 0
    assert "ClinVar evidence: 1" in result.output
    assert "Resource version: bundled" in result.output
    assert "Download reference genomes: no" in result.output


def test_setup_accepts_explicit_options() -> None:
    runner = CliRunner()

    with runner.isolated_filesystem():
        result = runner.invoke(
            main,
            [
                "resources",
                "setup",
                "--output-dir",
                "resources",
                "--clinvar-evidence",
                "3",
                "--resource-version",
                "latest",
                "--download-reference-genomes",
            ],
        )

    assert result.exit_code == 0
    assert "ClinVar evidence: 3" in result.output
    assert "Resource version: latest" in result.output
    assert "Download reference genomes: yes" in result.output

def test_setup_rejects_invalid_clinvar_evidence() -> None:
    runner = CliRunner()

    result = runner.invoke(
        main,
        [
            "resources",
            "setup",
            "--output-dir",
            "resources",
            "--clinvar-evidence",
            "6",
        ],
    )

    assert result.exit_code != 0
    assert "Invalid value for '--clinvar-evidence'" in result.output

def test_setup_rejects_invalid_resource_version() -> None:
    runner = CliRunner()

    result = runner.invoke(
        main,
        [
            "resources",
            "setup",
            "--output-dir",
            "resources",
            "--resource-version",
            "unknown",
        ],
    )

    assert result.exit_code != 0
    assert "Invalid value for '--resource-version'" in result.output