# tests/variant_confirmation/test_candidate_vcf_normalizer.py

import gzip
import os
import shutil
import subprocess

from pathlib import Path

import pytest

from sftool.variant_confirmation.candidate_vcf import (
    CandidateVcfCommandError,
    CandidateVcfNormalizer,
    CandidateVcfReferenceError,
    CandidateVcfResultError,
    CandidateVcfWriter,
    InvalidCandidateCollectionError,
)
from sftool.variant_confirmation.models import (
    VariantCandidate,
)


def build_candidate(
        candidate_id: str = "candidate_1",
        chromosome: str = "2",
        position: int = 100,
        reference: str = "A",
        alternate: str = "G",
        assembly: str = "GRCh38",
) -> VariantCandidate:
    """
    Build a diagnostic variant candidate for normalizer unit tests.
    """
    return VariantCandidate(
        candidate_id=candidate_id,
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
        assembly=assembly,
    )


def write_raw_vcf(
        output_directory: Path,
        candidates: list[VariantCandidate],
) -> Path:
    """
    Generate the raw candidate VCF consumed by CandidateVcfNormalizer.
    """
    raw_vcf_path = (
            output_directory
            / "diagnostic_candidates.raw.vcf"
    )

    CandidateVcfWriter().write(
        candidates=candidates,
        output_path=raw_vcf_path,
    )

    return raw_vcf_path


class FakeBcftoolsRunner:
    """
    Simulate the filesystem effects of bcftools norm, sort and index.

    The runner records all commands and creates the files expected by
    CandidateVcfNormalizer without invoking a real bcftools executable.
    """

    def __init__(
            self,
            normalized_records: list[str],
            fail_on: str | None = None,
    ):
        self.normalized_records = normalized_records
        self.fail_on = fail_on
        self.commands: list[list[str]] = []
        self.reference_compatible_vcf: str | None = None

    def __call__(
            self,
            command,
            check,
            capture_output,
            text,
    ):
        self.commands.append(
            list(command)
        )

        subcommand = command[1]

        if subcommand == self.fail_on:
            raise subprocess.CalledProcessError(
                returncode=1,
                cmd=command,
                stderr=f"simulated {subcommand} failure",
            )

        if subcommand == "norm":
            input_path = Path(
                command[-1]
            )

            self.reference_compatible_vcf = input_path.read_text(
                encoding="utf-8"
            )

            output_path = Path(
                command[
                    command.index("--output") + 1
                    ]
            )

            self._write_normalized_vcf(
                output_path
            )

        elif subcommand == "sort":
            input_path = Path(
                command[-1]
            )

            output_path = Path(
                command[
                    command.index("--output-file") + 1
                    ]
            )

            shutil.copyfile(
                input_path,
                output_path,
            )

        elif subcommand == "index":
            vcf_path = Path(
                command[-1]
            )

            Path(
                f"{vcf_path}.tbi"
            ).write_bytes(
                b"fake tabix index"
            )

        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout="",
            stderr="",
        )

    def _write_normalized_vcf(
            self,
            output_path: Path,
    ) -> None:
        with gzip.open(
                output_path,
                "wt",
                encoding="utf-8",
        ) as output:
            output.write(
                "##fileformat=VCFv4.2\n"
            )

            output.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            )

            for record in self.normalized_records:
                output.write(
                    record
                )

                output.write(
                    "\n"
                )


@pytest.fixture
def fake_bcftools(
        tmp_path,
) -> Path:
    """
    Create an executable placeholder accepted by the normalizer constructor.
    """
    executable_path = (
            tmp_path
            / "bcftools"
    )

    executable_path.write_text(
        "#!/bin/sh\nexit 0\n",
        encoding="utf-8",
    )

    executable_path.chmod(
        executable_path.stat().st_mode
        | 0o111
    )

    return executable_path


@pytest.fixture
def reference_fasta(
        tmp_path,
) -> Path:
    """
    Create a minimal indexed reference containing chr1, chr2 and chrM.
    """
    fasta_path = (
            tmp_path
            / "reference.fa"
    )

    fasta_path.write_text(
        ">chr1\n"
        "AAAAAAAAAAAAAAAAAAAA\n"
        ">chr2\n"
        "CCCCCCCCCCCCCCCCCCCC\n"
        ">chrM\n"
        "GGGGGGGGGGGGGGGGGGGG\n",
        encoding="utf-8",
    )

    fasta_index_path = Path(
        f"{fasta_path}.fai"
    )

    fasta_index_path.write_text(
        "chr1\t20\t6\t20\t21\n"
        "chr2\t20\t33\t20\t21\n"
        "chrM\t20\t60\t20\t21\n",
        encoding="utf-8",
    )

    return fasta_path


# ---------------------------------------------------------------------------
# Normalizes a candidate, creates the VCF and index, and updates the candidate.
# ---------------------------------------------------------------------------
def test_normalize_candidate_vcf(
        tmp_path,
        fake_bcftools,
        reference_fasta,
):
    candidate = build_candidate(
        chromosome="2",
        position=101,
        reference="AA",
        alternate="A",
    )

    raw_vcf_path = write_raw_vcf(
        output_directory=tmp_path,
        candidates=[candidate],
    )

    output_path = (
            tmp_path
            / "diagnostic_candidates.normalized.vcf.gz"
    )

    runner = FakeBcftoolsRunner(
        normalized_records=[
            (
                "chr2\t100\tcandidate_1\tAA\tA\t.\t.\t"
                "ASSEMBLY=GRCh38"
            ),
        ]
    )

    normalizer = CandidateVcfNormalizer(
        bcftools_path=fake_bcftools,
        runner=runner,
    )

    result = normalizer.normalize(
        raw_vcf_path=raw_vcf_path,
        candidates=[candidate],
        reference_fasta_path=reference_fasta,
        output_path=output_path,
    )

    assert result == output_path

    assert output_path.exists()

    assert Path(
        f"{output_path}.tbi"
    ).exists()

    assert [
               command[1]
               for command in runner.commands
           ] == [
               "norm",
               "sort",
               "index",
           ]

    assert candidate.chromosome == "2"
    assert candidate.position == 101
    assert candidate.reference == "AA"
    assert candidate.alternate == "A"

    assert candidate.is_normalized() is True

    assert candidate.get_normalized_variant() == (
        "chr2:100:AA:A"
    )

    assert candidate.get_matching_key() == (
        "chr2",
        100,
        "AA",
        "A",
    )


# ---------------------------------------------------------------------------
# Harmonizes chromosome names with the contigs declared by the reference FASTA.
# ---------------------------------------------------------------------------
def test_normalize_harmonizes_chromosome_name(
        tmp_path,
        fake_bcftools,
        reference_fasta,
):
    candidate = build_candidate(
        chromosome="2",
    )

    raw_vcf_path = write_raw_vcf(
        output_directory=tmp_path,
        candidates=[candidate],
    )

    runner = FakeBcftoolsRunner(
        normalized_records=[
            "chr2\t100\tcandidate_1\tA\tG\t.\t.\t.",
        ]
    )

    CandidateVcfNormalizer(
        bcftools_path=fake_bcftools,
        runner=runner,
    ).normalize(
        raw_vcf_path=raw_vcf_path,
        candidates=[candidate],
        reference_fasta_path=reference_fasta,
        output_path=tmp_path / "normalized.vcf.gz",
    )

    assert runner.reference_compatible_vcf is not None

    assert "##contig=<ID=2>" not in (
        runner.reference_compatible_vcf
    )

    assert "##contig=<ID=chr2>" in (
        runner.reference_compatible_vcf
    )

    assert (
            "chr2\t100\tcandidate_1\tA\tG\t.\t.\t"
            in runner.reference_compatible_vcf
    )

    assert candidate.normalized_chromosome == "chr2"


# ---------------------------------------------------------------------------
# Updates multiple candidates by candidate ID rather than by record order.
# ---------------------------------------------------------------------------
def test_normalize_multiple_candidates_by_id(
        tmp_path,
        fake_bcftools,
        reference_fasta,
):
    candidate_1 = build_candidate(
        candidate_id="candidate_1",
        chromosome="1",
        position=100,
        reference="A",
        alternate="G",
    )

    candidate_2 = build_candidate(
        candidate_id="candidate_2",
        chromosome="2",
        position=200,
        reference="C",
        alternate="T",
    )

    candidates = [
        candidate_1,
        candidate_2,
    ]

    raw_vcf_path = write_raw_vcf(
        output_directory=tmp_path,
        candidates=candidates,
    )

    runner = FakeBcftoolsRunner(
        normalized_records=[
            (
                "chr2\t200\tcandidate_2\tC\tT\t.\t.\t."
            ),
            (
                "chr1\t100\tcandidate_1\tA\tG\t.\t.\t."
            ),
        ]
    )

    CandidateVcfNormalizer(
        bcftools_path=fake_bcftools,
        runner=runner,
    ).normalize(
        raw_vcf_path=raw_vcf_path,
        candidates=candidates,
        reference_fasta_path=reference_fasta,
        output_path=tmp_path / "normalized.vcf.gz",
    )

    assert candidate_1.get_matching_key() == (
        "chr1",
        100,
        "A",
        "G",
    )

    assert candidate_2.get_matching_key() == (
        "chr2",
        200,
        "C",
        "T",
    )


# ---------------------------------------------------------------------------
# Uses the standard output filename when normalizing into a directory.
# ---------------------------------------------------------------------------
def test_normalize_to_directory_uses_standard_filename(
        tmp_path,
        fake_bcftools,
        reference_fasta,
):
    candidate = build_candidate(
        chromosome="2",
    )

    raw_vcf_path = write_raw_vcf(
        output_directory=tmp_path,
        candidates=[candidate],
    )

    runner = FakeBcftoolsRunner(
        normalized_records=[
            "chr2\t100\tcandidate_1\tA\tG\t.\t.\t.",
        ]
    )

    output_directory = (
            tmp_path
            / "sample_1"
            / "variant_confirmation"
    )

    result = CandidateVcfNormalizer(
        bcftools_path=fake_bcftools,
        runner=runner,
    ).normalize_to_directory(
        raw_vcf_path=raw_vcf_path,
        candidates=[candidate],
        reference_fasta_path=reference_fasta,
        output_directory=output_directory,
    )

    expected_path = (
            output_directory
            / "diagnostic_candidates.normalized.vcf.gz"
    )

    assert result == expected_path

    assert expected_path.exists()

    assert Path(
        f"{expected_path}.tbi"
    ).exists()


# ---------------------------------------------------------------------------
# Rejects a reference FASTA when its required .fai index is missing.
# ---------------------------------------------------------------------------
def test_normalize_requires_reference_fasta_index(
        tmp_path,
        fake_bcftools,
):
    candidate = build_candidate()

    raw_vcf_path = write_raw_vcf(
        output_directory=tmp_path,
        candidates=[candidate],
    )

    fasta_path = (
            tmp_path
            / "reference.fa"
    )

    fasta_path.write_text(
        ">chr2\n"
        "AAAAAAAAAAAAAAAAAAAA\n",
        encoding="utf-8",
    )

    with pytest.raises(
            CandidateVcfReferenceError,
            match="Reference FASTA index not found",
    ):
        CandidateVcfNormalizer(
            bcftools_path=fake_bcftools,
        ).normalize(
            raw_vcf_path=raw_vcf_path,
            candidates=[candidate],
            reference_fasta_path=fasta_path,
            output_path=tmp_path / "normalized.vcf.gz",
        )


# ---------------------------------------------------------------------------
# Rejects candidates whose chromosome cannot be mapped to the reference FASTA.
# ---------------------------------------------------------------------------
def test_normalize_rejects_chromosome_not_in_reference(
        tmp_path,
        fake_bcftools,
        reference_fasta,
):
    candidate = build_candidate(
        chromosome="22",
    )

    raw_vcf_path = write_raw_vcf(
        output_directory=tmp_path,
        candidates=[candidate],
    )

    with pytest.raises(
            CandidateVcfReferenceError,
            match="is not present in the reference FASTA",
    ):
        CandidateVcfNormalizer(
            bcftools_path=fake_bcftools,
        ).normalize(
            raw_vcf_path=raw_vcf_path,
            candidates=[candidate],
            reference_fasta_path=reference_fasta,
            output_path=tmp_path / "normalized.vcf.gz",
        )


# ---------------------------------------------------------------------------
# Rejects duplicate candidate IDs before executing any bcftools command.
# ---------------------------------------------------------------------------
def test_normalize_rejects_duplicate_candidate_ids(
        tmp_path,
        fake_bcftools,
        reference_fasta,
):
    candidates = [
        build_candidate(
            candidate_id="candidate_1",
            alternate="G",
        ),
        build_candidate(
            candidate_id="candidate_1",
            alternate="T",
        ),
    ]

    raw_vcf_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    raw_vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )

    runner = FakeBcftoolsRunner(
        normalized_records=[]
    )

    with pytest.raises(
            InvalidCandidateCollectionError,
            match="Candidate IDs must be unique",
    ):
        CandidateVcfNormalizer(
            bcftools_path=fake_bcftools,
            runner=runner,
        ).normalize(
            raw_vcf_path=raw_vcf_path,
            candidates=candidates,
            reference_fasta_path=reference_fasta,
            output_path=tmp_path / "normalized.vcf.gz",
        )

    assert runner.commands == []


# ---------------------------------------------------------------------------
# Rejects normalized records whose candidate ID is unknown to the normalizer.
# ---------------------------------------------------------------------------
def test_normalize_rejects_unknown_candidate_id(
        tmp_path,
        fake_bcftools,
        reference_fasta,
):
    candidate = build_candidate(
        chromosome="2",
    )

    raw_vcf_path = write_raw_vcf(
        output_directory=tmp_path,
        candidates=[candidate],
    )

    runner = FakeBcftoolsRunner(
        normalized_records=[
            (
                "chr2\t100\tunknown_candidate\tA\tG\t"
                ".\t.\t."
            ),
        ]
    )

    output_path = (
            tmp_path
            / "normalized.vcf.gz"
    )

    with pytest.raises(
            CandidateVcfResultError,
            match="unknown candidate ID",
    ):
        CandidateVcfNormalizer(
            bcftools_path=fake_bcftools,
            runner=runner,
        ).normalize(
            raw_vcf_path=raw_vcf_path,
            candidates=[candidate],
            reference_fasta_path=reference_fasta,
            output_path=output_path,
        )

    assert candidate.is_normalized() is False
    assert not output_path.exists()
    assert not Path(f"{output_path}.tbi").exists()


# ---------------------------------------------------------------------------
# Wraps bcftools failures and removes every temporary normalization file.
# ---------------------------------------------------------------------------
def test_normalize_wraps_command_error_and_cleans_temporary_files(
        tmp_path,
        fake_bcftools,
        reference_fasta,
):
    candidate = build_candidate(
        chromosome="2",
    )

    raw_vcf_path = write_raw_vcf(
        output_directory=tmp_path,
        candidates=[candidate],
    )

    runner = FakeBcftoolsRunner(
        normalized_records=[
            "chr2\t100\tcandidate_1\tA\tG\t.\t.\t.",
        ],
        fail_on="sort",
    )

    output_path = (
            tmp_path
            / "diagnostic_candidates.normalized.vcf.gz"
    )

    with pytest.raises(
            CandidateVcfCommandError,
            match="simulated sort failure",
    ):
        CandidateVcfNormalizer(
            bcftools_path=fake_bcftools,
            runner=runner,
        ).normalize(
            raw_vcf_path=raw_vcf_path,
            candidates=[candidate],
            reference_fasta_path=reference_fasta,
            output_path=output_path,
        )

    assert not output_path.exists()

    assert not Path(
        f"{output_path}.tbi"
    ).exists()

    assert not (
            tmp_path
            / ".diagnostic_candidates.reference-compatible.vcf"
    ).exists()

    assert not (
            tmp_path
            / ".diagnostic_candidates.normalized.unsorted.vcf.gz"
    ).exists()

    assert not (
            tmp_path
            / (
                ".diagnostic_candidates.normalized.vcf.gz"
                ".tmp.vcf.gz"
            )
    ).exists()

    assert candidate.is_normalized() is False