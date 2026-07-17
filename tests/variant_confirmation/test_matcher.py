# tests/variant_confirmation/test_matcher.py

import gzip
import shutil
import subprocess

from pathlib import Path

import pytest
import vcfpy

from sftool.variant_confirmation.matcher import (
    CandidateMatcher,
    InvalidMatcherInputError,
    PatientVcfError,
)
from sftool.variant_confirmation.models import (
    VariantCandidate,
)


def build_normalized_candidate(
        candidate_id: str = "candidate_1",
        chromosome: str = "2",
        position: int = 100,
        reference: str = "A",
        alternate: str = "G",
        assembly: str = "GRCh38",
) -> VariantCandidate:
    """
    Build a normalized diagnostic candidate for matcher unit tests.
    """
    candidate = VariantCandidate(
        candidate_id=candidate_id,
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
        assembly=assembly,
    )

    candidate.set_normalized_coordinates(
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
    )

    return candidate


def write_patient_vcf(
        output_path: Path,
        records: list[str],
        sample_names: list[str] | None = None,
) -> Path:
    """
    Write a minimal compressed patient VCF for matcher unit tests.
    """
    if sample_names is None:
        sample_names = [
            "patient"
        ]

    header = (
            "##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,'
            'Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,'
            'Description="Read depth">\n'
            '##FILTER=<ID=LowQual,Description="Low-quality variant">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(sample_names)
            + "\n"
    )

    with gzip.open(
            output_path,
            "wt",
            encoding="utf-8",
    ) as output:
        output.write(
            header
        )

        for record in records:
            output.write(
                record
            )

            output.write(
                "\n"
            )

    Path(
        f"{output_path}.tbi"
    ).write_bytes(
        b"fake patient VCF index"
    )

    return output_path


def read_vcf_records(
        vcf_path: Path,
) -> list[vcfpy.Record]:
    """
    Read all records from a generated VCF.
    """
    reader = vcfpy.Reader.from_path(
        str(vcf_path)
    )

    try:
        return list(
            reader
        )

    finally:
        reader.close()


class FakeBcftoolsRunner:
    """
    Simulate the bcftools commands used by CandidateMatcher.

    Regional extraction copies the patient VCF records to an uncompressed VCF.
    Compression creates a gzip VCF and indexing creates a placeholder TBI file.
    """

    def __init__(self):
        self.commands: list[list[str]] = []

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

        if subcommand == "view":
            output_path = Path(
                command[
                    command.index("--output") + 1
                    ]
            )

            input_path = Path(
                command[-1]
            )

            output_type = command[
                command.index("--output-type") + 1
                ]

            if output_type == "v":
                self._extract_regions(
                    input_path=input_path,
                    output_path=output_path,
                )

            elif output_type == "z":
                self._compress_vcf(
                    input_path=input_path,
                    output_path=output_path,
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

    @staticmethod
    def _extract_regions(
            input_path: Path,
            output_path: Path,
    ) -> None:
        """
        Copy the compressed patient VCF to an uncompressed regional VCF.
        """
        with gzip.open(
                input_path,
                "rt",
                encoding="utf-8",
        ) as source:
            with output_path.open(
                    "w",
                    encoding="utf-8",
            ) as destination:
                shutil.copyfileobj(
                    source,
                    destination,
                )

    @staticmethod
    def _compress_vcf(
            input_path: Path,
            output_path: Path,
    ) -> None:
        """
        Compress an uncompressed VCF using gzip for unit testing.
        """
        with input_path.open(
                "rt",
                encoding="utf-8",
        ) as source:
            with gzip.open(
                    output_path,
                    "wt",
                    encoding="utf-8",
            ) as destination:
                shutil.copyfileobj(
                    source,
                    destination,
                )


@pytest.fixture
def fake_bcftools(
        tmp_path,
) -> Path:
    """
    Create an executable placeholder accepted by CandidateMatcher.
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


# ---------------------------------------------------------------------------
# Finds an exact candidate and retains its informational patient VCF fields.
# ---------------------------------------------------------------------------
def test_match_exact_candidate(
        tmp_path,
        fake_bcftools,
):
    """
    Verify that an exact CHROM, POS, REF and ALT match is detected.

    The test also checks that genotype, quality, filters and sample FORMAT
    fields are retained for later reporting.
    """
    candidate = build_normalized_candidate()

    patient_vcf_path = write_patient_vcf(
        output_path=tmp_path / "patient.norm.vcf.gz",
        records=[
            "2\t100\t.\tA\tG\t99\tPASS\t.\tGT:DP\t0/1:42"
        ],
    )

    runner = FakeBcftoolsRunner()

    matcher = CandidateMatcher(
        bcftools_path=fake_bcftools,
        runner=runner,
    )

    output_path = (
            tmp_path
            / "diagnostic_candidates.matches.vcf.gz"
    )

    result = matcher.match(
        candidates=[
            candidate
        ],
        patient_vcf_path=patient_vcf_path,
        output_path=output_path,
    )

    assert result.vcf_path == output_path
    assert output_path.is_file()
    assert Path(f"{output_path}.tbi").is_file()

    assert len(result.matches) == 1

    variant_match = result.matches[0]

    assert variant_match.candidate_id == "candidate_1"
    assert variant_match.found is True
    assert variant_match.genotype == "0/1"
    assert variant_match.quality == 99.0
    assert variant_match.filters == [
        "PASS"
    ]
    assert variant_match.sample_format == {
        "GT": "0/1",
        "DP": 42,
    }
    assert variant_match.warnings == []

    output_records = read_vcf_records(
        output_path
    )

    assert len(output_records) == 1
    assert output_records[0].CHROM == "2"
    assert output_records[0].POS == 100
    assert output_records[0].REF == "A"
    assert output_records[0].ALT[0].value == "G"


# ---------------------------------------------------------------------------
# Does not match a patient record when only its genomic position is identical.
# ---------------------------------------------------------------------------
def test_same_position_with_different_allele_is_not_a_match(
        tmp_path,
        fake_bcftools,
):
    """
    Verify that sharing the same genomic position is not sufficient.

    A patient record with a different ALT allele must not be associated with
    the diagnostic candidate.
    """
    candidate = build_normalized_candidate(
        alternate="G"
    )

    patient_vcf_path = write_patient_vcf(
        output_path=tmp_path / "patient.norm.vcf.gz",
        records=[
            "2\t100\t.\tA\tT\t80\tPASS\t.\tGT:DP\t0/1:30"
        ],
    )

    matcher = CandidateMatcher(
        bcftools_path=fake_bcftools,
        runner=FakeBcftoolsRunner(),
    )

    result = matcher.match(
        candidates=[
            candidate
        ],
        patient_vcf_path=patient_vcf_path,
        output_path=(
                tmp_path
                / "diagnostic_candidates.matches.vcf.gz"
        ),
    )

    variant_match = result.matches[0]

    assert variant_match.found is False
    assert variant_match.genotype is None
    assert variant_match.sample_format == {}
    assert variant_match.warnings == [
        "Diagnostic variant candidate was not found in the "
        "normalized patient VCF"
    ]

    output_records = read_vcf_records(
        result.vcf_path
    )

    assert output_records == []


# ---------------------------------------------------------------------------
# Returns one VariantMatch per candidate while preserving candidate order.
# ---------------------------------------------------------------------------
def test_match_multiple_candidates(
        tmp_path,
        fake_bcftools,
):
    """
    Verify matching when more than one diagnostic candidate is provided.

    The result must contain one VariantMatch per candidate in input order,
    including candidates that are not present in the patient VCF.
    """
    first_candidate = build_normalized_candidate(
        candidate_id="candidate_1",
        position=100,
        alternate="G",
    )

    second_candidate = build_normalized_candidate(
        candidate_id="candidate_2",
        position=200,
        reference="C",
        alternate="T",
    )

    patient_vcf_path = write_patient_vcf(
        output_path=tmp_path / "patient.norm.vcf.gz",
        records=[
            "2\t200\t.\tC\tT\t75\tLowQual\t.\tGT:DP\t1/1:18"
        ],
    )

    matcher = CandidateMatcher(
        bcftools_path=fake_bcftools,
        runner=FakeBcftoolsRunner(),
    )

    result = matcher.match(
        candidates=[
            first_candidate,
            second_candidate,
        ],
        patient_vcf_path=patient_vcf_path,
        output_path=(
                tmp_path
                / "diagnostic_candidates.matches.vcf.gz"
        ),
    )

    assert [
               variant_match.candidate_id
               for variant_match in result.matches
           ] == [
               "candidate_1",
               "candidate_2",
           ]

    first_match = result.matches[0]
    second_match = result.matches[1]

    assert first_match.found is False
    assert second_match.found is True

    assert second_match.genotype == "1/1"
    assert second_match.filters == [
        "LowQual"
    ]
    assert second_match.sample_format["DP"] == 18
    assert second_match.warnings == [
        "Matched variant does not pass the patient VCF filters: LowQual"
    ]


# ---------------------------------------------------------------------------
# Rejects a patient VCF containing more than one sample.
# ---------------------------------------------------------------------------
def test_reject_multisample_patient_vcf(
        tmp_path,
        fake_bcftools,
):
    """
    Verify that CandidateMatcher only accepts single-sample patient VCFs.

    A VCF declaring multiple samples must raise PatientVcfError before
    candidate matching starts.
    """
    candidate = build_normalized_candidate()

    patient_vcf_path = write_patient_vcf(
        output_path=tmp_path / "patient.norm.vcf.gz",
        sample_names=[
            "patient_1",
            "patient_2",
        ],
        records=[
            (
                "2\t100\t.\tA\tG\t99\tPASS\t.\tGT:DP\t"
                "0/1:42\t1/1:35"
            )
        ],
    )

    matcher = CandidateMatcher(
        bcftools_path=fake_bcftools,
        runner=FakeBcftoolsRunner(),
    )

    with pytest.raises(
            PatientVcfError,
            match="must contain exactly one sample",
    ):
        matcher.match(
            candidates=[
                candidate
            ],
            patient_vcf_path=patient_vcf_path,
            output_path=(
                    tmp_path
                    / "diagnostic_candidates.matches.vcf.gz"
            ),
        )


# ---------------------------------------------------------------------------
# Rejects candidates that have not passed through candidate VCF normalization.
# ---------------------------------------------------------------------------
def test_reject_unnormalized_candidate(
        tmp_path,
        fake_bcftools,
):
    """
    Verify that only normalized diagnostic candidates can be matched.

    A candidate without normalized genomic coordinates must raise
    InvalidMatcherInputError.
    """
    candidate = VariantCandidate(
        candidate_id="candidate_1",
        chromosome="2",
        position=100,
        reference="A",
        alternate="G",
        assembly="GRCh38",
    )

    patient_vcf_path = write_patient_vcf(
        output_path=tmp_path / "patient.norm.vcf.gz",
        records=[],
    )

    matcher = CandidateMatcher(
        bcftools_path=fake_bcftools,
        runner=FakeBcftoolsRunner(),
    )

    with pytest.raises(
            InvalidMatcherInputError,
            match="has not been normalized",
    ):
        matcher.match(
            candidates=[
                candidate
            ],
            patient_vcf_path=patient_vcf_path,
            output_path=(
                    tmp_path
                    / "diagnostic_candidates.matches.vcf.gz"
            ),
        )