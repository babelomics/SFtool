# tests/variant_confirmation/test_candidate_vcf.py

from pathlib import Path

import pytest

from sftool.variant_confirmation.candidate_vcf import (
    CandidateVcfWriteError,
    CandidateVcfWriter,
    InvalidCandidateCollectionError,
)
from sftool.variant_confirmation.models import (
    VariantCandidate,
)


def build_candidate(
        candidate_id: str = "candidate_1",
        chromosome: str = "2",
        position: int = 47476391,
        reference: str = "C",
        alternate: str = "A",
        assembly: str = "GRCh38",
) -> VariantCandidate:
    """
    Build a candidate for candidate VCF unit tests.
    """
    return VariantCandidate(
        candidate_id=candidate_id,
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
        assembly=assembly,
    )


def read_lines(
        path: Path,
) -> list[str]:
    return path.read_text(
        encoding="utf-8"
    ).splitlines()


def get_records(
        path: Path,
) -> list[str]:
    return [
        line
        for line in read_lines(path)
        if not line.startswith("#")
    ]


def test_write_single_candidate_vcf(
        tmp_path,
):
    candidate = build_candidate()

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    writer = CandidateVcfWriter()

    result = writer.write(
        candidates=[candidate],
        output_path=output_path,
    )

    assert result == output_path
    assert output_path.exists()

    assert read_lines(output_path) == [
        "##fileformat=VCFv4.2",
        "##source=SFtool",
        "##reference=GRCh38",
        (
            '##INFO=<ID=ASSEMBLY,Number=1,Type=String,'
            'Description="Reference genome assembly">'
        ),
        (
            '##INFO=<ID=ORIGINAL_VARIANT,Number=1,Type=String,'
            'Description="Original genomic candidate before normalization">'
        ),
        "##contig=<ID=2>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        (
            "2\t47476391\tcandidate_1\tC\tA\t.\t.\t"
            "ASSEMBLY=GRCh38;ORIGINAL_VARIANT=2:47476391:C:A"
        ),
    ]


def test_write_multiple_candidates(
        tmp_path,
):
    candidates = [
        build_candidate(
            candidate_id="candidate_1",
            alternate="A",
        ),
        build_candidate(
            candidate_id="candidate_2",
            alternate="G",
        ),
    ]

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    writer = CandidateVcfWriter()

    writer.write(
        candidates=candidates,
        output_path=output_path,
    )

    records = get_records(
        output_path
    )

    assert records == [
        (
            "2\t47476391\tcandidate_1\tC\tA\t.\t.\t"
            "ASSEMBLY=GRCh38;ORIGINAL_VARIANT=2:47476391:C:A"
        ),
        (
            "2\t47476391\tcandidate_2\tC\tG\t.\t.\t"
            "ASSEMBLY=GRCh38;ORIGINAL_VARIANT=2:47476391:C:G"
        ),
    ]


def test_preserves_candidate_order(
        tmp_path,
):
    candidates = [
        build_candidate(
            candidate_id="candidate_2",
            chromosome="2",
            position=200,
        ),
        build_candidate(
            candidate_id="candidate_1",
            chromosome="1",
            position=100,
        ),
    ]

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    CandidateVcfWriter().write(
        candidates=candidates,
        output_path=output_path,
    )

    records = get_records(
        output_path
    )

    assert records[0].startswith(
        "2\t200\tcandidate_2"
    )

    assert records[1].startswith(
        "1\t100\tcandidate_1"
    )


@pytest.mark.parametrize(
    "assembly",
    [
        "GRCh37",
        "GRCh38",
    ],
)
def test_supports_both_reference_assemblies(
        tmp_path,
        assembly,
):
    candidate = build_candidate(
        assembly=assembly
    )

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    CandidateVcfWriter().write(
        candidates=[candidate],
        output_path=output_path,
    )

    lines = read_lines(
        output_path
    )

    assert f"##reference={assembly}" in lines
    assert (
            f"ASSEMBLY={assembly};"
            in get_records(output_path)[0]
    )


def test_write_to_directory_uses_standard_filename(
        tmp_path,
):
    output_directory = (
            tmp_path
            / "sample_1"
            / "variant_confirmation"
    )

    result = CandidateVcfWriter().write_to_directory(
        candidates=[build_candidate()],
        output_directory=output_directory,
    )

    assert result == (
            output_directory
            / "diagnostic_candidates.raw.vcf"
    )

    assert result.exists()


def test_creates_missing_output_directories(
        tmp_path,
):
    output_path = (
            tmp_path
            / "run"
            / "sample"
            / "variant_confirmation"
            / "diagnostic_candidates.raw.vcf"
    )

    CandidateVcfWriter().write(
        candidates=[build_candidate()],
        output_path=output_path,
    )

    assert output_path.exists()


def test_declares_each_contig_once(
        tmp_path,
):
    candidates = [
        build_candidate(
            candidate_id="candidate_1",
            chromosome="2",
            position=100,
        ),
        build_candidate(
            candidate_id="candidate_2",
            chromosome="2",
            position=200,
        ),
        build_candidate(
            candidate_id="candidate_3",
            chromosome="X",
            position=300,
        ),
    ]

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    CandidateVcfWriter().write(
        candidates=candidates,
        output_path=output_path,
    )

    contig_lines = [
        line
        for line in read_lines(output_path)
        if line.startswith("##contig=")
    ]

    assert contig_lines == [
        "##contig=<ID=2>",
        "##contig=<ID=X>",
    ]


def test_does_not_modify_candidate_coordinates(
        tmp_path,
):
    candidate = build_candidate(
        chromosome="chr2",
        position=47476391,
        reference="AC",
        alternate="A",
    )

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    CandidateVcfWriter().write(
        candidates=[candidate],
        output_path=output_path,
    )

    assert candidate.chromosome == "chr2"
    assert candidate.position == 47476391
    assert candidate.reference == "AC"
    assert candidate.alternate == "A"

    assert candidate.is_normalized() is False

    assert get_records(output_path) == [
        (
            "chr2\t47476391\tcandidate_1\tAC\tA\t.\t.\t"
            "ASSEMBLY=GRCh38;"
            "ORIGINAL_VARIANT=chr2:47476391:AC:A"
        )
    ]


def test_rejects_empty_candidate_collection(
        tmp_path,
):
    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    with pytest.raises(
            InvalidCandidateCollectionError,
            match="At least one",
    ):
        CandidateVcfWriter().write(
            candidates=[],
            output_path=output_path,
        )

    assert not output_path.exists()


@pytest.mark.parametrize(
    "invalid_candidates",
    [
        None,
        "candidate",
        {"candidate_1": "invalid"},
        iter([]),
    ],
)
def test_rejects_non_sequence_candidate_collection(
        tmp_path,
        invalid_candidates,
):
    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    with pytest.raises(
            InvalidCandidateCollectionError,
            match="must be a sequence",
    ):
        CandidateVcfWriter().write(
            candidates=invalid_candidates,
            output_path=output_path,
        )


def test_rejects_non_candidate_object(
        tmp_path,
):
    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    with pytest.raises(
            InvalidCandidateCollectionError,
            match="must be a VariantCandidate",
    ):
        CandidateVcfWriter().write(
            candidates=[
                build_candidate(),
                "invalid",
            ],
            output_path=output_path,
        )


def test_rejects_candidates_from_different_assemblies(
        tmp_path,
):
    candidates = [
        build_candidate(
            candidate_id="candidate_1",
            assembly="GRCh37",
        ),
        build_candidate(
            candidate_id="candidate_2",
            assembly="GRCh38",
        ),
    ]

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    with pytest.raises(
            InvalidCandidateCollectionError,
            match="same reference assembly",
    ):
        CandidateVcfWriter().write(
            candidates=candidates,
            output_path=output_path,
        )


def test_rejects_duplicate_candidate_ids(
        tmp_path,
):
    candidates = [
        build_candidate(
            candidate_id="candidate_1",
            alternate="A",
        ),
        build_candidate(
            candidate_id="candidate_1",
            alternate="G",
        ),
    ]

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    with pytest.raises(
            InvalidCandidateCollectionError,
            match="Candidate IDs must be unique",
    ):
        CandidateVcfWriter().write(
            candidates=candidates,
            output_path=output_path,
        )


def test_rejects_duplicate_genomic_candidates(
        tmp_path,
):
    candidates = [
        build_candidate(
            candidate_id="candidate_1",
        ),
        build_candidate(
            candidate_id="candidate_2",
        ),
    ]

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    with pytest.raises(
            InvalidCandidateCollectionError,
            match="Duplicate genomic candidates",
    ):
        CandidateVcfWriter().write(
            candidates=candidates,
            output_path=output_path,
        )


@pytest.mark.parametrize(
    (
            "reference",
            "alternate",
    ),
    [
        ("A", "<DEL>"),
        ("N", "A[2:100["),
        ("A*", "G"),
    ],
)
def test_rejects_non_sequence_alleles(
        tmp_path,
        reference,
        alternate,
):
    candidate = build_candidate(
        reference=reference,
        alternate=alternate,
    )

    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    with pytest.raises(
            InvalidCandidateCollectionError,
            match="invalid .* allele",
    ):
        CandidateVcfWriter().write(
            candidates=[candidate],
            output_path=output_path,
        )


@pytest.mark.parametrize(
    "invalid_output_path",
    [
        None,
        123,
        [],
    ],
)
def test_rejects_invalid_output_path_type(
        invalid_output_path,
):
    with pytest.raises(
            TypeError,
            match="output_path",
    ):
        CandidateVcfWriter().write(
            candidates=[build_candidate()],
            output_path=invalid_output_path,
        )


def test_requires_vcf_extension(
        tmp_path,
):
    output_path = (
            tmp_path
            / "diagnostic_candidates.txt"
    )

    with pytest.raises(
            ValueError,
            match=r"\.vcf extension",
    ):
        CandidateVcfWriter().write(
            candidates=[build_candidate()],
            output_path=output_path,
        )


def test_rejects_filename_with_directory_components(
        tmp_path,
):
    with pytest.raises(
            ValueError,
            match="must not contain directory components",
    ):
        CandidateVcfWriter().write_to_directory(
            candidates=[build_candidate()],
            output_directory=tmp_path,
            filename="nested/candidates.vcf",
        )


def test_existing_file_is_replaced(
        tmp_path,
):
    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    output_path.write_text(
        "old content\n",
        encoding="utf-8",
    )

    CandidateVcfWriter().write(
        candidates=[build_candidate()],
        output_path=output_path,
    )

    assert "old content" not in output_path.read_text(
        encoding="utf-8"
    )

    assert output_path.read_text(
        encoding="utf-8"
    ).startswith(
        "##fileformat=VCFv4.2\n"
    )


def test_write_error_is_wrapped(
        tmp_path,
        monkeypatch,
):
    output_path = (
            tmp_path
            / "diagnostic_candidates.raw.vcf"
    )

    def raise_os_error(
            source,
            destination,
    ):
        raise OSError(
            "simulated write failure"
        )

    monkeypatch.setattr(
        "sftool.variant_confirmation.candidate_vcf.os.replace",
        raise_os_error,
    )

    with pytest.raises(
            CandidateVcfWriteError,
            match="simulated write failure",
    ):
        CandidateVcfWriter().write(
            candidates=[build_candidate()],
            output_path=output_path,
        )

    assert not output_path.exists()