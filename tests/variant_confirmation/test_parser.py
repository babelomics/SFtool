# tests/variant_confirmation/test_parser.py

import pytest

from sftool.variant_confirmation.models import (
    VariantConfirmationRequest,
)
from sftool.variant_confirmation.parser import (
    UnsupportedVariantRepresentationError,
    VariantRepresentationParser,
)


@pytest.fixture
def parser():
    return VariantRepresentationParser()


@pytest.mark.parametrize(
    "variant",
    [
        "2:47476391:C:T",
        "chr2:47476391:C:T",
        "X:12345:A:G",
        "chrX:12345:A:G",
        "MT:123:A:G",
        "chrM:123:A:G",
    ],
)
def test_detect_genomic_variants(parser, variant):
    assert parser.detect(variant) == "genomic"

@pytest.mark.parametrize(
    "variant",
    [
        "chr2:47476391:CT:C",
        "chr2:47476391:C:CT",
        "2:47476391:ACGT:A",
        "2:47476391:A:ACGT",
        "8:342345233:TATC:T",
        "8:342345233:T:TGGA",
    ],
)
def test_detect_genomic_indels(parser, variant):
    assert parser.detect(variant) == "genomic"


@pytest.mark.parametrize(
    "variant",
    [
        "NM_000251.3:c.2030C>A",
        "NM_000251.3:c.2030_2032del",
        "NM_000251.3:c.942+3A>T",
        "NM_000251.3:c.-15C>T",
        "NM_000251.3:c.2030_2031insA",
        "NM_000251.3:c.2030_2032delinsTT",
    ],
)
def test_detect_hgvsc_variants(parser, variant):
    assert parser.detect(variant) == "hgvsc"


@pytest.mark.parametrize(
    "variant",
    [
        "NC_000002.12:g.47476391C>T",
        "NC_000002.12:g.47476391_47476393del",
        "NG_007110.2:g.5000A>G",
    ],
)
def test_detect_hgvsg_variants(parser, variant):
    assert parser.detect(variant) == "hgvsg"


@pytest.mark.parametrize(
    "variant",
    [
        "NP_000242.1:p.Pro677Gln",
        "NP_000242.1:p.(Pro677Gln)",
        "NP_000242.1:p.Pro677GlnfsTer12",
        "NP_000242.1:p.Pro677*",
        "NP_000242.1:p.?",
    ],
)
def test_detect_hgvsp_variants(parser, variant):
    assert parser.detect(variant) == "hgvsp"


def test_parse_sets_request_representation_type(parser):
    request = VariantConfirmationRequest(
        "NM_000251.3:c.2030C>A"
    )

    returned_request = parser.parse(request)

    assert returned_request is request
    assert request.representation_type == "hgvsc"


@pytest.mark.parametrize(
    ("variant", "expected_type"),
    [
        ("2:47476391:C:T", "genomic"),
        ("NM_000251.3:c.2030C>A", "hgvsc"),
        ("NC_000002.12:g.47476391C>T", "hgvsg"),
        ("NP_000242.1:p.Pro677Gln", "hgvsp"),
    ],
)
def test_parse_assigns_expected_representation_type(
        parser,
        variant,
        expected_type,
):
    request = VariantConfirmationRequest(variant)

    parser.parse(request)

    assert request.representation_type == expected_type


def test_parse_rejects_non_request_object(parser):
    with pytest.raises(
            TypeError,
            match="request must be a VariantConfirmationRequest",
    ):
        parser.parse("NM_000251.3:c.2030C>A")



@pytest.mark.parametrize(
    "variant",
    [
        "prefix-chr2:47476391:C:T",
        "chr2:47476391:C:T-extra",
    ],
)
def test_detect_rejects_partial_matches(parser, variant):
    with pytest.raises(
            UnsupportedVariantRepresentationError
    ):
        parser.detect(variant)


@pytest.mark.parametrize(
    "variant",
    [
        "NM_000251.3:c.",
        "NC_000002.12:g.",
        "NP_000242.1:p.",
    ],
)
def test_detect_rejects_missing_hgvs_description(
        parser,
        variant,
):
    with pytest.raises(
            UnsupportedVariantRepresentationError
    ):
        parser.detect(variant)


@pytest.mark.parametrize(
    "variant",
    [
        "chr2:0:C:T",
        "chr2:-1:C:T",
        "chr2:position:C:T",
    ],
)
def test_detect_rejects_invalid_genomic_position(
        parser,
        variant,
):
    with pytest.raises(
            UnsupportedVariantRepresentationError
    ):
        parser.detect(variant)


@pytest.mark.parametrize(
    "variant",
    [
        "chr0:47476391:C:T",
        "chr23:47476391:C:T",
        "chr25:47476391:C:T",
        "chrZ:47476391:C:T",
    ],
)
def test_detect_rejects_invalid_chromosome(
        parser,
        variant,
):
    with pytest.raises(
            UnsupportedVariantRepresentationError
    ):
        parser.detect(variant)


@pytest.mark.parametrize(
    "variant",
    [
        "chr2:47476391::T",
        "chr2:47476391:C:",
        "chr2:47476391:-:T",
        "chr2:47476391:C:-",
        "chr2:47476391:C:T:G",

    ],
)
def test_detect_rejects_invalid_genomic_alleles(
        parser,
        variant,
):
    with pytest.raises(
            UnsupportedVariantRepresentationError
    ):
        parser.detect(variant)


@pytest.mark.parametrize(
    "variant",
    [
        "8:342345233:ATC:",
        "8:342345233::GGA",
    ],
)
def test_detect_rejects_indels_without_anchor_base(
        parser,
        variant,
):
    with pytest.raises(
            UnsupportedVariantRepresentationError,
            match=(
                    "REF and ALT must both contain at least one nucleotide"
            ),
    ):
        parser.detect(variant)

@pytest.mark.parametrize(
    ("variant", "expected_type"),
    [
        ("CHR2:47476391:c:t", "genomic"),
        ("nm_000251.3:C.2030c>a", "hgvsc"),
        ("nc_000002.12:G.47476391c>t", "hgvsg"),
        ("np_000242.1:P.pro677gln", "hgvsp"),
    ],
)
def test_detection_is_case_insensitive(
        parser,
        variant,
        expected_type,
):
    assert parser.detect(variant) == expected_type


def test_detect_strips_outer_whitespace(parser):
    assert (
            parser.detect("  NM_000251.3:c.2030C>A  ")
            == "hgvsc"
    )


@pytest.mark.parametrize(
    "variant",
    [
        "NM_000251.3: c.2030C>A",
        "NM_000251.3:c.2030 C>A",
        "chr2: 47476391:C:T",
        "chr2:47476391:C: T",
    ],
)
def test_detect_rejects_internal_whitespace(
        parser,
        variant,
):
    with pytest.raises(
            UnsupportedVariantRepresentationError
    ):
        parser.detect(variant)


@pytest.mark.parametrize(
    "variant",
    [
        None,
        123,
        ["NM_000251.3:c.2030C>A"],
        {"variant": "NM_000251.3:c.2030C>A"},
    ],
)
def test_detect_rejects_non_string_values(parser, variant):
    with pytest.raises(
            TypeError,
            match="variant must be a string",
    ):
        parser.detect(variant)


@pytest.mark.parametrize(
    "variant",
    [
        "",
        " ",
        "\t",
        "\n",
    ],
)
def test_detect_rejects_empty_strings(parser, variant):
    with pytest.raises(
            ValueError,
            match="variant must be a non-empty string",
    ):
        parser.detect(variant)