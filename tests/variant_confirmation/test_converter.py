# tests/variant_confirmation/test_converter.py

from unittest.mock import Mock

import pytest
import requests

from sftool.variant_confirmation.converter import (
    GeneBeAssemblyMismatchError,
    GeneBeNoCandidatesError,
    GeneBeRequestError,
    GeneBeResponseError,
    GeneBeVariantConverter,
)
from sftool.variant_confirmation.models import (
    VariantConfirmationRequest,
)


def build_request(
        variant: str,
        representation_type: str,
) -> VariantConfirmationRequest:
    """
    Build a parsed VariantConfirmationRequest for converter tests.
    """
    request = VariantConfirmationRequest(
        variant=variant
    )

    request.set_representation_type(
        representation_type
    )

    return request


def build_mock_converter(
        payload,
        timeout: float = 30.0,
):
    """
    Build a converter whose HTTP session returns a controlled response.
    """
    session = Mock()
    response = Mock()

    response.raise_for_status.return_value = None
    response.json.return_value = payload
    response.status_code = 200
    response.text = ""

    session.post.return_value = response

    converter = GeneBeVariantConverter(
        session=session,
        timeout=timeout,
    )

    return converter, session, response


def test_convert_hgvsc_to_single_candidate():
    request = build_request(
        variant="NM_000251.3:c.2030C>A",
        representation_type="hgvsc",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 47476391,
                    "ref": "C",
                    "alt": "A",
                }
            ]
        }
    ]

    converter, session, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    assert len(candidates) == 1

    candidate = candidates[0]

    assert candidate.candidate_id == "candidate_1"
    assert candidate.chromosome == "2"
    assert candidate.position == 47476391
    assert candidate.reference == "C"
    assert candidate.alternate == "A"
    assert candidate.assembly == "GRCh38"
    assert candidate.get_genomic_variant() == "2:47476391:C:A"
    assert candidate.conversion_warnings == []

    session.post.assert_called_once_with(
        converter.API_URL,
        params={"genome": "hg38"},
        json=["NM_000251.3:c.2030C>A"],
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
        },
        timeout=30.0,
    )


def test_convert_hgvsg_to_single_candidate():
    request = build_request(
        variant="NC_000002.12:g.47476391C>T",
        representation_type="hgvsg",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 47476391,
                    "ref": "C",
                    "alt": "T",
                }
            ]
        }
    ]

    converter, _, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    assert len(candidates) == 1
    assert candidates[0].get_genomic_variant() == (
        "2:47476391:C:T"
    )


def test_convert_genomic_representation():
    request = build_request(
        variant="2:47476391:C:T",
        representation_type="genomic",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 47476391,
                    "ref": "C",
                    "alt": "T",
                }
            ]
        }
    ]

    converter, session, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    assert len(candidates) == 1
    assert candidates[0].get_genomic_variant() == (
        "2:47476391:C:T"
    )

    assert session.post.call_args.kwargs["json"] == [
        "2:47476391:C:T"
    ]


def test_convert_hgvsp_to_multiple_candidates():
    request = build_request(
        variant="NP_000242.1:p.Pro677His",
        representation_type="hgvsp",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 47476391,
                    "ref": "C",
                    "alt": "A",
                },
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 47476391,
                    "ref": "C",
                    "alt": "G",
                },
            ]
        }
    ]

    converter, _, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    assert len(candidates) == 2

    assert candidates[0].candidate_id == "candidate_1"
    assert candidates[1].candidate_id == "candidate_2"

    assert candidates[0].get_genomic_variant() == (
        "2:47476391:C:A"
    )
    assert candidates[1].get_genomic_variant() == (
        "2:47476391:C:G"
    )

    expected_warning = (
        "Protein representation resolved to multiple "
        "genomic candidates"
    )

    assert candidates[0].conversion_warnings == [
        expected_warning
    ]
    assert candidates[1].conversion_warnings == [
        expected_warning
    ]


def test_multiple_non_protein_candidates_generate_warning():
    request = build_request(
        variant="NM_000251.3:c.2030C>A",
        representation_type="hgvsc",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 47476391,
                    "ref": "C",
                    "alt": "A",
                },
                {
                    "genome": "hg38",
                    "chr": "3",
                    "pos": 100,
                    "ref": "G",
                    "alt": "T",
                },
            ]
        }
    ]

    converter, _, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    expected_warning = (
        "Variant representation resolved to multiple "
        "genomic candidates"
    )

    assert len(candidates) == 2
    assert candidates[0].conversion_warnings == [
        expected_warning
    ]
    assert candidates[1].conversion_warnings == [
        expected_warning
    ]


@pytest.mark.parametrize(
    (
            "assembly",
            "genebe_genome",
            "returned_genome",
    ),
    [
        ("GRCh37", "hg19", "hg19"),
        ("GRCh38", "hg38", "hg38"),
    ],
)
def test_assembly_mapping(
        assembly,
        genebe_genome,
        returned_genome,
):
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": returned_genome,
                    "chr": "1",
                    "pos": 100,
                    "ref": "A",
                    "alt": "G",
                }
            ]
        }
    ]

    converter, session, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly=assembly,
    )

    assert candidates[0].assembly == assembly

    assert session.post.call_args.kwargs["params"] == {
        "genome": genebe_genome
    }


def test_duplicate_candidates_are_removed():
    request = build_request(
        variant="NP_000242.1:p.Pro677His",
        representation_type="hgvsp",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 100,
                    "ref": "C",
                    "alt": "A",
                },
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 100,
                    "ref": "C",
                    "alt": "A",
                },
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 100,
                    "ref": "C",
                    "alt": "G",
                },
            ]
        }
    ]

    converter, _, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    assert len(candidates) == 2

    assert candidates[0].candidate_id == "candidate_1"
    assert candidates[1].candidate_id == "candidate_2"

    assert candidates[0].alternate == "A"
    assert candidates[1].alternate == "G"

    assert candidates[0].conversion_warnings == [
        "Protein representation resolved to multiple "
        "genomic candidates",
        "GeneBe returned duplicate genomic candidates; "
        "duplicates were removed",
    ]


def test_duplicate_single_candidate_generates_only_duplicate_warning():
    request = build_request(
        variant="NM_000251.3:c.2030C>A",
        representation_type="hgvsc",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 100,
                    "ref": "C",
                    "alt": "A",
                },
                {
                    "genome": "hg38",
                    "chr": "2",
                    "pos": 100,
                    "ref": "C",
                    "alt": "A",
                },
            ]
        }
    ]

    converter, _, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    assert len(candidates) == 1

    assert candidates[0].conversion_warnings == [
        "GeneBe returned duplicate genomic candidates; "
        "duplicates were removed"
    ]


def test_alleles_are_converted_to_uppercase():
    request = build_request(
        variant="1:100:a:g",
        representation_type="genomic",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "1",
                    "pos": 100,
                    "ref": "a",
                    "alt": "g",
                }
            ]
        }
    ]

    converter, _, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    assert candidates[0].reference == "A"
    assert candidates[0].alternate == "G"


def test_chromosome_is_preserved_as_returned():
    request = build_request(
        variant="chr1:100:A:G",
        representation_type="genomic",
    )

    payload = [
        {
            "variants": [
                {
                    "genome": "hg38",
                    "chr": "chr1",
                    "pos": 100,
                    "ref": "A",
                    "alt": "G",
                }
            ]
        }
    ]

    converter, _, _ = build_mock_converter(
        payload
    )

    candidates = converter.convert(
        request=request,
        assembly="GRCh38",
    )

    assert candidates[0].chromosome == "chr1"


def test_rejects_request_with_wrong_type():
    converter = GeneBeVariantConverter(
        session=Mock()
    )

    with pytest.raises(
            TypeError,
            match="request must be a VariantConfirmationRequest",
    ):
        converter.convert(
            request="NM_000251.3:c.2030C>A",
            assembly="GRCh38",
        )


def test_rejects_unparsed_request():
    request = VariantConfirmationRequest(
        variant="NM_000251.3:c.2030C>A"
    )

    converter = GeneBeVariantConverter(
        session=Mock()
    )

    with pytest.raises(
            ValueError,
            match="has not been parsed",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_rejects_unsupported_assembly():
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter = GeneBeVariantConverter(
        session=Mock()
    )

    with pytest.raises(
            ValueError,
            match="Unsupported reference genome",
    ):
        converter.convert(
            request=request,
            assembly="GRCh36",
        )


@pytest.mark.parametrize(
    "invalid_timeout",
    [
        0,
        -1,
        -0.5,
    ],
)
def test_rejects_non_positive_timeout(
        invalid_timeout,
):
    with pytest.raises(
            ValueError,
            match="greater than zero",
    ):
        GeneBeVariantConverter(
            timeout=invalid_timeout
        )


@pytest.mark.parametrize(
    "invalid_timeout",
    [
        None,
        "30",
        True,
    ],
)
def test_rejects_invalid_timeout_type(
        invalid_timeout,
):
    with pytest.raises(
            TypeError,
            match="positive number",
    ):
        GeneBeVariantConverter(
            timeout=invalid_timeout
        )


def test_timeout_is_converted_to_request_error():
    request = build_request(
        variant="NM_000251.3:c.2030C>A",
        representation_type="hgvsc",
    )

    session = Mock()
    session.post.side_effect = requests.Timeout(
        "request timeout"
    )

    converter = GeneBeVariantConverter(
        session=session
    )

    with pytest.raises(
            GeneBeRequestError,
            match="timed out",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_connection_error_is_converted_to_request_error():
    request = build_request(
        variant="NM_000251.3:c.2030C>A",
        representation_type="hgvsc",
    )

    session = Mock()
    session.post.side_effect = requests.ConnectionError(
        "connection refused"
    )

    converter = GeneBeVariantConverter(
        session=session
    )

    with pytest.raises(
            GeneBeRequestError,
            match="Could not connect",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


@pytest.mark.parametrize(
    "status_code",
    [
        400,
        429,
        500,
    ],
)
def test_http_error_is_converted_to_request_error(
        status_code,
):
    request = build_request(
        variant="NM_000251.3:c.2030C>A",
        representation_type="hgvsc",
    )

    session = Mock()
    response = Mock()

    response.status_code = status_code
    response.text = "GeneBe service error"

    http_error = requests.HTTPError(
        f"HTTP {status_code}",
        response=response,
    )

    response.raise_for_status.side_effect = http_error
    session.post.return_value = response

    converter = GeneBeVariantConverter(
        session=session
    )

    with pytest.raises(
            GeneBeRequestError,
            match=str(status_code),
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_invalid_json_is_rejected():
    request = build_request(
        variant="NM_000251.3:c.2030C>A",
        representation_type="hgvsc",
    )

    session = Mock()
    response = Mock()

    response.raise_for_status.return_value = None
    response.json.side_effect = ValueError(
        "invalid JSON"
    )

    session.post.return_value = response

    converter = GeneBeVariantConverter(
        session=session
    )

    with pytest.raises(
            GeneBeResponseError,
            match="non-JSON",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


@pytest.mark.parametrize(
    "payload",
    [
        {},
        "invalid",
        None,
    ],
)
def test_top_level_response_must_be_list(
        payload,
):
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        payload
    )

    with pytest.raises(
            GeneBeResponseError,
            match="must be a list",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


@pytest.mark.parametrize(
    "payload",
    [
        [],
        [
            {"variants": []},
            {"variants": []},
        ],
    ],
)
def test_response_must_contain_exactly_one_result(
        payload,
):
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        payload
    )

    with pytest.raises(
            GeneBeResponseError,
            match="exactly one conversion result",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_conversion_result_must_be_dictionary():
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        ["invalid"]
    )

    with pytest.raises(
            GeneBeResponseError,
            match="must be a dictionary",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_missing_variants_field_is_rejected():
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        [{"other": []}]
    )

    with pytest.raises(
            GeneBeResponseError,
            match="does not contain 'variants'",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_genebe_error_message_without_variants_is_conversion_failure():
    request = build_request(
        variant="NM_INVALID:c.1A>G",
        representation_type="hgvsc",
    )

    converter, _, _ = build_mock_converter(
        [
            {
                "error": "Unable to parse variant"
            }
        ]
    )

    with pytest.raises(
            GeneBeNoCandidatesError,
            match="Unable to parse variant",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_variants_field_must_be_list():
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": {}
            }
        ]
    )

    with pytest.raises(
            GeneBeResponseError,
            match="'variants' field must be a list",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_empty_variants_list_is_conversion_failure():
    request = build_request(
        variant="NM_INVALID:c.1A>G",
        representation_type="hgvsc",
    )

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": []
            }
        ]
    )

    with pytest.raises(
            GeneBeNoCandidatesError,
            match="no genomic candidates",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_candidate_must_be_dictionary():
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": [
                    "invalid"
                ]
            }
        ]
    )

    with pytest.raises(
            GeneBeResponseError,
            match="candidate 1 must be a dictionary",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


@pytest.mark.parametrize(
    "missing_field",
    [
        "genome",
        "chr",
        "pos",
        "ref",
        "alt",
    ],
)
def test_missing_candidate_field_is_rejected(
        missing_field,
):
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    candidate = {
        "genome": "hg38",
        "chr": "1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
    }

    del candidate[missing_field]

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": [
                    candidate
                ]
            }
        ]
    )

    with pytest.raises(
            GeneBeResponseError,
            match="missing required fields",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


@pytest.mark.parametrize(
    "position",
    [
        0,
        -1,
        1.5,
        "100",
        True,
        None,
    ],
)
def test_invalid_candidate_position_is_rejected(
        position,
):
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": [
                    {
                        "genome": "hg38",
                        "chr": "1",
                        "pos": position,
                        "ref": "A",
                        "alt": "G",
                    }
                ]
            }
        ]
    )

    with pytest.raises(
            GeneBeResponseError,
            match="positive integer",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


@pytest.mark.parametrize(
    (
            "field_name",
            "invalid_value",
    ),
    [
        ("genome", ""),
        ("genome", None),
        ("chr", ""),
        ("chr", None),
        ("ref", ""),
        ("ref", None),
        ("alt", ""),
        ("alt", None),
    ],
)
def test_invalid_string_candidate_field_is_rejected(
        field_name,
        invalid_value,
):
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    candidate = {
        "genome": "hg38",
        "chr": "1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
    }

    candidate[field_name] = invalid_value

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": [
                    candidate
                ]
            }
        ]
    )

    with pytest.raises(
            GeneBeResponseError,
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_unknown_genebe_genome_is_rejected():
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": [
                    {
                        "genome": "t2t",
                        "chr": "1",
                        "pos": 100,
                        "ref": "A",
                        "alt": "G",
                    }
                ]
            }
        ]
    )

    with pytest.raises(
            GeneBeResponseError,
            match="unsupported genome",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


def test_assembly_mismatch_is_rejected():
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": [
                    {
                        "genome": "hg19",
                        "chr": "1",
                        "pos": 100,
                        "ref": "A",
                        "alt": "G",
                    }
                ]
            }
        ]
    )

    with pytest.raises(
            GeneBeAssemblyMismatchError,
            match="hg19",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )


@pytest.mark.parametrize(
    "field_name",
    [
        "ref",
        "alt",
    ],
)
def test_dot_allele_is_rejected(
        field_name,
):
    request = build_request(
        variant="1:100:A:G",
        representation_type="genomic",
    )

    candidate = {
        "genome": "hg38",
        "chr": "1",
        "pos": 100,
        "ref": "A",
        "alt": "G",
    }

    candidate[field_name] = "."

    converter, _, _ = build_mock_converter(
        [
            {
                "variants": [
                    candidate
                ]
            }
        ]
    )

    with pytest.raises(
            GeneBeResponseError,
            match="invalid allele",
    ):
        converter.convert(
            request=request,
            assembly="GRCh38",
        )