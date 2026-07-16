import pytest

from sftool.variant_confirmation.models import (
    VariantConfirmationRequest,
    VariantCandidate,
    VariantConfirmationResult,
    VariantMatch

)
from sftool.utils.errors import ValidationError



def test_request_from_dict():
    request = VariantConfirmationRequest.from_dict(
        {
            "variant": "NM_000251.3:c.2030C>A"
        }
    )

    assert request.variant == "NM_000251.3:c.2030C>A"
    assert request.representation_type is None


def test_candidate_requires_normalization_before_matching():
    candidate = VariantCandidate(
        candidate_id="candidate_001",
        chromosome="2",
        position=47476391,
        reference="C",
        alternate="T",
        assembly="GRCh38",
    )

    with pytest.raises(ValueError):
        candidate.get_matching_key()


def test_result_rejects_match_for_unknown_candidate():
    request = VariantConfirmationRequest(
        "NM_000251.3:c.2030C>A"
    )

    result = VariantConfirmationResult(request)

    match = VariantMatch(
        candidate_id="candidate_999",
        chromosome="2",
        position=47476391,
        reference="C",
        alternate="T",
    )

    with pytest.raises(ValueError):
        result.add_match(match)


def test_result_requires_one_match_per_candidate():
    request = VariantConfirmationRequest(
        "NM_000251.3:c.2030C>A"
    )

    result = VariantConfirmationResult(request)

    result.add_candidate(
        VariantCandidate(
            candidate_id="candidate_001",
            chromosome="2",
            position=47476391,
            reference="C",
            alternate="T",
            assembly="GRCh38",
        )
    )

    with pytest.raises(ValueError):
        result.to_dict()


def test_ambiguous_result():
    request = VariantConfirmationRequest(
        "NP_000242.1:p.Ser677Phe",
        representation_type="hgvsp",
    )

    result = VariantConfirmationResult(request)

    # Add two candidates and their respective matches.

    assert result.is_ambiguous() is True