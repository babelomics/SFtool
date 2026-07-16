import pytest

from sftool.core.bootstrap import (
    validate_execution_block,
    validate_variant_confirmation_requests
)
from sftool.utils.errors import ValidationError

### Tests for sftool.core.bootstrap.validate_execution_block

@pytest.mark.parametrize(
    "modes,num_samples",
    [
        (
                ["secondary_findings_discovery"],
                1,
        ),
        (
                ["variant_confirmation"],
                1,
        ),
        (
                [
                    "secondary_findings_discovery",
                    "variant_confirmation",
                ],
                1,
        ),
        (
                ["secondary_findings_discovery"],
                2,
        ),
        (
                [
                    "secondary_findings_discovery",
                    "variant_confirmation",
                ],
                2,
        ),
    ],
)

def test_validate_execution_block_accepts_supported_modes(
        modes,
        num_samples,
):
    execution = {
        "modes": modes,
        "reference_genome": "GRCh38",
    }

    validate_execution_block(
        execution,
        num_samples=num_samples,
    )


def test_validate_execution_block_rejects_empty_modes():
    with pytest.raises(
            ValidationError,
            match="at least one workflow",
    ):
        validate_execution_block(
            {"modes": []},
            num_samples=1,
        )


def test_two_samples_require_secondary_findings():
    with pytest.raises(
            ValidationError,
            match="Two-sample executions require",
    ):
        validate_execution_block(
            {
                "modes": [
                    "variant_confirmation"
                ]
            },
            num_samples=2,
        )


# Variant confirmation test per sample

def test_valid_variant_confirmation_request():
    samples = [
        {
            "sample_id": "proband",
            "variant_confirmation": {
                "variant": "NM_000251.3:c.2030C>A"
            },
        }
    ]

    validate_variant_confirmation_requests(
        samples,
        modes=["variant_confirmation"],
    )


def test_confirmation_requires_at_least_one_request():
    samples = [
        {
            "sample_id": "proband",
        }
    ]

    with pytest.raises(
            ValidationError,
            match="At least one sample",
    ):
        validate_variant_confirmation_requests(
            samples,
            modes=["variant_confirmation"],
        )

def test_confirmation_request_rejects_extra_fields():
    samples = [
        {
            "sample_id": "proband",
            "variant_confirmation": {
                "variant": "NM_000251.3:c.2030C>A",
                "enabled": True,
            },
        }
    ]

    with pytest.raises(
            ValidationError,
            match="exactly one field",
    ):
        validate_variant_confirmation_requests(
            samples,
            modes=["variant_confirmation"],
        )


@pytest.mark.parametrize(
    "variant",
    [
        "",
        " ",
        "   ",
    ],
)

def test_confirmation_request_rejects_empty_variant(
        variant,
):
    samples = [
        {
            "sample_id": "proband",
            "variant_confirmation": {
                "variant": variant,
            },
        }
    ]

    with pytest.raises(
            ValidationError,
            match="non-empty string",
    ):
        validate_variant_confirmation_requests(
            samples,
            modes=["variant_confirmation"],
        )