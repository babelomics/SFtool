"""
Utilities for loading and validating SFtool resource specifications.
"""

from __future__ import annotations

import json
from importlib import resources
from typing import Any


BUNDLED_RESOURCE_PACKAGE = "sftool.data.resources"
BUNDLED_RESOURCE_FILENAME = "bundled_resources.json"
SUPPORTED_BUNDLED_SCHEMA_VERSIONS = {1}


class ResourceSpecificationError(ValueError):
    """
    Raised when a bundled resource specification is missing or invalid.
    """


def load_bundled_resources() -> dict[str, Any]:
    """
    Load and validate the resource specification bundled with SFtool.

    Returns
    -------
    dict[str, Any]
        Parsed bundled resource specification.

    Raises
    ------
    ResourceSpecificationError
        If the package resource is missing, malformed, or unsupported.
    """

    try:
        resource_file = resources.files(
            BUNDLED_RESOURCE_PACKAGE
        ).joinpath(BUNDLED_RESOURCE_FILENAME)
    except (ModuleNotFoundError, AttributeError) as exc:
        raise ResourceSpecificationError(
            "Unable to locate the bundled SFtool resource package."
        ) from exc

    try:
        raw_content = resource_file.read_text(encoding="utf-8")
    except (FileNotFoundError, OSError) as exc:
        raise ResourceSpecificationError(
            "Unable to read the bundled SFtool resource specification: "
            f"{BUNDLED_RESOURCE_FILENAME}"
        ) from exc

    try:
        specification = json.loads(raw_content)
    except json.JSONDecodeError as exc:
        raise ResourceSpecificationError(
            "The bundled SFtool resource specification contains "
            f"invalid JSON: {exc}"
        ) from exc

    validate_bundled_resources(specification)

    return specification


def validate_bundled_resources(
        specification: dict[str, Any],
) -> None:
    """
    Validate the basic bundled-resource specification structure.

    This performs structural validation only. Network availability and
    remote resource existence are checked during resource setup.
    """

    if not isinstance(specification, dict):
        raise ResourceSpecificationError(
            "The bundled resource specification must be a JSON object."
        )

    schema_version = specification.get("schema_version")

    if schema_version not in SUPPORTED_BUNDLED_SCHEMA_VERSIONS:
        raise ResourceSpecificationError(
            "Unsupported bundled resource schema version: "
            f"{schema_version!r}. Supported versions: "
            f"{sorted(SUPPORTED_BUNDLED_SCHEMA_VERSIONS)}"
        )

    required_sections = {
        "clinvar",
        "hpo",
        "pharmcat",
        "reference_genomes",
    }

    missing_sections = sorted(
        required_sections.difference(specification)
    )

    if missing_sections:
        raise ResourceSpecificationError(
            "The bundled resource specification is missing required "
            f"sections: {', '.join(missing_sections)}"
        )

    _validate_clinvar_specification(specification["clinvar"])
    _validate_download_specification(
        name="hpo",
        specification=specification["hpo"],
    )
    _validate_download_specification(
        name="pharmcat",
        specification=specification["pharmcat"],
    )
    _validate_reference_genomes(
        specification["reference_genomes"]
    )

def _validate_clinvar_specification(
        specification: Any,
) -> None:
    if not isinstance(specification, dict):
        raise ResourceSpecificationError(
            "The 'clinvar' resource definition must be a JSON object."
        )

    required_fields = {
        "version",
        "archive_base_url",
        "variant_summary_filename",
        "submission_summary_filename",
    }

    _require_non_empty_string_fields(
        section_name="clinvar",
        specification=specification,
        required_fields=required_fields,
    )

    version = specification["version"]

    for field_name in (
            "variant_summary_filename",
            "submission_summary_filename",
    ):
        template = specification[field_name]

        try:
            rendered_filename = template.format(version=version)
        except (KeyError, ValueError) as exc:
            raise ResourceSpecificationError(
                f"Invalid ClinVar filename template "
                f"'{field_name}': {template!r}"
            ) from exc

        if not rendered_filename:
            raise ResourceSpecificationError(
                f"ClinVar filename template '{field_name}' "
                "produced an empty filename."
            )


def _validate_download_specification(
        *,
        name: str,
        specification: Any,
) -> None:
    if not isinstance(specification, dict):
        raise ResourceSpecificationError(
            f"The '{name}' resource definition must be a JSON object."
        )

    _require_non_empty_string_fields(
        section_name=name,
        specification=specification,
        required_fields={
            "version",
            "url",
            "filename",
        },
    )


def _validate_reference_genomes(
        specification: Any,
) -> None:
    if not isinstance(specification, dict):
        raise ResourceSpecificationError(
            "The 'reference_genomes' definition must be a JSON object."
        )

    for assembly in ("GRCh37", "GRCh38"):
        if assembly not in specification:
            raise ResourceSpecificationError(
                f"The reference-genome specification is missing assembly: {assembly}"
            )

        genome_specification = specification[assembly]

        if not isinstance(genome_specification, dict):
            raise ResourceSpecificationError(
                f"The '{assembly}' reference-genome definition must be a JSON object."
            )

        # Required keys
        for field in ("url", "filename"):
            if field not in genome_specification:
                raise ResourceSpecificationError(
                    f"Resource section 'reference_genomes.{assembly}' "
                    f"is missing required field: {field}"
                )

        url = genome_specification["url"]
        filename = genome_specification["filename"]

        # Both absent -> valid (resource not yet defined)
        if url is None and filename is None:
            continue

        # One present but the other absent -> invalid
        if (url is None) != (filename is None):
            raise ResourceSpecificationError(
                f"Resource section 'reference_genomes.{assembly}' "
                "must define both 'url' and 'filename' or neither."
            )

        if not isinstance(url, str) or not url.strip():
            raise ResourceSpecificationError(
                f"Invalid url for reference_genomes.{assembly}"
            )

        if not isinstance(filename, str) or not filename.strip():
            raise ResourceSpecificationError(
                f"Invalid filename for reference_genomes.{assembly}"
            )


def _require_non_empty_string_fields(
        *,
        section_name: str,
        specification: dict[str, Any],
        required_fields: set[str],
) -> None:
    missing_fields = sorted(
        required_fields.difference(specification)
    )

    if missing_fields:
        raise ResourceSpecificationError(
            f"Resource section '{section_name}' is missing required "
            f"fields: {', '.join(missing_fields)}"
        )

    invalid_fields = sorted(
        field_name
        for field_name in required_fields
        if not isinstance(specification[field_name], str)
        or not specification[field_name].strip()
    )

    if invalid_fields:
        raise ResourceSpecificationError(
            f"Resource section '{section_name}' contains invalid or "
            f"empty fields: {', '.join(invalid_fields)}"
        )