"""
Utilities for loading and validating SFtool resource specifications.
"""

from __future__ import annotations

import json
from importlib import resources
from typing import Any
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import hashlib
from collections.abc import Mapping, Callable
import shutil


BUNDLED_RESOURCE_PACKAGE = "sftool.data.resources"
BUNDLED_RESOURCE_FILENAME = "bundled_resources.json"
SUPPORTED_BUNDLED_SCHEMA_VERSIONS = {1}


EXPECTED_HPO_COLUMNS = [
    "ncbi_gene_id",
    "gene_symbol",
    "hpo_id",
    "hpo_name",
    "frequency",
    "disease_id",
]

class ResourceSpecificationError(ValueError):
    """
    Raised when a bundled resource specification is missing or invalid.
    """

class ResourceOperationError(RuntimeError):
    """
    Raised when a resource filesystem or download operation fails.
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
            "installed_filename",
        },
    )

    render_resource_filename(
        specification["installed_filename"],
        version=specification["version"],
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

from pathlib import Path


def ensure_directory(path: Path) -> Path:
    """
    Create a directory and its parents if necessary.

    Parameters
    ----------
    path
        Directory to create.

    Returns
    -------
    Path
        The normalized directory path.

    Raises
    ------
    ResourceOperationError
        If the path exists but is not a directory, or if directory creation
        fails.
    """

    directory = Path(path)

    if directory.exists() and not directory.is_dir():
        raise ResourceOperationError(
            f"Resource directory path exists but is not a directory: "
            f"{directory}"
        )

    try:
        directory.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise ResourceOperationError(
            f"Could not create resource directory: {directory}"
        ) from error

    return directory

def render_resource_filename(
        filename_template: str,
        *,
        version: str,
) -> str:
    """
    Render a resource filename template using its version.

    Parameters
    ----------
    filename_template
        Filename template, potentially containing ``{version}``.
    version
        Resource version used to render the template.

    Returns
    -------
    str
        Rendered filename.

    Raises
    ------
    ResourceSpecificationError
        If the template cannot be rendered or produces an invalid filename.
    """

    if not isinstance(filename_template, str) or not filename_template.strip():
        raise ResourceSpecificationError(
            "Resource filename template must be a non-empty string."
        )

    if not isinstance(version, str) or not version.strip():
        raise ResourceSpecificationError(
            "Resource version must be a non-empty string."
        )

    try:
        filename = filename_template.format(version=version)
    except (KeyError, ValueError, IndexError) as error:
        raise ResourceSpecificationError(
            f"Could not render resource filename template: "
            f"{filename_template!r}"
        ) from error

    if not filename.strip():
        raise ResourceSpecificationError(
            "Rendered resource filename must not be empty."
        )

    if Path(filename).name != filename:
        raise ResourceSpecificationError(
            f"Rendered resource filename must not contain a directory path: "
            f"{filename!r}"
        )

    return filename


from urllib.parse import urljoin


def resolve_resource_url(
        base_url: str,
        filename: str,
) -> str:
    """
    Build the complete URL for a resource file.

    Parameters
    ----------
    base_url
        Base URL containing the resource.
    filename
        Resource filename.

    Returns
    -------
    str
        Complete resource URL.

    Raises
    ------
    ResourceSpecificationError
        If either argument is invalid.
    """

    if not isinstance(base_url, str) or not base_url.strip():
        raise ResourceSpecificationError(
            "Resource base URL must be a non-empty string."
        )

    if not isinstance(filename, str) or not filename.strip():
        raise ResourceSpecificationError(
            "Resource filename must be a non-empty string."
        )

    normalized_base_url = base_url.rstrip("/") + "/"

    return urljoin(normalized_base_url, filename)


def download_file(
        url: str,
        destination: Path,
        *,
        overwrite: bool = False,
        chunk_size: int = 1024 * 1024,
        timeout: float = 60.0,
) -> Path:
    """
    Download a resource to a local file.

    The file is first written to a temporary ``.part`` path and is moved to
    the final destination only after the download succeeds.

    Parameters
    ----------
    url
        HTTP or HTTPS resource URL.
    destination
        Final local path.
    overwrite
        Replace an existing destination file when true.
    chunk_size
        Number of bytes copied per iteration.
    timeout
        Network timeout in seconds.

    Returns
    -------
    Path
        Final downloaded file path.

    Raises
    ------
    ResourceOperationError
        If the URL is invalid, the destination already exists, or the
        download fails.
    """

    if not isinstance(url, str) or not url.strip():
        raise ResourceOperationError(
            "Download URL must be a non-empty string."
        )

    if chunk_size <= 0:
        raise ValueError("chunk_size must be greater than zero.")

    if timeout <= 0:
        raise ValueError("timeout must be greater than zero.")

    destination = Path(destination)

    if destination.exists() and not overwrite:
        return destination

    ensure_directory(destination.parent)

    temporary_path = destination.with_name(
        f"{destination.name}.part"
    )

    if temporary_path.exists():
        try:
            temporary_path.unlink()
        except OSError as error:
            raise ResourceOperationError(
                f"Could not remove incomplete download: {temporary_path}"
            ) from error

    request = Request(
        url,
        headers={
            "User-Agent": "SFtool resource downloader",
        },
    )

    try:
        with urlopen(request, timeout=timeout) as response:
            with temporary_path.open("wb") as output_handle:
                shutil.copyfileobj(
                    response,
                    output_handle,
                    length=chunk_size,
                )

        temporary_path.replace(destination)

    except (HTTPError, URLError, TimeoutError, OSError) as error:
        try:
            temporary_path.unlink(missing_ok=True)
        except OSError:
            pass

        raise ResourceOperationError(
            f"Could not download resource from {url!r} "
            f"to {destination}"
        ) from error

    return destination

def calculate_sha256(
        path: Path,
        *,
        chunk_size: int = 1024 * 1024,
) -> str:
    """
    Calculate the SHA-256 digest of a file.

    Parameters
    ----------
    path
        File to hash.
    chunk_size
        Number of bytes read per iteration.

    Returns
    -------
    str
        Lowercase hexadecimal SHA-256 digest.

    Raises
    ------
    ResourceOperationError
        If the path does not exist, is not a file, or cannot be read.
    """

    file_path = Path(path)

    if not file_path.is_file():
        raise ResourceOperationError(
            f"Cannot calculate checksum because the resource file "
            f"does not exist: {file_path}"
        )

    if chunk_size <= 0:
        raise ValueError("chunk_size must be greater than zero.")

    digest = hashlib.sha256()

    try:
        with file_path.open("rb") as handle:
            while chunk := handle.read(chunk_size):
                digest.update(chunk)
    except OSError as error:
        raise ResourceOperationError(
            f"Could not read resource file: {file_path}"
        ) from error

    return digest.hexdigest()

def write_json(
        data: Mapping[str, Any],
        destination: Path,
) -> Path:
    """
    Write a JSON object atomically.

    Parameters
    ----------
    data
        JSON-compatible mapping.
    destination
        Output JSON path.

    Returns
    -------
    Path
        Written JSON path.

    Raises
    ------
    ResourceOperationError
        If serialization or writing fails.
    """

    destination = Path(destination)
    ensure_directory(destination.parent)

    temporary_path = destination.with_name(
        f"{destination.name}.tmp"
    )

    try:
        with temporary_path.open("w", encoding="utf-8") as handle:
            json.dump(
                data,
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")

        temporary_path.replace(destination)

    except (OSError, TypeError, ValueError) as error:
        try:
            temporary_path.unlink(missing_ok=True)
        except OSError:
            pass

        raise ResourceOperationError(
            f"Could not write JSON resource file: {destination}"
        ) from error

    return destination

def read_json(path: Path) -> dict[str, Any]:
    """
    Read a JSON object from a local file.

    Parameters
    ----------
    path
        JSON file path.

    Returns
    -------
    dict[str, Any]
        Parsed JSON object.

    Raises
    ------
    ResourceOperationError
        If the file cannot be read, contains invalid JSON, or does not contain
        a JSON object.
    """

    json_path = Path(path)

    if not json_path.is_file():
        raise ResourceOperationError(
            f"JSON resource file does not exist: {json_path}"
        )

    try:
        with json_path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ResourceOperationError(
            f"Could not read JSON resource file: {json_path}"
        ) from error

    if not isinstance(data, dict):
        raise ResourceOperationError(
            f"JSON resource file must contain an object: {json_path}"
        )

    return data

def download_versioned_resource(
        *,
        output_root: Path,
        resource_name: str,
        specification: Mapping[str, Any],
        validator: Callable[[Path], None] | None = None,
) -> dict[str, str]:
    version = specification["version"]

    filename = render_resource_filename(
        specification["installed_filename"],
        version=version,
    )

    destination = (
            Path(output_root)
            / resource_name
            / filename
    )

    path = download_file(
        url=specification["url"],
        destination=destination,
    )

    if validator is not None:
        validator(path)

    return {
        "version": version,
        "source_url": specification["url"],
        "path": path.relative_to(output_root).as_posix(),
        "sha256": calculate_sha256(path),
    }

def validate_hpo_gene_to_phenotype(path: Path) -> None:
    path = Path(path)

    try:
        with path.open("r", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\r\n")
    except OSError as error:
        raise ResourceOperationError(
            f"Could not read HPO resource: {path}"
        ) from error

    if not header:
        raise ResourceOperationError(
            f"HPO resource is empty: {path}"
        )

    actual_columns = header.split("\t")

    if actual_columns != EXPECTED_HPO_COLUMNS:
        raise ResourceOperationError(
            "HPO resource has an unexpected header. "
            f"Expected {EXPECTED_HPO_COLUMNS}, "
            f"found {actual_columns}: {path}"
        )

def validate_vcf_resource(path: Path) -> None:
    try:
        with path.open("r", encoding="utf-8") as handle:
            first_line = handle.readline().strip()
    except OSError as error:
        raise ResourceOperationError(
            f"Could not read downloaded VCF resource: {path}"
        ) from error

    if not first_line.startswith("##fileformat=VCF"):
        raise ResourceOperationError(
            f"Downloaded file is not a valid VCF resource: {path}"
        )