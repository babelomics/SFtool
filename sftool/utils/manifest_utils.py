from __future__ import annotations

from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from sftool.utils.resource_utils import (
    ResourceOperationError,
    calculate_sha256,
    read_json,
    write_json,
)


INSTALLED_MANIFEST_SCHEMA_VERSION = 1
INSTALLED_MANIFEST_FILENAME = "resources.json"


def get_sftool_version() -> str:
    try:
        return version("sftool")
    except PackageNotFoundError:
        return "unknown"


def utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )



def relative_resource_path(
        path: Path | str,
        output_root: Path,
) -> str:
    root = output_root.resolve()
    resource_path = Path(path).resolve()

    try:
        relative_path = resource_path.relative_to(root)
    except ValueError as error:
        raise ResourceOperationError(
            f"Resource path is outside the installation root: "
            f"{resource_path}"
        ) from error

    return relative_path.as_posix()


def describe_installed_file(
        path: Path | str,
        *,
        output_root: Path,
) -> dict[str, str]:
    absolute_path = Path(path).resolve()
    root = output_root.resolve()

    try:
        relative_path = absolute_path.relative_to(root)
    except ValueError as exc:
        raise ResourceOperationError(
            f"Resource file is outside the installation root: "
            f"{absolute_path}"
        ) from exc

    if not absolute_path.is_file():
        raise ResourceOperationError(
            f"Installed resource file does not exist: {absolute_path}"
        )

    return {
        "path": relative_path.as_posix(),
        "sha256": calculate_sha256(absolute_path),
    }


def build_catalog_manifest(
        catalog_resources: dict[str, Any],
        *,
        output_root: Path,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {}

    for catalog_name, assemblies in catalog_resources.items():
        manifest[catalog_name] = {}

        for assembly, files in assemblies.items():
            manifest[catalog_name][assembly] = {
                "bed": describe_installed_file(
                    files["bed"],
                    output_root=output_root,
                ),
                "json": describe_installed_file(
                    files["json"],
                    output_root=output_root,
                ),
            }

    return manifest

def build_clinvar_manifest(
        clinvar_resources: dict[str, Any],
        *,
        output_root: Path,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {
        "version": clinvar_resources["version"],
        "files": {},
    }

    for catalog_name, assemblies in clinvar_resources["files"].items():
        manifest["files"][catalog_name] = {}

        for assembly, evidence_files in assemblies.items():
            manifest["files"][catalog_name][assembly] = {}

            for evidence_level, path in evidence_files.items():
                manifest["files"][catalog_name][assembly][evidence_level] = (
                    describe_installed_file(
                        path,
                        output_root=output_root,
                    )
                )

    return manifest

def build_reference_manifest(
        reference_resources: dict[str, Any],
        *,
        output_root: Path,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {}

    for assembly, files in reference_resources.items():
        manifest[assembly] = {
            file_type: describe_installed_file(
                path,
                output_root=output_root,
            )
            for file_type, path in files.items()
        }

    return manifest

def build_installed_manifest(
        *,
        output_root: Path,
        resource_version: str,
        catalog_resources: dict[str, Any],
        clinvar_resources: dict[str, Any],
        hpo_resource: dict[str, Any],
        pharmcat_resource: dict[str, Any],
        reference_resources: dict[str, Any] | None = None,
        existing_manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    now = utc_timestamp()

    created_at = (
        existing_manifest.get("created_at", now)
        if existing_manifest
        else now
    )

    datasets: dict[str, Any] = {
        "catalogs": build_catalog_manifest(
            catalog_resources,
            output_root=output_root,
        ),
        "clinvar": build_clinvar_manifest(
            clinvar_resources,
            output_root=output_root,
        ),
        "hpo": {
            "version": hpo_resource["version"],
            "file": describe_installed_file(
                hpo_resource["path"],
                output_root=output_root,
            ),
        },
        "pharmcat": {
            "version": pharmcat_resource["version"],
            "positions_vcf": describe_installed_file(
                pharmcat_resource["positions_vcf"],
                output_root=output_root,
            ),
        },
    }

    if reference_resources:
        datasets["reference_genomes"] = build_reference_manifest(
            reference_resources,
            output_root=output_root,
        )

    return {
        "schema_version": INSTALLED_MANIFEST_SCHEMA_VERSION,
        "resource_version": resource_version,
        "sftool_version": get_sftool_version(),
        "created_at": created_at,
        "updated_at": now,
        "datasets": datasets,
    }

def validate_existing_manifest_compatibility(
        existing_manifest: dict[str, Any],
        *,
        resource_version: str,
) -> None:
    existing_schema_version = existing_manifest.get(
        "schema_version"
    )

    if existing_schema_version != INSTALLED_MANIFEST_SCHEMA_VERSION:
        raise ResourceOperationError(
            "The existing resource manifest uses an unsupported "
            f"schema version: {existing_schema_version!r}"
        )

    existing_resource_version = existing_manifest.get(
        "resource_version"
    )

    if existing_resource_version != resource_version:
        raise ResourceOperationError(
            "The resource directory was installed using resource "
            f"version {existing_resource_version!r}, but the requested "
            f"version is {resource_version!r}."
        )

def write_installed_manifest(
        *,
        output_root: Path,
        resource_version: str,
        bundled_resources: dict[str, Any],
        catalog_resources: dict[str, Any],
        clinvar_resources: dict[str, Any],
        hpo_resource: dict[str, Any],
        pharmcat_resource: dict[str, Any],
        reference_resources: dict[str, Any] | None = None,
) -> Path:
    manifest_path = output_root / "resources.json"

    existing_manifest: dict[str, Any] | None = None

    if manifest_path.exists():
        existing_manifest = read_json(manifest_path)

        validate_existing_manifest_compatibility(
            existing_manifest,
            resource_version=resource_version,
        )

    manifest = build_installed_manifest(
        output_root=output_root,
        resource_version=resource_version,
        bundled_resources=bundled_resources,
        catalog_resources=catalog_resources,
        clinvar_resources=clinvar_resources,
        hpo_resource=hpo_resource,
        pharmcat_resource=pharmcat_resource,
        reference_resources=reference_resources,
        existing_manifest=existing_manifest,
    )

    return write_json(
        manifest,
        manifest_path,
    )