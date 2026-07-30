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


INSTALLED_MANIFEST_SCHEMA_VERSION = 2
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
        source_url: str | None = None,
) -> dict[str, str]:
    root = output_root.resolve()
    resource_path = Path(path)

    if resource_path.is_absolute():
        absolute_path = resource_path.resolve()
    else:
        absolute_path = (root / resource_path).resolve()

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

    descriptor = {
        "path": relative_path.as_posix(),
        "sha256": calculate_sha256(absolute_path),
    }

    if source_url is not None:
        descriptor["source_url"] = source_url

    return descriptor

def build_catalog_manifest(
        catalog_resources: dict[str, Any],
        *,
        output_root: Path,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {
        "assemblies": {},
    }

    assemblies = catalog_resources.get("assemblies", {})

    for assembly, catalogs in assemblies.items():
        manifest["assemblies"][assembly] = {}

        for catalog_name, resources in catalogs.items():
            manifest["assemblies"][assembly][catalog_name] = {
                "version": resources["version"],
                "bed": describe_installed_file(
                    resources["bed"],
                    output_root=output_root,
                ),
                "chr_bed": describe_installed_file(
                    resources["chr_bed"],
                    output_root=output_root,
                ),
                "json": describe_installed_file(
                    resources["json"],
                    output_root=output_root,
                ),
            }

    if "RR_STR" in catalog_resources:
        manifest["RR_STR"] = {
            "version": catalog_resources["RR_STR"]["version"],
            "csv": describe_installed_file(
                catalog_resources["RR_STR"]["csv"],
                output_root=output_root,
            ),
        }

    return manifest

def describe_file_tree(
        resources: dict[str, Any],
        *,
        output_root: Path,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {}

    for key, value in resources.items():
        if isinstance(value, dict):
            manifest[key] = describe_file_tree(
                value,
                output_root=output_root,
            )
        else:
            manifest[key] = describe_installed_file(
                value,
                output_root=output_root,
            )

    return manifest

def build_clinvar_manifest(
        clinvar_resources: dict[str, Any],
        clinvar_databases: dict[str, Any],
        filtered_clinvar_databases: dict[str, Any],
        *,
        output_root: Path,
) -> dict[str, Any]:
    return {
        "version": clinvar_resources["version"],
        "source_files": {
            "variant_summary": describe_installed_file(
                clinvar_resources["variant_summary"],
                output_root=output_root,
                source_url=clinvar_resources[
                    "variant_summary_url"
                ],
            ),
            "submission_summary": describe_installed_file(
                clinvar_resources["submission_summary"],
                output_root=output_root,
                source_url=clinvar_resources[
                    "submission_summary_url"
                ],
            ),
        },
        "databases": describe_file_tree(
            clinvar_databases,
            output_root=output_root,
        ),
        "filtered_databases": describe_file_tree(
            filtered_clinvar_databases,
            output_root=output_root,
        ),
    }

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
        resource_version_policy: str,
        catalog_resources: dict[str, Any],
        clinvar_resources: dict[str, Any],
        clinvar_databases: dict[str, Any],
        filtered_clinvar_databases: dict[str, Any],
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
            clinvar_databases,
            filtered_clinvar_databases,
            output_root=output_root,
        ),
        "hpo": {
            "version": hpo_resource["version"],
            "file": describe_installed_file(
                hpo_resource["path"],
                output_root=output_root,
                source_url=hpo_resource["source_url"],
            ),
        },
        "pharmcat": {
            "version": pharmcat_resource["version"],
            "positions_vcf": describe_installed_file(
                pharmcat_resource["path"],
                output_root=output_root,
                source_url=pharmcat_resource["source_url"],
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
        "resource_version_policy": resource_version_policy,
        "sftool_version": get_sftool_version(),
        "created_at": created_at,
        "updated_at": now,
        "datasets": datasets,
    }

def validate_existing_manifest_compatibility(
        existing_manifest: dict[str, Any],
        *,
        resource_version_policy: str,
) -> None:
    existing_schema_version = existing_manifest.get(
        "schema_version"
    )

    if existing_schema_version != INSTALLED_MANIFEST_SCHEMA_VERSION:
        raise ResourceOperationError(
            "The existing resource manifest uses an unsupported "
            f"schema version: {existing_schema_version!r}"
        )

    existing_resource_version_policy = existing_manifest.get(
        "resource_version_policy"
    )

    if (
            existing_resource_version_policy
            != resource_version_policy
    ):
        raise ResourceOperationError(
            "The resource directory was installed using resource "
            f"version policy "
            f"{existing_resource_version_policy!r}, but the requested "
            f"policy is {resource_version_policy!r}."
        )

def write_installed_manifest(
        *,
        output_root: Path,
        resource_version_policy: str,
        catalog_resources: dict[str, Any],
        clinvar_resources: dict[str, Any],
        clinvar_databases: dict[str, Any],
        filtered_clinvar_databases: dict[str, Any],
        hpo_resource: dict[str, Any],
        pharmcat_resource: dict[str, Any],
        reference_resources: dict[str, Any] | None = None,
) -> Path:
    manifest_path = output_root / INSTALLED_MANIFEST_FILENAME

    existing_manifest: dict[str, Any] | None = None

    if manifest_path.exists():
        existing_manifest = read_json(manifest_path)

        validate_existing_manifest_compatibility(
            existing_manifest,
            resource_version_policy=resource_version_policy,
        )

    manifest = build_installed_manifest(
        output_root=output_root,
        resource_version_policy=resource_version_policy,
        catalog_resources=catalog_resources,
        clinvar_resources=clinvar_resources,
        clinvar_databases=clinvar_databases,
        filtered_clinvar_databases=filtered_clinvar_databases,
        hpo_resource=hpo_resource,
        pharmcat_resource=pharmcat_resource,
        reference_resources=reference_resources,
        existing_manifest=existing_manifest,
    )

    return write_json(
        manifest,
        manifest_path,
    )