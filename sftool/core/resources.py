from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from sftool.utils.errors import ResourceManifestError

from sftool.utils.checksums import (
    ChecksumError,
    calculate_sha256,
)

SUPPORTED_MANIFEST_SCHEMA_VERSIONS = {2}
SUPPORTED_ASSEMBLIES = {"GRCh37", "GRCh38"}
SUPPORTED_CATALOGS = {"PR", "RR"}
SUPPORTED_CLINVAR_EVIDENCE_LEVELS = {1, 2, 3, 4}

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")

class RuntimeResources:
    """
    Installed resources resolved for one SFtool execution.

    The original manifest remains available for provenance, while resource
    paths are exposed as pathlib.Path objects.
    """

    def __init__(
            self,
            *,
            manifest_path: Path,
            manifest: dict[str, Any],
            installed: dict[str, Any],
            execution: dict[str, Any],
    ):
        self.manifest_path = manifest_path
        self.manifest = manifest
        self.installed = installed
        self.execution = execution

    @property
    def root(self) -> Path:
        return self.installed["root"]

    @property
    def reference_genome(self) -> Path:
        return self.execution["reference_genome"]["fasta"]

    @property
    def reference_genome_index(self) -> Path:
        return self.execution["reference_genome"]["fai"]

    @property
    def reference_genome_source(self) -> str:
        return self.execution["reference_genome"]["source"]

    @property
    def hpo_file(self) -> Path:
        return self.execution["hpo"]["file"]["path"]

    @property
    def hpo_version(self) -> str:
        return self.execution["hpo"]["version"]

    @property
    def pharmcat_positions_vcf(self) -> Path:
        return self.execution[
            "pharmcat"
        ]["positions_vcf"]["path"]

    @property
    def pharmcat_version(self) -> str:
        return self.execution["pharmcat"]["version"]

    @property
    def rr_str_catalog(self) -> Path | None:
        rr_str = self.execution.get("rr_str")

        if rr_str is None:
            return None

        return rr_str["csv"]["path"]

def load_resource_manifest(
        manifest_path: str | Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Load, validate and resolve an installed SFtool resource manifest.

    Returns:
        raw_manifest:
            Original JSON-compatible manifest.

        installed_resources:
            Runtime-normalized resource structure containing absolute Paths.
    """
    path = Path(manifest_path).expanduser().resolve()

    if not path.exists():
        raise ResourceManifestError(
            f"Resource manifest not found: {path}"
        )

    if not path.is_file():
        raise ResourceManifestError(
            f"Resource manifest is not a file: {path}"
        )

    try:
        with path.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
    except json.JSONDecodeError as exc:
        raise ResourceManifestError(
            f"Invalid JSON in resource manifest {path}: {exc}"
        ) from exc
    except OSError as exc:
        raise ResourceManifestError(
            f"Could not read resource manifest {path}: {exc}"
        ) from exc

    _validate_manifest_structure(manifest, manifest_path=path)

    root = _resolve_resources_root(
        manifest=manifest,
        manifest_path=path,
    )

    installed = _resolve_installed_resources(
        manifest=manifest,
        root=root,
    )

    return manifest, installed

def _validate_resource_checksum(
        descriptor: dict[str, Any],
        *,
        label: str,
) -> None:
    """
    Verify one resolved manifest file descriptor.
    """
    path = descriptor["path"]
    expected_checksum = descriptor["sha256"]

    try:
        actual_checksum = calculate_sha256(path)
    except ChecksumError as exc:
        raise ResourceManifestError(
            f"Could not verify installed resource "
            f"{label}: {path}. "
            "Ensure the file is readable or reinstall the "
            "resource bundle with 'sftool resources setup'."
        ) from exc

    if actual_checksum != expected_checksum:
        raise ResourceManifestError(
            f"Checksum mismatch for installed resource "
            f"{label}: {path}. "
            f"Expected SHA-256 {expected_checksum}, "
            f"found {actual_checksum}. "
            "The resource file may be corrupted or modified. "
            "Run 'sftool resources setup' again using the "
            "same resource directory."
        )

def _validate_required_resource_checksums(
        installed: dict[str, Any],
) -> None:
    catalogs = installed["catalogs"]

    _validate_resource_checksum(
        catalogs["RR_STR"]["csv"],
        label="datasets.catalogs.RR_STR.csv",
    )

    for assembly in sorted(SUPPORTED_ASSEMBLIES):
        assembly_catalogs = (
            catalogs["assemblies"][assembly]
        )

        for catalog_name in sorted(SUPPORTED_CATALOGS):
            catalog = assembly_catalogs[catalog_name]

            for file_type in ("bed", "chr_bed", "json"):
                _validate_resource_checksum(
                    catalog[file_type],
                    label=(
                        "datasets.catalogs.assemblies."
                        f"{assembly}.{catalog_name}."
                        f"{file_type}"
                    ),
                )

    clinvar = installed["clinvar"]

    for source_name in (
            "variant_summary",
            "submission_summary",
    ):
        _validate_resource_checksum(
            clinvar["source_files"][source_name],
            label=(
                "datasets.clinvar.source_files."
                f"{source_name}"
            ),
        )

    for assembly in sorted(SUPPORTED_ASSEMBLIES):
        _validate_resource_checksum(
            clinvar["databases"][assembly],
            label=(
                f"datasets.clinvar.databases.{assembly}"
            ),
        )

        for catalog_name in sorted(SUPPORTED_CATALOGS):
            evidence_resources = (
                clinvar["filtered_databases"]
                [assembly][catalog_name]
            )

            for evidence in sorted(
                    SUPPORTED_CLINVAR_EVIDENCE_LEVELS
            ):
                _validate_resource_checksum(
                    evidence_resources[evidence],
                    label=(
                        "datasets.clinvar."
                        "filtered_databases."
                        f"{assembly}.{catalog_name}."
                        f"{evidence}"
                    ),
                )

    _validate_resource_checksum(
        installed["hpo"]["file"],
        label="datasets.hpo.file",
    )

    _validate_resource_checksum(
        installed["pharmcat"]["positions_vcf"],
        label="datasets.pharmcat.positions_vcf",
    )

def _validate_optional_reference_checksums(
        installed: dict[str, Any],
) -> None:
    references = installed["reference_genomes"]

    for assembly in sorted(SUPPORTED_ASSEMBLIES):
        reference = references[assembly]

        if reference is None:
            continue

        _validate_resource_checksum(
            reference["fasta"],
            label=(
                "datasets.reference_genomes."
                f"{assembly}.fasta"
            ),
        )

        _validate_resource_checksum(
            reference["fai"],
            label=(
                "datasets.reference_genomes."
                f"{assembly}.fai"
            ),
        )

def _validate_installed_resource_checksums(
        installed: dict[str, Any],
) -> None:
    """
    Verify all files registered in an installed bundle.
    """
    _validate_required_resource_checksums(
        installed
    )
    _validate_optional_reference_checksums(
        installed
    )

def validate_resource_bundle(
        manifest_path: str | Path,
        *,
        verify_checksums: bool = True,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Load and validate an installed SFtool resource bundle.

    Validation includes:

    - manifest JSON and schema compatibility;
    - required manifest structure;
    - installation-root resolution;
    - relative-path containment;
    - presence of every manifest-referenced file;
    - SHA-256 integrity verification when enabled.

    Returns:
        A tuple containing the raw manifest and the
        runtime-normalized installed-resource structure.

    Raises:
        ResourceManifestError:
            If the resource bundle is incompatible,
            incomplete or corrupted.
    """
    manifest, installed = load_resource_manifest(
        manifest_path
    )

    if verify_checksums:
        _validate_installed_resource_checksums(
            installed
        )

    return manifest, installed

def _resolve_resources_root(
        *,
        manifest: dict[str, Any],
        manifest_path: Path,
) -> Path:
    resources_root_value = manifest.get("resources_root", ".")

    if not isinstance(resources_root_value, str):
        raise ResourceManifestError(
            "Manifest field 'resources_root' must be a string"
        )

    configured_root = Path(resources_root_value).expanduser()

    if configured_root.is_absolute():
        root = configured_root.resolve()
    else:
        root = (manifest_path.parent / configured_root).resolve()

    if not root.exists():
        raise ResourceManifestError(
            f"Resource installation root does not exist: {root}"
        )

    if not root.is_dir():
        raise ResourceManifestError(
            f"Resource installation root is not a directory: {root}"
        )

    return root

def _validate_manifest_structure(
        manifest: Any,
        *,
        manifest_path: Path,
) -> None:
    if not isinstance(manifest, dict):
        raise ResourceManifestError(
            f"Resource manifest must contain a JSON object: {manifest_path}"
        )

    schema_version = manifest.get("schema_version")

    if schema_version not in SUPPORTED_MANIFEST_SCHEMA_VERSIONS:
        supported = ", ".join(
            str(version)
            for version in sorted(SUPPORTED_MANIFEST_SCHEMA_VERSIONS)
        )

        raise ResourceManifestError(
            f"Unsupported resource manifest schema version "
            f"{schema_version!r} in {manifest_path}. "
            f"Supported version(s): {supported}."
        )

    datasets = manifest.get("datasets")

    if not isinstance(datasets, dict):
        raise ResourceManifestError(
            "Manifest field 'datasets' must be an object"
        )

    required_datasets = {
        "catalogs",
        "clinvar",
        "hpo",
        "pharmcat",
    }

    missing = required_datasets - datasets.keys()

    if missing:
        raise ResourceManifestError(
            "Resource manifest is missing required dataset block(s): "
            + ", ".join(sorted(missing))
        )

    _validate_catalog_manifest(datasets["catalogs"])
    _validate_clinvar_manifest(datasets["clinvar"])
    _validate_single_file_dataset("hpo", datasets.get("hpo"), descriptor_key="file")
    _validate_single_file_dataset("pharmcat", datasets.get("pharmcat"), descriptor_key="positions_vcf")
    _validate_reference_manifest(
        datasets.get("reference_genomes")
    )

def _validate_file_descriptor(
        descriptor: Any,
        *,
        label: str,
) -> None:
    if not isinstance(descriptor, dict):
        raise ResourceManifestError(
            f"{label} must be an object"
        )

    path_value = descriptor.get("path")

    if not isinstance(path_value, str) or not path_value.strip():
        raise ResourceManifestError(
            f"{label}.path must be a non-empty string"
        )

    checksum = descriptor.get("sha256")

    if not isinstance(checksum, str):
        raise ResourceManifestError(
            f"{label}.sha256 must be a string"
        )

    if not _SHA256_PATTERN.fullmatch(checksum):
        raise ResourceManifestError(
            f"{label}.sha256 is not a valid lowercase SHA-256 digest"
        )

    source_url = descriptor.get("source_url")

    if source_url is not None and (
            not isinstance(source_url, str)
            or not source_url.strip()
    ):
        raise ResourceManifestError(
            f"{label}.source_url must be a non-empty string when provided"
        )

def _resolve_descriptor(
        descriptor: dict[str, Any],
        *,
        root: Path,
        label: str,
        require_file: bool = True,
) -> dict[str, Any]:
    _validate_file_descriptor(
        descriptor,
        label=label,
    )

    relative_path = Path(descriptor["path"])

    if relative_path.is_absolute():
        raise ResourceManifestError(
            f"{label}.path must be relative to the resource root: "
            f"{relative_path}"
        )

    absolute_path = (root / relative_path).resolve()

    try:
        absolute_path.relative_to(root)
    except ValueError as exc:
        raise ResourceManifestError(
            f"{label}.path escapes the resource installation root: "
            f"{relative_path}"
        ) from exc

    if require_file:
        if not absolute_path.exists():
            raise ResourceManifestError(
                f"Installed resource file does not exist for {label}: "
                f"{absolute_path}"
            )

        if not absolute_path.is_file():
            raise ResourceManifestError(
                f"Installed resource path is not a file for {label}: "
                f"{absolute_path}"
            )

    resolved = dict(descriptor)
    resolved["path"] = absolute_path

    return resolved

def _validate_catalog_manifest(catalogs: Any) -> None:
    if not isinstance(catalogs, dict):
        raise ResourceManifestError(
            "datasets.catalogs must be an object"
        )

    rr_str = catalogs.get("RR_STR")

    if not isinstance(rr_str, dict):
        raise ResourceManifestError(
            "datasets.catalogs.RR_STR must be an object"
        )

    _validate_version(
        rr_str.get("version"),
        label="datasets.catalogs.RR_STR.version",
    )

    _validate_file_descriptor(
        rr_str.get("csv"),
        label="datasets.catalogs.RR_STR.csv",
    )

    assemblies = catalogs.get("assemblies")

    if not isinstance(assemblies, dict):
        raise ResourceManifestError(
            "datasets.catalogs.assemblies must be an object"
        )

    for assembly in SUPPORTED_ASSEMBLIES:
        assembly_catalogs = assemblies.get(assembly)

        if not isinstance(assembly_catalogs, dict):
            raise ResourceManifestError(
                f"datasets.catalogs.assemblies.{assembly} "
                "must be an object"
            )

        for catalog_name in SUPPORTED_CATALOGS:
            catalog = assembly_catalogs.get(catalog_name)

            if not isinstance(catalog, dict):
                raise ResourceManifestError(
                    f"datasets.catalogs.assemblies."
                    f"{assembly}.{catalog_name} must be an object"
                )

            _validate_version(
                catalog.get("version"),
                label=(
                    f"datasets.catalogs.assemblies."
                    f"{assembly}.{catalog_name}.version"
                ),
            )

            for file_type in ("bed", "chr_bed", "json"):
                _validate_file_descriptor(
                    catalog.get(file_type),
                    label=(
                        f"datasets.catalogs.assemblies."
                        f"{assembly}.{catalog_name}.{file_type}"
                    ),
                )

def _validate_version(value: Any, *, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ResourceManifestError(
            f"{label} must be a non-empty string"
        )

def _resolve_catalogs(
        catalogs: dict[str, Any],
        *,
        root: Path,
) -> dict[str, Any]:
    resolved: dict[str, Any] = {
        "RR_STR": {
            "version": catalogs["RR_STR"]["version"],
            "csv": _resolve_descriptor(
                catalogs["RR_STR"]["csv"],
                root=root,
                label="datasets.catalogs.RR_STR.csv",
            ),
        },
        "assemblies": {},
    }

    for assembly, assembly_catalogs in catalogs["assemblies"].items():
        resolved["assemblies"][assembly] = {}

        for catalog_name, catalog in assembly_catalogs.items():
            resolved_catalog = {
                "version": catalog["version"],
            }

            for file_type in ("bed", "chr_bed", "json"):
                resolved_catalog[file_type] = _resolve_descriptor(
                    catalog[file_type],
                    root=root,
                    label=(
                        f"datasets.catalogs.assemblies."
                        f"{assembly}.{catalog_name}.{file_type}"
                    ),
                )

            resolved["assemblies"][assembly][catalog_name] = (
                resolved_catalog
            )

    return resolved

def _validate_clinvar_manifest(clinvar: Any) -> None:
    if not isinstance(clinvar, dict):
        raise ResourceManifestError(
            "datasets.clinvar must be an object"
        )

    _validate_version(
        clinvar.get("version"),
        label="datasets.clinvar.version",
    )

    source_files = clinvar.get("source_files")

    if not isinstance(source_files, dict):
        raise ResourceManifestError(
            "datasets.clinvar.source_files must be an object"
        )

    for source_name in (
            "variant_summary",
            "submission_summary",
    ):
        _validate_file_descriptor(
            source_files.get(source_name),
            label=f"datasets.clinvar.source_files.{source_name}",
        )

    databases = clinvar.get("databases")

    if not isinstance(databases, dict):
        raise ResourceManifestError(
            "datasets.clinvar.databases must be an object"
        )

    filtered_databases = clinvar.get("filtered_databases")

    if not isinstance(filtered_databases, dict):
        raise ResourceManifestError(
            "datasets.clinvar.filtered_databases must be an object"
        )

    for assembly in SUPPORTED_ASSEMBLIES:
        _validate_file_descriptor(
            databases.get(assembly),
            label=f"datasets.clinvar.databases.{assembly}",
        )

        assembly_filtered_databases = filtered_databases.get(assembly)

        if not isinstance(assembly_filtered_databases, dict):
            raise ResourceManifestError(
                f"datasets.clinvar.filtered_databases."
                f"{assembly} must be an object"
            )

        for catalog_name in SUPPORTED_CATALOGS:
            evidence_files = assembly_filtered_databases.get(
                catalog_name
            )

            if not isinstance(evidence_files, dict):
                raise ResourceManifestError(
                    f"datasets.clinvar.filtered_databases."
                    f"{assembly}.{catalog_name} must be an object"
                )

            actual_levels = set()

            for level_key, descriptor in evidence_files.items():
                try:
                    level = int(level_key)
                except (TypeError, ValueError) as exc:
                    raise ResourceManifestError(
                        f"Invalid ClinVar evidence key "
                        f"{level_key!r} for "
                        f"{assembly}/{catalog_name}"
                    ) from exc

                actual_levels.add(level)

                _validate_file_descriptor(
                    descriptor,
                    label=(
                        f"datasets.clinvar.filtered_databases."
                        f"{assembly}.{catalog_name}.{level_key}"
                    ),
                )

            if actual_levels != SUPPORTED_CLINVAR_EVIDENCE_LEVELS:
                raise ResourceManifestError(
                    f"ClinVar resources for "
                    f"{assembly}/{catalog_name} must provide "
                    f"evidence levels "
                    f"{sorted(SUPPORTED_CLINVAR_EVIDENCE_LEVELS)}; "
                    f"found {sorted(actual_levels)}"
                )


def _resolve_clinvar(
        clinvar: dict[str, Any],
        *,
        root: Path,
) -> dict[str, Any]:
    resolved: dict[str, Any] = {
        "version": clinvar["version"],
        "source_files": {},
        "databases": {},
        "filtered_databases": {},
    }

    for source_name, descriptor in clinvar["source_files"].items():
        resolved["source_files"][source_name] = _resolve_descriptor(
            descriptor,
            root=root,
            label=f"datasets.clinvar.source_files.{source_name}",
        )

    for assembly, descriptor in clinvar["databases"].items():
        resolved["databases"][assembly] = _resolve_descriptor(
            descriptor,
            root=root,
            label=f"datasets.clinvar.databases.{assembly}",
        )

    for assembly, assembly_databases in (
            clinvar["filtered_databases"].items()
    ):
        resolved["filtered_databases"][assembly] = {}

        for catalog_name, evidence_files in (
                assembly_databases.items()
        ):
            resolved["filtered_databases"][assembly][catalog_name] = {
                int(level): _resolve_descriptor(
                    descriptor,
                    root=root,
                    label=(
                        f"datasets.clinvar.filtered_databases."
                        f"{assembly}.{catalog_name}.{level}"
                    ),
                )
                for level, descriptor in evidence_files.items()
            }

    return resolved

def _validate_single_file_dataset(
        dataset_name: str,
        dataset: Any,
        *,
        descriptor_key: str,
) -> None:
    if not isinstance(dataset, dict):
        raise ResourceManifestError(
            f"datasets.{dataset_name} must be an object"
        )

    _validate_version(
        dataset.get("version"),
        label=f"datasets.{dataset_name}.version",
    )

    _validate_file_descriptor(
        dataset.get(descriptor_key),
        label=f"datasets.{dataset_name}.{descriptor_key}",
    )

def _resolve_hpo(
        hpo: dict[str, Any],
        *,
        root: Path,
) -> dict[str, Any]:
    return {
        "version": hpo["version"],
        "file": _resolve_descriptor(
            hpo["file"],
            root=root,
            label="datasets.hpo.file",
        ),
    }

def _resolve_pharmcat(
        pharmcat: dict[str, Any],
        *,
        root: Path,
) -> dict[str, Any]:
    return {
        "version": pharmcat["version"],
        "positions_vcf": _resolve_descriptor(
            pharmcat["positions_vcf"],
            root=root,
            label="datasets.pharmcat.positions_vcf",
        ),
    }

def _validate_reference_manifest(
        references: Any,
) -> None:
    if references is None:
        return

    if not isinstance(references, dict):
        raise ResourceManifestError(
            "datasets.reference_genomes must be "
            "an object or null"
        )

    unexpected = set(references) - SUPPORTED_ASSEMBLIES

    if unexpected:
        raise ResourceManifestError(
            "Unsupported reference genome assembly key(s): "
            + ", ".join(sorted(unexpected))
        )

    for assembly in SUPPORTED_ASSEMBLIES:
        reference = references.get(assembly)

        if reference is None:
            continue

        if not isinstance(reference, dict):
            raise ResourceManifestError(
                f"datasets.reference_genomes.{assembly} "
                "must be an object or null"
            )

        _validate_file_descriptor(
            reference.get("fasta"),
            label=(
                f"datasets.reference_genomes."
                f"{assembly}.fasta"
            ),
        )

        _validate_file_descriptor(
            reference.get("fai"),
            label=(
                f"datasets.reference_genomes."
                f"{assembly}.fai"
            ),
        )

def _resolve_reference_genomes(
        references: dict[str, Any] | None,
        *,
        root: Path,
) -> dict[str, dict[str, Any] | None]:
    references = references or {}

    resolved: dict[str, dict[str, Any] | None] = {}

    for assembly in SUPPORTED_ASSEMBLIES:
        reference = references.get(assembly)

        if reference is None:
            resolved[assembly] = None
            continue

        resolved[assembly] = {
            "fasta": _resolve_descriptor(
                reference["fasta"],
                root=root,
                label=(
                    f"datasets.reference_genomes."
                    f"{assembly}.fasta"
                ),
            ),
            "fai": _resolve_descriptor(
                reference["fai"],
                root=root,
                label=(
                    f"datasets.reference_genomes."
                    f"{assembly}.fai"
                ),
            ),
        }

    return resolved

def _resolve_installed_resources(
        *,
        manifest: dict[str, Any],
        root: Path,
) -> dict[str, Any]:
    datasets = manifest["datasets"]

    return {
        "root": root,
        "catalogs": _resolve_catalogs(
            datasets["catalogs"],
            root=root,
        ),
        "clinvar": _resolve_clinvar(
            datasets["clinvar"],
            root=root,
        ),
        "hpo": _resolve_hpo(
            datasets["hpo"],
            root=root,
        ),
        "pharmcat": _resolve_pharmcat(
            datasets["pharmcat"],
            root=root,
        ),
        "reference_genomes": _resolve_reference_genomes(
            datasets.get("reference_genomes"),
            root=root,
        ),
    }


def _validate_user_reference(
        fasta_path: Path,
        *,
        assembly: str,
) -> dict[str, Any]:
    fasta = fasta_path.expanduser().resolve()

    if not fasta.exists():
        raise ResourceManifestError(
            f"Configured {assembly} reference genome "
            f"does not exist: {fasta}"
        )

    if not fasta.is_file():
        raise ResourceManifestError(
            f"Configured {assembly} reference genome "
            f"is not a file: {fasta}"
        )

    fai = Path(f"{fasta}.fai")

    if not fai.exists():
        raise ResourceManifestError(
            f"FASTA index not found for configured {assembly} "
            f"reference genome: {fai}. "
            "Index the FASTA with 'samtools faidx' or provide "
            "a reference installed by 'sftool resources setup'."
        )

    if not fai.is_file():
        raise ResourceManifestError(
            f"Configured {assembly} FASTA index "
            f"is not a file: {fai}"
        )

    return {
        "fasta": fasta,
        "fai": fai.resolve(),
        "source": "config",
        "fasta_sha256": None,
        "fai_sha256": None,
    }

def _resolve_execution_reference(
        *,
        assembly: str,
        configured_genomes: dict[str, Path | None],
        installed_references: dict[
            str,
            dict[str, Any] | None
        ],
) -> dict[str, Any]:
    configured_reference = configured_genomes.get(assembly)

    if configured_reference is not None:
        return _validate_user_reference(
            configured_reference,
            assembly=assembly,
        )

    installed_reference = installed_references.get(assembly)

    if installed_reference is not None:
        return {
            "fasta": installed_reference["fasta"]["path"],
            "fai": installed_reference["fai"]["path"],
            "source": "manifest",
            "fasta_sha256": (
                installed_reference["fasta"]["sha256"]
            ),
            "fai_sha256": (
                installed_reference["fai"]["sha256"]
            ),
        }

    raise ResourceManifestError(
        f"No reference genome is available for {assembly}. "
        f"Set references.genomes.{assembly} in the runtime "
        "configuration or install reference genomes with "
        "'sftool resources setup --download-reference-genomes'."
    )

def resolve_execution_resources(
        *,
        assembly: str,
        clinvar_evidence: int,
        categories: set[str],
        configured_genomes: dict[str, Path | None],
        installed: dict[str, Any],
) -> dict[str, Any]:
    if assembly not in SUPPORTED_ASSEMBLIES:
        raise ResourceManifestError(
            f"Unsupported assembly: {assembly}"
        )

    if (
            not isinstance(clinvar_evidence, int)
            or isinstance(clinvar_evidence, bool)
            or clinvar_evidence
            not in SUPPORTED_CLINVAR_EVIDENCE_LEVELS
    ):
        supported = ", ".join(
            str(level)
            for level in sorted(
                SUPPORTED_CLINVAR_EVIDENCE_LEVELS
            )
        )

        raise ResourceManifestError(
            f"Unsupported ClinVar evidence level: "
            f"{clinvar_evidence!r}. "
            f"Supported levels are: {supported}."
        )

    unsupported_categories = (
            set(categories) - SUPPORTED_CATALOGS
    )

    if unsupported_categories:
        raise ResourceManifestError(
            "Unsupported resource category or categories: "
            + ", ".join(sorted(unsupported_categories))
        )

    assembly_catalogs = (
        installed["catalogs"]
        ["assemblies"]
        .get(assembly)
    )

    if assembly_catalogs is None:
        raise ResourceManifestError(
            f"No installed catalogs are available for "
            f"{assembly}"
        )

    clinvar_filtered_databases = (
        installed["clinvar"]
        ["filtered_databases"]
        .get(assembly)
    )

    if clinvar_filtered_databases is None:
        raise ResourceManifestError(
            f"No installed filtered ClinVar resources "
            f"are available for {assembly}"
        )

    selected_catalogs: dict[str, Any] = {}
    selected_clinvar: dict[str, Any] = {}

    for category in sorted(categories):
        catalog = assembly_catalogs.get(category)

        if catalog is None:
            raise ResourceManifestError(
                f"No installed {category} catalog is "
                f"available for {assembly}"
            )

        selected_catalogs[category] = catalog

        evidence_files = (
            clinvar_filtered_databases.get(category)
        )

        if evidence_files is None:
            raise ResourceManifestError(
                f"No installed ClinVar resources are "
                f"available for assembly={assembly}, "
                f"catalog={category}"
            )

        clinvar_file = evidence_files.get(
            clinvar_evidence
        )

        if clinvar_file is None:
            raise ResourceManifestError(
                f"No installed ClinVar resource is "
                f"available for assembly={assembly}, "
                f"catalog={category}, "
                f"evidence={clinvar_evidence}"
            )

        selected_clinvar[category] = clinvar_file

    rr_str = None

    if "RR" in categories:
        rr_str = installed["catalogs"].get("RR_STR")

        if rr_str is None:
            raise ResourceManifestError(
                "The installed resource bundle does not "
                "contain the RR-STR catalog required for "
                "RR analysis"
            )

    return {
        "assembly": assembly,
        "reference_genome": _resolve_execution_reference(
            assembly=assembly,
            configured_genomes=configured_genomes,
            installed_references=(
                installed["reference_genomes"]
            ),
        ),
        "catalogs": selected_catalogs,
        "clinvar": selected_clinvar,
        "rr_str": rr_str,
        "hpo": installed["hpo"],
        "pharmcat": installed["pharmcat"],
    }

def get_required_resource_categories(
        samples_info: dict[str, Any],
) -> set[str]:
    execution = samples_info["execution"]

    if "secondary_findings_discovery" not in execution["modes"]:
        return set()

    required = {
        category
        for sample in samples_info["samples"]
        for category in sample.get("categories", [])
        if category in SUPPORTED_CATALOGS
    }

    return required

def clinvar_file(
        self,
        category: str,
) -> Path:
    descriptor = self.execution[
        "clinvar"
    ].get(category)

    if descriptor is None:
        raise ResourceManifestError(
            f"No ClinVar resource was selected for "
            f"category {category}"
        )

    return descriptor["path"]

def catalog_json(
        self,
        category: str,
) -> Path:
    catalog = self.execution[
        "catalogs"
    ].get(category)

    if catalog is None:
        raise ResourceManifestError(
            f"No catalog resource was selected for "
            f"category {category}"
        )

    return catalog["json"]["path"]

def catalog_bed(
        self,
        category: str,
) -> Path:
    catalog = self.execution[
        "catalogs"
    ].get(category)

    if catalog is None:
        raise ResourceManifestError(
            f"No catalog resource was selected for "
            f"category {category}"
        )

    return catalog["bed"]["path"]

def catalog_chr_bed(
        self,
        category: str,
) -> Path:
    catalog = self.execution[
        "catalogs"
    ].get(category)

    if catalog is None:
        raise ResourceManifestError(
            f"No catalog resource was selected for "
            f"category {category}"
        )

    return catalog["chr_bed"]["path"]