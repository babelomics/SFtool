"""
Runtime configuration wrappers.

This module contains executable configuration, credentials,
analysis thresholds, the installed resource manifest location,
and optional reference-genome overrides. Persistent biological
resources are resolved from the installed manifest through
RuntimeResources.
"""
from __future__ import annotations

from typing import Dict, Any
from pathlib import Path

OBSOLETE_CONFIG_FIELDS = {
    "catalogs": (
        "Catalog resources are installed by `sftool resources setup` "
        "and resolved from resources.manifest."
    ),
    "clinvar": (
        "ClinVar resources are installed by `sftool resources setup` "
        "and resolved from resources.manifest."
    ),
}

OBSOLETE_REFERENCE_FIELDS = {
    "gene_to_phenotype_file",
    "pharmcat_positions",
    "reproductive_risk_geneset_STR",
}

class Config:
    """
    Run-level configuration facade.
    """

    def __init__(self, cfg: dict[str, Any]):

        _reject_obsolete_configuration(cfg)

        self.version = cfg.get("version", {})
        self.paths = PathsConfig(cfg.get("paths", {}))
        self.resources = ResourceConfig(cfg.get("resources", {}))
        self.references = ReferenceDataConfig(
            cfg.get("references", {})
        )
        self.smaca_thresholds = SMAcaConfig(
            cfg.get("smaca_thresholds", {})
        )
        self.genebe_credentials = GeneBeConfig(
            cfg.get("genebe_credentials", {})
        )

class ResourceConfig:
    """
    Installed SFtool resource manifest configuration.
    """

    def __init__(self, cfg: dict[str, Any]):
        manifest = cfg.get("manifest")

        if not isinstance(manifest, str) or not manifest.strip():
            raise ValueError(
                "resources.manifest must be a non-empty path"
            )

        self.manifest: Path = (
            Path(manifest)
            .expanduser()
            .resolve()
        )

        self.raw_manifest: dict[str, Any] | None = None
        self.installed: dict[str, Any] | None = None

class GeneBeConfig:
    """
    GeneBe-related configuration.

    Structure preserved as-is.
    """
    def __init__(self, cfg: dict):
        self.api_key = cfg["api_key"]
        self.username = cfg["username"]

class SMAcaConfig:
    """
    SMAca thresholds and parameters.
    """
    def __init__(self, cfg: dict):
        self.cv_fail: float = cfg.get("cv_fail")
        self.cv_warn: float = cfg.get("cv_warn")
        self.low_cov_absolute: float = cfg.get("low_cov_absolute")
        self.low_cov_relative: float = cfg.get("low_cov_relative")


class ReferenceDataConfig:
    """
    User-provided reference genome overrides.

    Missing and null entries mean that SFtool should try the installed
    resource manifest.
    """

    def __init__(self, cfg: dict[str, Any]):
        _reject_obsolete_reference_fields(cfg)

        genomes = cfg.get("genomes", {})

        if genomes is None:
            genomes = {}

        if not isinstance(genomes, dict):
            raise ValueError(
                "references.genomes must be an object"
            )

        unexpected = set(genomes) - {
            "GRCh37",
            "GRCh38",
        }

        if unexpected:
            raise ValueError(
                "references.genomes contains unsupported "
                "assembly key(s): "
                + ", ".join(sorted(unexpected))
            )

        self.genomes: dict[str, Path | None] = {}

        for assembly in ("GRCh37", "GRCh38"):
            value = genomes.get(assembly)

            if value is None:
                self.genomes[assembly] = None
                continue

            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"references.genomes.{assembly} must be "
                    "a non-empty path or null"
                )

            self.genomes[assembly] = (
                Path(value)
                .expanduser()
                .resolve()
            )
class PathsConfig:
    """
    Runtime executable paths.
    """
    def __init__(self, cfg: dict):
        self.bcftools = cfg["bcftools"]
        self.java = cfg["java"]
        self.genebe = cfg["genebe"]
        self.bgzip = cfg["bgzip"]
        self.python = cfg["python"]
        self.pharmCAT = cfg["pharmCAT"]


def _reject_obsolete_configuration(
        cfg: dict[str, Any],
) -> None:
    obsolete = sorted(
        key for key in OBSOLETE_CONFIG_FIELDS
        if key in cfg
    )

    if not obsolete:
        return

    details = "; ".join(
        f"{key}: {OBSOLETE_CONFIG_FIELDS[key]}"
        for key in obsolete
    )

    raise ValueError(
        "Obsolete runtime configuration field(s): "
        f"{', '.join(obsolete)}. {details}"
    )

def _reject_obsolete_reference_fields(
        cfg: dict[str, Any],
) -> None:
    obsolete = sorted(
        key
        for key in OBSOLETE_REFERENCE_FIELDS
        if key in cfg
    )

    if not obsolete:
        return

    details = "; ".join(
        f"{field}: {OBSOLETE_REFERENCE_FIELDS[field]}"
        for field in obsolete
    )

    raise ValueError(
        "Obsolete references configuration field(s): "
        f"{', '.join(obsolete)}. {details}"
    )