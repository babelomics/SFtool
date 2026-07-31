"""
Tool- and resource-level configuration wrappers.

Each class maps 1:1 to an existing top-level block
in config_example.json. Internal structures are preserved
exactly as provided.
"""
from __future__ import annotations

from typing import Dict, Any
from importlib.resources import files
from pathlib import Path


class Config:
    """
    Run-level configuration facade.
    """

    def __init__(self, cfg: dict[str, Any]):
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

class CatalogConfig:
    """
    Catalog definitions (PR / RR / PGx).
    """
    def __init__(self, cfg: dict):
        categories_dir = files("sftool.data.categories")


        self.personal_risk_geneset = cfg.get(
            "personal_risk_geneset",
            categories_dir / "PR" / "PR_risk_genes_ACMG_SF_v3.1.csv"
        )

        self.reproductive_risk_geneset = cfg.get(
            "reproductive_risk_geneset",
            categories_dir / "RR" / "RR_risk_genes_ACMG_CS_v2021.csv"
        )

        self.reproductive_risk_geneset_STR = cfg.get(
            "reproductive_risk_geneset_STR",
            categories_dir / "RR" / "RR_risk_genes_STR_ACMG_CS_v2021.csv"
        )


class ClinVarConfig:
    """
    ClinVar database configuration.
    """
    def __init__(self, cfg: dict):
        self.db_path = cfg["db_path"]
        self.version = cfg["version"]


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
    Tool and resource paths.

    This maps to the top-level 'paths' block in config.json.
    """
    def __init__(self, cfg: dict):
        self.bcftools = cfg["bcftools"]
        self.java = cfg["java"]
        self.genebe = cfg["genebe"]
        self.bgzip = cfg["bgzip"]
        self.python = cfg["python"]
        self.pharmCAT = cfg["pharmCAT"]

