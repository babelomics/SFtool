"""
Tool- and resource-level configuration wrappers.

Each class maps 1:1 to an existing top-level block
in config_example.json. Internal structures are preserved
exactly as provided.
"""

from typing import Dict

class Config:
    """
    Run-level configuration facade.

    Aggregates all config block wrappers.
    """

    def __init__(self, cfg: dict):
        self.version = cfg.get("version",{})
        self.paths = PathsConfig(cfg.get("paths", {}))
        self.references = ReferenceDataConfig(cfg.get("references", {}))
        self.catalogs = CatalogConfig(cfg.get("catalogs", {}))
        self.clinvar = ClinVarConfig(cfg.get("clinvar", {}))
        self.smaca_thresholds = SMAcaConfig(cfg.get("smaca_thresholds", {}))
        self.genebe_credentials = GeneBeConfig(cfg.get("genebe_credentials", {}))



class CatalogConfig:
    """
    Catalog definitions (PR / RR / PGx).
    """
    def __init__(self, cfg: dict):
        self.personal_risk_geneset = cfg["personal_risk_geneset"]
        self.reproductive_risk_geneset = cfg["reproductive_risk_geneset"]
        self.reproductive_risk_geneset_STR = cfg["reproductive_risk_geneset_STR"]


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
    Reference-related resources, including:
      - reference genomes
      - gene-to-phenotype mappings
      - pharmcat_positions_vcf
    """
    def __init__(self, cfg: dict):
        self.genomes: Dict = cfg.get("genomes", {})
        self.gene_to_phenotype_file: Optional[str] = cfg.get(
            "gene_to_phenotype_file"
        )
        self.pharmCAT_positions_vcf: Optional[str] = cfg.get(
            "pharmCAT_positions_vcf"
        )

class PathsConfig:
    """
    Tool and resource paths.

    This maps to the top-level 'paths' block in config.json.
    """
    def __init__(self, cfg: dict):
        self.categories = cfg["categories"]
        self.bcftools = cfg["bcftools"]
        self.java = cfg["java"]
        self.genebe = cfg["genebe"]
        self.bgzip = cfg["bgzip"]
        self.python = cfg["python"]
        self.pharmCAT = cfg["pharmCAT"]

