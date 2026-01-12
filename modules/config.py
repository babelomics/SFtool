"""
Tool- and resource-level configuration wrappers.

Each class maps 1:1 to an existing top-level block
in config_example.json. Internal structures are preserved
exactly as provided.
"""


class CatalogConfig:
    """
    Catalog definitions (PR / RR / PGx).
    """
    def __init__(self, cfg: dict):
        self.cfg = cfg


class ClinVarConfig:
    """
    ClinVar database configuration.
    """
    def __init__(self, cfg: dict):
        self.cfg = cfg

    @property
    def db_path(self) -> str:
        return self.cfg["db_path"]

    @property
    def version(self) -> str:
        return self.cfg["version"]


class GeneBeConfig:
    """
    GeneBe-related configuration.

    Structure preserved as-is.
    """
    def __init__(self, cfg: dict):
        self.cfg = cfg

    @property
    def api_key(self) -> str:
        return self.cfg["api_key"]

    @property
    def username(self) -> str:
        return self.cfg["username"]


class SMAcaConfig:
    """
    SMAca thresholds and parameters.
    """
    def __init__(self, cfg: dict):
        self.cfg = cfg

        # Explicit attributes
        self.cv_fail: float = cfg.get("cv_fail")
        self.cv_warn: float = cfg.get("cv_warn")
        self.low_cov_absolute: float = cfg.get("low_cov_absolute")
        self.low_cov_relative: float = cfg.get("low_cov_relative")


class ReferenceDataConfig:
    """
    Reference-related resources, including:
      - reference genomes
      - gene-to-phenotype mappings
    """
    def __init__(self, cfg: dict):
        self.cfg = cfg

class PathsConfig:
    """
    Tool and resource paths.

    This maps to the top-level 'paths' block in config.json.
    """
    def __init__(self, cfg: dict):
        self.cfg = cfg

    @property
    def categories(self) -> str:
        """
        Base directory for PR / RR / PGx category BED/JSON files.
        """
        return self.cfg["categories"]

    @property
    def bcftools(self) -> str:
        """
        Path to bcftools executable.
        """
        return self.cfg["bcftools"]

    @property
    def java(self) -> str:
        return self.cfg["java"]

    @property
    def genebe(self) -> str:
        return self.cfg["genebe"]

    @property
    def htslib(self) -> str:
        return self.cfg["htslib"]

    @property
    def python(self) -> str:
        return self.cfg["python"]

    @property
    def pharmCAT(self) -> str:
        return self.cfg["pharmCAT"]