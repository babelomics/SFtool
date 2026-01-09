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


class SMAcaConfig:
    """
    SMAca thresholds and parameters.
    """
    def __init__(self, cfg: dict):
        self.cfg = cfg


class ReferenceDataConfig:
    """
    Reference-related resources, including:
      - reference genomes
      - gene-to-phenotype mappings
    """
    def __init__(self, cfg: dict):
        self.cfg = cfg
