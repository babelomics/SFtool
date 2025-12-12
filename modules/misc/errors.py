# modules/misc/errors.py

"""
Central definitions of custom exception types used across SFtool.
Each exception class represents a category of failure in the system.
"""


class ValidationError(Exception):
    """
    Raised when input JSON files (samples_info.json, config.json)
    are malformed, missing required fields, or violate biological rules.
    """
    pass


class RuntimeDependencyError(Exception):
    """
    Raised when the runtime environment is invalid:
    - missing binaries (bcftools, htslib, java, GeneBe, PharmCAT, etc.)
    - missing reference genomes or index files
    - missing catalog/category files
    - errors creating output directories
    """
    pass


class BioinformaticsProcessingError(Exception):
    """
    Raised when internal processing steps fail:
    - VCF normalization
    - BED generation
    - BED/VCF intersection
    - ClinVar parsing
    - subprocess command failures
    - corrupted temporary files
    """
    pass


class ModuleExecutionError(Exception):
    """
    Raised by individual modules (PR, RR, PGx) when:
    - required internal data is missing
    - a module-specific computation fails
    - GeneBe or PharmCAT return unexpected results
    - logic errors occur inside the module
    """
    pass
