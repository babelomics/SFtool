# modules/misc/errors.py

"""
Central definitions of custom exception types used across SFtool.
Each exception class represents a category of failure in the system.
"""


class BootstrapError(Exception):
    """
    Raised when execution bootstrap fails
    """
    pass


class ValidationError(BootstrapError):
    """
    Raised when Validation fails
    """
    pass


class RuntimeDependencyError(Exception):
    """
    Raised when the runtime environment is invalid:
    - missing third-party tools (bcftools, htslib, java, GeneBe, PharmCAT, etc.)
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

class ResourceManifestError(Exception):
    """
    Raised when installed resources cannot be loaded or resolved.
    """
    pass
