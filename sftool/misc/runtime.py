"""
Runtime dependency checks for SFtool.

This module performs run-level preflight checks to ensure that
required external tools and resources exist before execution.
"""

import os
from sftool.misc.errors import RuntimeDependencyError
from sftool.core.config import PathsConfig


# =====================================================
# Utility helpers
# =====================================================
def check_file_exists(label, path):
    """
    Checks that a file exists. Raises RuntimeDependencyError if missing.
    Appropriate for JAR files or any resource that does NOT need
    executable permissions.
    """
    if not isinstance(path, str) or path.strip() == "":
        raise RuntimeDependencyError(f"Invalid path for {label}: empty or not a string")

    if not os.path.exists(path):
        raise RuntimeDependencyError(f"{label} not found at: {path}")


def check_executable(label, path):
    """
    Checks that a binary exists and IS executable.
    Appropriate for bcftools, java, python, etc.
    """
    check_file_exists(label, path)

    if not os.access(path, os.X_OK):
        raise RuntimeDependencyError(
            f"{label} exists but is NOT executable: {path}"
        )


# =====================================================
# Runtime dependency checking
# =====================================================
def check_runtime_dependencies(config_data):
    """
    Checks system-level runtime dependencies required to RUN SFtool.
    Checks included here (adapted from the original SFtool):
      - GeneBeClient.jar exists
      - bcftools binary exists AND is executable
      - Java binary exists AND is executable
      - pharmcat.jar exists
      - Python interpreter exists AND is executable

    Any missing or invalid dependency raises RuntimeDependencyError.
    SFtool.py is responsible for catching and reporting errors.
    """

def check_runtime_dependencies(paths_cfg: PathsConfig):
    """
    Check that required runtime dependencies exist.

    This function preserves the semantics of the legacy runtime.py
    and only adapts the input interface to PathsConfig.
    """

    # --------------------------------------------------------------
    # Executables (must be executable)
    # --------------------------------------------------------------
    check_executable("python interpreter", paths_cfg.python)
    check_executable("bcftools binary", paths_cfg.bcftools)
    check_executable("java binary", paths_cfg.java)
    check_executable("bgzip binary", paths_cfg.bgzip)

    # --------------------------------------------------------------
    # Files / resources (existence only)
    # --------------------------------------------------------------
    check_file_exists("GeneBe JAR file", paths_cfg.genebe)
    check_file_exists("PharmCAT JAR file", paths_cfg.pharmCAT)

