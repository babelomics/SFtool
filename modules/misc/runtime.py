

import os
from modules.misc.errors import RuntimeDependencyError


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

    paths = config_data.get("paths", {})

    # -----------------------------------------
    # GeneBeClient.jar
    # -----------------------------------------
    genebe_jar = paths.get("genebe")
    check_file_exists("GeneBeClient.jar", genebe_jar)

    # -----------------------------------------
    # bcftools executable
    # -----------------------------------------
    bcftools_bin = paths.get("bcftools")
    check_executable("bcftools binary", bcftools_bin)

    # -----------------------------------------
    # Java binary
    # -----------------------------------------
    java_bin = paths.get("java")
    check_executable("Java binary", java_bin)

    # -----------------------------------------
    # PharmCAT JAR file
    # -----------------------------------------
    pharmcat_jar = paths.get("pharmCAT")
    check_file_exists("PharmCAT JAR file", pharmcat_jar)

    # -----------------------------------------
    # Python interpreter
    # -----------------------------------------
    python_bin = paths.get("python")
    check_executable("Python interpreter", python_bin)

    # All checks passed
    return True
