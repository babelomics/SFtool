"""
Runtime dependency checks for SFtool.

This module performs run-level preflight checks to ensure that
required external tools and resources exist before execution.
"""

import os
import re
import subprocess
from dataclasses import dataclass
from packaging.version import Version, InvalidVersion
from sftool.utils.errors import RuntimeDependencyError
from sftool.core.config import PathsConfig
from importlib.resources import files
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib



@dataclass
class ToolVersion:
    name: str
    path: str
    version: str
    minimum: str
    ok: bool



def load_runtime_requirements() -> dict[str, str]:
    requirements_file = files("sftool.data") / "runtime_requirements.toml"

    with requirements_file.open("rb") as fh:
        data = tomllib.load(fh)

    return {
        tool: values["minimum"]
        for tool, values in data.items()
    }


def run_command(cmd: list[str]) -> str:
    try:
        result = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        return (result.stdout or result.stderr).strip()
    except subprocess.CalledProcessError as e:
        raise RuntimeDependencyError(
            f"Command failed: {' '.join(map(str, cmd))}\n{e.stderr or e.stdout}"
        ) from e



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

def normalize_version(version: str) -> str:
    """
    Normalize non-standard version strings before comparison.

    Example:
    0.1.0-a.22 -> 0.1.0a22
    """
    return (
        version.strip()
        .replace("-a.", "a")
        .replace("-a", "a")
    )

def is_version_supported(version: str, minimum: str) -> bool:
    return Version(normalize_version(version)) >= Version(normalize_version(minimum))

def inspect_tool(name: str, path: str, version: str, minimum: str) -> ToolVersion:
    return ToolVersion(
        name=name,
        path=path,
        version=version,
        minimum=minimum,
        ok=is_version_supported(version, minimum),
    )

def get_python_version(python_path: str) -> str:
    output = run_command([python_path, "--version"])
    match = re.search(r"Python\s+([\d.]+)", output)
    if not match:
        raise RuntimeDependencyError(f"Could not parse Python version from: {output}")
    return match.group(1)


def get_java_version(java_path: str) -> str:
    output = run_command([java_path, "-version"])

    # openjdk version "21.0.2" 2024-01-16
    # java version "21.0.0" ...
    match = re.search(r'version\s+"([^"]+)"', output)
    if not match:
        raise RuntimeDependencyError(f"Could not parse Java version from: {output}")

    return match.group(1)


def get_bcftools_version(bcftools_path: str) -> str:
    output = run_command([bcftools_path, "--version"])

    # bcftools 1.21
    match = re.search(r"bcftools\s+([\w.\-]+)", output)
    if not match:
        raise RuntimeDependencyError(f"Could not parse bcftools version from: {output}")

    return match.group(1)


def get_bgzip_version(bgzip_path: str) -> str:
    output = run_command([bgzip_path, "--version"])

    # bgzip (htslib) 1.21
    # or bgzip 1.21
    match = re.search(r"(?:bgzip\s+\(htslib\)\s+|bgzip\s+)([\w.\-]+)", output)
    if not match:
        raise RuntimeDependencyError(f"Could not parse bgzip version from: {output}")

    return match.group(1)


def get_genebe_version(java_path: str, genebe_path: str) -> str:
    output = run_command([
        java_path,
        "-jar",
        genebe_path,
        "version",
    ])

    # current parsing in sample_tables.py expects:
    # version: ... ::
    match = re.search(r"version:\s*(.*?)\s*::", output)
    if not match:
        raise RuntimeDependencyError(f"Could not parse GeneBe version from: {output}")

    return match.group(1)


def get_pharmcat_version(java_path: str, pharmcat_path: str) -> str:
    output = run_command([
        java_path,
        "-jar",
        pharmcat_path,
        "-version",
    ])

    # PharmCAT often returns something like "PharmCAT 3.2.0" or just "3.2.0"
    match = re.search(r"(\d+\.\d+\.\d+)", output)
    if not match:
        raise RuntimeDependencyError(f"Could not parse PharmCAT version from: {output}")

    return match.group(1)


def get_runtime_versions(paths_cfg: PathsConfig) -> dict[str, ToolVersion]:
    requirements = load_runtime_requirements()

    return {
        "python": inspect_tool(
            "python",
            paths_cfg.python,
            get_python_version(paths_cfg.python),
            requirements["python"],
        ),
        "java": inspect_tool(
            "java",
            paths_cfg.java,
            get_java_version(paths_cfg.java),
            requirements["java"],
        ),
        "bcftools": inspect_tool(
            "bcftools",
            paths_cfg.bcftools,
            get_bcftools_version(paths_cfg.bcftools),
            requirements["bcftools"],
        ),
        "bgzip": inspect_tool(
            "bgzip",
            paths_cfg.bgzip,
            get_bgzip_version(paths_cfg.bgzip),
            requirements["bgzip"],
        ),
        "genebe": inspect_tool(
            "genebe",
            paths_cfg.genebe,
            get_genebe_version(paths_cfg.java, paths_cfg.genebe),
            requirements["genebe"],
        ),
        "pharmcat": inspect_tool(
            "pharmcat",
            paths_cfg.pharmCAT,
            get_pharmcat_version(paths_cfg.java, paths_cfg.pharmCAT),
            requirements["pharmcat"],
        ),
    }
# =====================================================
# Runtime dependency checking
# =====================================================

def check_runtime_dependencies(paths_cfg: PathsConfig):
    """
    Check that required runtime dependencies exist, are executable and their versions.

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

    # --------------------------------------------------------------
    # Version checks
    # --------------------------------------------------------------
    versions = get_runtime_versions(paths_cfg)

    failed = [
        tool
        for tool in versions.values()
        if not tool.ok
    ]

    if failed:
        details = "\n".join(
            f"- {tool.name}: found {tool.version}, required >= {tool.minimum}"
            for tool in failed
        )

        raise RuntimeDependencyError(
            "Unsupported runtime dependency versions:\n"
            f"{details}"
        )

    return versions