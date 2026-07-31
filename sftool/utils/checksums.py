# sftool/utils/checksums.py

from __future__ import annotations

import hashlib
from pathlib import Path


class ChecksumError(RuntimeError):
    """
    Raised when a file checksum cannot be calculated.
    """


def calculate_sha256(
        path: Path | str,
        *,
        chunk_size: int = 1024 * 1024,
) -> str:
    """
    Calculate the lowercase SHA-256 digest of a file.

    Args:
        path:
            File whose digest will be calculated.
        chunk_size:
            Number of bytes read per iteration.

    Returns:
        Lowercase hexadecimal SHA-256 digest.

    Raises:
        ValueError:
            If chunk_size is not greater than zero.
        ChecksumError:
            If the path does not exist, is not a regular file,
            or cannot be read.
    """
    file_path = Path(path)

    if chunk_size <= 0:
        raise ValueError(
            "chunk_size must be greater than zero"
        )

    if not file_path.exists():
        raise ChecksumError(
            f"Cannot calculate SHA-256 because the file "
            f"does not exist: {file_path}"
        )

    if not file_path.is_file():
        raise ChecksumError(
            f"Cannot calculate SHA-256 because the path "
            f"is not a file: {file_path}"
        )

    digest = hashlib.sha256()

    try:
        with file_path.open("rb") as handle:
            while chunk := handle.read(chunk_size):
                digest.update(chunk)
    except OSError as exc:
        raise ChecksumError(
            f"Could not read file while calculating SHA-256: "
            f"{file_path}"
        ) from exc

    return digest.hexdigest()