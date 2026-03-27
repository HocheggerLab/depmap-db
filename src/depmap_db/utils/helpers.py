"""Helper utility functions for common operations."""

import hashlib
from pathlib import Path


def validate_file_checksum(
    file_path: Path, expected_checksum: str, algorithm: str = "md5"
) -> bool:
    """Validate a file's checksum.

    Args:
        file_path: Path to the file to validate
        expected_checksum: Expected checksum value
        algorithm: Hashing algorithm to use (md5, sha256, etc.)

    Returns:
        True if checksum matches, False otherwise
    """
    if not file_path.exists():
        return False

    # Get the appropriate hash function
    try:
        hash_func = getattr(hashlib, algorithm.lower())()
    except AttributeError as err:
        raise ValueError(f"Unsupported hash algorithm: {algorithm}") from err

    # Calculate file checksum
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_func.update(chunk)

    calculated_checksum: str = hash_func.hexdigest()
    return calculated_checksum.lower() == expected_checksum.lower()


def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format.

    Args:
        size_bytes: Size in bytes

    Returns:
        Formatted size string (e.g., "1.5 GB", "256 MB")
    """
    if size_bytes == 0:
        return "0 B"

    units = ["B", "KB", "MB", "GB", "TB"]
    unit_index = 0
    size = float(size_bytes)

    while size >= 1024.0 and unit_index < len(units) - 1:
        size /= 1024.0
        unit_index += 1

    # Format with appropriate precision
    if unit_index == 0:
        return f"{int(size)} {units[unit_index]}"
    else:
        return f"{size:.1f} {units[unit_index]}"


def sanitize_filename(filename: str) -> str:
    """Sanitize a filename by removing or replacing invalid characters.

    Args:
        filename: Original filename

    Returns:
        Sanitized filename safe for filesystem use
    """
    # Characters to remove or replace
    invalid_chars = '<>:"/\\|?*'

    # Replace invalid characters with underscores
    sanitized = filename
    for char in invalid_chars:
        sanitized = sanitized.replace(char, "_")

    # Remove leading/trailing whitespace and dots
    sanitized = sanitized.strip(". ")

    # Ensure filename is not empty
    if not sanitized:
        sanitized = "unnamed_file"

    return sanitized


def get_file_extension(filename: str) -> str:
    """Get the file extension from a filename.

    Args:
        filename: The filename

    Returns:
        File extension including the dot (e.g., ".csv")
    """
    path = Path(filename)
    return path.suffix.lower()


def is_valid_depmap_file(filename: str) -> bool:
    """Check if a filename appears to be a valid DepMap data file.

    Args:
        filename: The filename to check

    Returns:
        True if filename appears valid, False otherwise
    """
    from .constants import SUPPORTED_FILE_EXTENSIONS

    # Check file extension
    extension = get_file_extension(filename)
    if extension not in SUPPORTED_FILE_EXTENSIONS:
        return False

    # Check for common DepMap naming patterns
    depmap_patterns = [
        "CRISPR",
        "Model",
        "Gene",
        "Omics",
        "Achilles",
        "Expression",
        "Mutation",
        "Copy",
        "Screen",
    ]

    filename_upper = filename.upper()
    return any(
        pattern.upper() in filename_upper for pattern in depmap_patterns
    )


def extract_release_version(filename: str) -> str | None:
    """Extract DepMap release version from filename if present.

    Args:
        filename: The filename to analyze

    Returns:
        Release version string if found, None otherwise
    """
    # Common DepMap release version patterns: 25Q2, 24Q4, etc.
    import re

    patterns = [
        r"(\d{2}Q[1-4])",  # 25Q2, 24Q4, etc.
        r"(\d{4}_Q[1-4])",  # 2025_Q2, etc.
        r"v(\d+\.\d+)",  # v1.0, v2.1, etc.
    ]

    for pattern in patterns:
        match = re.search(pattern, filename, re.IGNORECASE)
        if match:
            return match.group(1)

    return None


def create_backup_filename(file_path: Path) -> Path:
    """Create a backup filename for an existing file.

    Args:
        file_path: Original file path

    Returns:
        Path for backup file
    """
    counter = 1
    backup_path = file_path.with_suffix(f".bak{file_path.suffix}")

    # Find an available backup filename
    while backup_path.exists():
        counter += 1
        backup_path = file_path.with_suffix(f".bak{counter}{file_path.suffix}")

    return backup_path
