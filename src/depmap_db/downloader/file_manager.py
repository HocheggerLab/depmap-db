"""File management for DepMap data downloads."""

import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from ..config import get_logger, get_settings
from ..utils.helpers import format_file_size
from .exceptions import FileSystemError

logger = get_logger(__name__)


@dataclass
class DownloadResult:
    """Result of a file download operation."""

    dataset_name: str
    filename: str
    file_path: Path
    file_size: int
    checksum: str
    download_time: datetime
    source_url: str

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for serialization."""
        result = asdict(self)
        result["file_path"] = str(result["file_path"])
        result["download_time"] = result["download_time"].isoformat()
        return result

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "DownloadResult":
        """Create from dictionary."""
        data["file_path"] = Path(data["file_path"])
        data["download_time"] = datetime.fromisoformat(data["download_time"])
        return cls(**data)


class FileManager:
    """Manages downloaded DepMap files and metadata."""

    def __init__(self, cache_dir: Path | None = None) -> None:
        """Initialize file manager.

        Args:
            cache_dir: Optional cache directory. If None, uses settings default.
        """
        self.settings = get_settings()
        self.cache_dir = cache_dir or self.settings.depmap.cache_dir
        self.metadata_file = self.cache_dir / "download_metadata.json"

        # Ensure cache directory exists
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Load existing metadata
        self._metadata: dict[str, DownloadResult] = {}
        self._load_metadata()

    def _load_metadata(self) -> None:
        """Load download metadata from file."""
        if self.metadata_file.exists():
            try:
                with self.metadata_file.open("r") as f:
                    data = json.load(f)

                for dataset_name, item_data in data.items():
                    try:
                        self._metadata[dataset_name] = (
                            DownloadResult.from_dict(item_data)
                        )
                    except (KeyError, TypeError, ValueError) as e:
                        logger.warning(
                            "Invalid metadata for %s: %s", dataset_name, e
                        )

                logger.info(
                    "Loaded metadata for %s files", len(self._metadata)
                )

            except (json.JSONDecodeError, OSError) as e:
                logger.error("Failed to load metadata file: %s", e)
                self._metadata = {}

    def _save_metadata(self) -> None:
        """Save download metadata to file."""
        try:
            metadata_dict = {
                name: result.to_dict()
                for name, result in self._metadata.items()
            }

            with self.metadata_file.open("w") as f:
                json.dump(metadata_dict, f, indent=2)

            logger.debug("Saved download metadata")

        except OSError as e:
            logger.error("Failed to save metadata: %s", e)
            raise FileSystemError(f"Failed to save metadata: {e}") from e

    def calculate_file_checksum(
        self, file_path: Path, algorithm: str = "md5"
    ) -> str:
        """Calculate checksum for a file.

        Args:
            file_path: Path to the file
            algorithm: Hash algorithm to use

        Returns:
            Hexadecimal checksum string
        """
        hash_func = getattr(hashlib, algorithm.lower())()

        with file_path.open("rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_func.update(chunk)

        result: str = hash_func.hexdigest()
        return result

    def register_download(
        self,
        dataset_name: str,
        file_path: Path,
        source_url: str,
        force_checksum: bool = True,
    ) -> DownloadResult:
        """Register a downloaded file with metadata.

        Args:
            dataset_name: Name of the dataset
            file_path: Path to the downloaded file
            source_url: URL the file was downloaded from
            force_checksum: Whether to force checksum calculation

        Returns:
            DownloadResult with file metadata
        """
        if not file_path.exists():
            raise FileSystemError(f"File does not exist: {file_path}")

        # Calculate file metadata
        file_size = file_path.stat().st_size
        checksum = (
            self.calculate_file_checksum(file_path) if force_checksum else ""
        )

        result = DownloadResult(
            dataset_name=dataset_name,
            filename=file_path.name,
            file_path=file_path,
            file_size=file_size,
            checksum=checksum,
            download_time=datetime.now(),
            source_url=source_url,
        )

        # Store metadata
        self._metadata[dataset_name] = result
        self._save_metadata()

        logger.info(
            "Registered download: %s (%s)",
            dataset_name,
            format_file_size(file_size),
        )
        return result

    def get_file_info(self, dataset_name: str) -> DownloadResult | None:
        """Get download information for a dataset.

        Args:
            dataset_name: Name of the dataset

        Returns:
            DownloadResult if found, None otherwise
        """
        return self._metadata.get(dataset_name)

    def is_file_cached(
        self, dataset_name: str, validate_integrity: bool = False
    ) -> bool:
        """Check if a file is cached and optionally validate integrity.

        Args:
            dataset_name: Name of the dataset
            validate_integrity: Whether to validate file checksum

        Returns:
            True if file is cached and valid, False otherwise
        """
        file_info = self.get_file_info(dataset_name)
        if not file_info:
            return False

        # Check if file exists
        if not file_info.file_path.exists():
            logger.warning("Cached file missing: %s", file_info.file_path)
            # Remove stale metadata
            del self._metadata[dataset_name]
            self._save_metadata()
            return False

        # Check file size
        current_size = file_info.file_path.stat().st_size
        if current_size != file_info.file_size:
            logger.warning(
                "File size mismatch for %s: expected %s, got %s",
                dataset_name,
                file_info.file_size,
                current_size,
            )
            return False

        # Validate checksum if requested and available
        if validate_integrity and file_info.checksum:
            try:
                current_checksum = self.calculate_file_checksum(
                    file_info.file_path
                )
                if current_checksum != file_info.checksum:
                    logger.warning("Checksum mismatch for %s", dataset_name)
                    return False
            except OSError as e:
                logger.error(
                    "Failed to validate checksum for %s: %s", dataset_name, e
                )
                return False

        return True

    def get_cached_files(
        self, validate_integrity: bool = False
    ) -> dict[str, DownloadResult]:
        """Get all cached files.

        Args:
            validate_integrity: Whether to validate file integrity

        Returns:
            Dictionary mapping dataset names to DownloadResult
        """
        cached = {}

        for dataset_name, file_info in self._metadata.items():
            if self.is_file_cached(dataset_name, validate_integrity):
                cached[dataset_name] = file_info

        return cached

    def get_cache_stats(self) -> dict[str, Any]:
        """Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        cached_files = self.get_cached_files()
        total_size = sum(info.file_size for info in cached_files.values())

        return {
            "cached_datasets": len(cached_files),
            "total_size_bytes": total_size,
            "total_size_formatted": format_file_size(total_size),
            "cache_directory": str(self.cache_dir),
            "datasets": list(cached_files.keys()),
        }

    def remove_cached_file(self, dataset_name: str) -> bool:
        """Remove a cached file and its metadata.

        Args:
            dataset_name: Name of the dataset to remove

        Returns:
            True if removed successfully, False if not found
        """
        file_info = self.get_file_info(dataset_name)
        if not file_info:
            return False

        try:
            # Remove file if it exists
            if file_info.file_path.exists():
                file_info.file_path.unlink()
                logger.info("Removed cached file: %s", file_info.file_path)

            # Remove metadata
            del self._metadata[dataset_name]
            self._save_metadata()

            return True

        except OSError as e:
            logger.error(
                "Failed to remove cached file %s: %s", dataset_name, e
            )
            return False

    def clear_cache(self, datasets: list[str] | None = None) -> int:
        """Clear cached files.

        Args:
            datasets: Optional list of specific datasets to clear. If None, clears all.

        Returns:
            Number of files removed
        """
        if datasets is None:
            datasets = list(self._metadata.keys())

        removed_count = 0
        for dataset_name in datasets:
            if self.remove_cached_file(dataset_name):
                removed_count += 1

        logger.info("Cleared %s cached files", removed_count)
        return removed_count

    def verify_all_files(self) -> dict[str, bool]:
        """Verify integrity of all cached files.

        Returns:
            Dictionary mapping dataset names to validation results
        """
        results: dict[str, bool] = {}

        for dataset_name in self._metadata:
            try:
                results[dataset_name] = self.is_file_cached(
                    dataset_name, validate_integrity=True
                )
            except (OSError, ValueError) as e:
                logger.error("Failed to verify %s: %s", dataset_name, e)
                results[dataset_name] = False

        return results

    def get_file_path(self, dataset_name: str) -> Path | None:
        """Get the file path for a cached dataset.

        Args:
            dataset_name: Name of the dataset

        Returns:
            Path to the file if cached, None otherwise
        """
        file_info = self.get_file_info(dataset_name)
        if file_info and self.is_file_cached(dataset_name):
            return file_info.file_path
        return None

    def cleanup_orphaned_files(self) -> int:
        """Remove files in cache directory that aren't tracked in metadata.

        Returns:
            Number of orphaned files removed
        """
        tracked_files = {
            info.file_path.name for info in self._metadata.values()
        }

        removed_count = 0

        # Find CSV files that aren't tracked
        for file_path in self.cache_dir.glob("*.csv"):
            if file_path.name not in tracked_files:
                try:
                    file_path.unlink()
                    logger.info("Removed orphaned file: %s", file_path.name)
                    removed_count += 1
                except OSError as e:
                    logger.error(
                        "Failed to remove orphaned file %s: %s",
                        file_path,
                        e,
                    )

        return removed_count
