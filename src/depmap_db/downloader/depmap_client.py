"""DepMap API client for downloading data files."""

import asyncio
from pathlib import Path
from urllib.parse import urljoin

import httpx
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

from ..config import get_logger, get_settings
from ..utils.constants import DEPMAP_FILES, SUPPORTED_DATASETS
from ..utils.helpers import format_file_size, is_valid_depmap_file
from .exceptions import DownloadError, NetworkError, ValidationError

logger = get_logger(__name__)


class DepMapClient:
    """Client for downloading DepMap data files."""

    def __init__(self) -> None:
        """Initialize the DepMap client."""
        self.settings = get_settings()
        self.base_url = self.settings.depmap.base_url
        self.cache_dir = self.settings.depmap.cache_dir
        self.timeout = httpx.Timeout(self.settings.depmap.timeout_seconds)

        # Ensure cache directory exists
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    async def get_available_files(self) -> dict[str, dict[str, object]]:
        """Get list of available files from DepMap portal.

        Returns:
            Dictionary mapping dataset names to file information
        """
        # For now, we'll use our static mapping
        # In a future version, we could scrape the actual portal
        available: dict[str, dict[str, object]] = {}

        for dataset_name, file_info in DEPMAP_FILES.items():
            available[dataset_name] = {
                "filename": file_info.filename,
                "description": file_info.description,
                "table_name": file_info.table_name,
                "priority": file_info.priority,
                "url": urljoin(self.base_url, file_info.filename),
            }

        return available

    async def check_file_exists(
        self, filename: str
    ) -> tuple[bool, int | None]:
        """Check if a file exists on the DepMap portal.

        Args:
            filename: Name of the file to check

        Returns:
            Tuple of (exists, file_size_bytes)
        """
        url = urljoin(self.base_url, filename)

        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.head(url, follow_redirects=True)

                if response.status_code == 200:
                    content_length = response.headers.get("content-length")
                    file_size = int(content_length) if content_length else None
                    return True, file_size
                elif response.status_code == 404:
                    return False, None
                else:
                    logger.warning(
                        "Unexpected status code %s for %s",
                        response.status_code,
                        filename,
                    )
                    return False, None

        except httpx.RequestError as e:
            logger.error("Network error checking %s: %s", filename, e)
            raise NetworkError(f"Failed to check file existence: {e}") from e

    async def download_file(
        self,
        filename: str,
        destination: Path | None = None,
        force_download: bool = False,
        show_progress: bool = True,
    ) -> Path:
        """Download a single file from DepMap portal.

        Args:
            filename: Name of the file to download
            destination: Optional destination path. If None, saves to cache directory
            force_download: If True, download even if file already exists
            show_progress: Whether to show download progress

        Returns:
            Path to the downloaded file

        Raises:
            DownloadError: If download fails
            ValidationError: If file validation fails
        """
        if destination is None:
            destination = self.cache_dir / filename
        else:
            destination = Path(destination)

        # Check if file already exists and is valid
        if destination.exists() and not force_download:
            logger.info("File already exists: %s", destination)
            return destination

        # Validate filename
        if not is_valid_depmap_file(filename):
            raise ValidationError(f"Invalid DepMap filename: {filename}")

        url = urljoin(self.base_url, filename)

        # Check if file exists remotely
        file_exists, file_size = await self.check_file_exists(filename)
        if not file_exists:
            raise DownloadError(f"File not found on DepMap portal: {filename}")

        # Create destination directory
        destination.parent.mkdir(parents=True, exist_ok=True)

        # Download with progress tracking
        try:
            if show_progress and file_size:
                await self._download_with_progress(url, destination, file_size)
            else:
                await self._download_simple(url, destination)

        except (httpx.HTTPError, OSError) as e:
            # Clean up partial download
            if destination.exists():
                destination.unlink()
            raise DownloadError(f"Failed to download {filename}: {e}") from e

        # Validate downloaded file
        if not destination.exists():
            raise DownloadError(
                f"Download completed but file not found: {destination}"
            )

        actual_size = destination.stat().st_size
        if file_size and actual_size != file_size:
            destination.unlink()
            raise ValidationError(
                f"File size mismatch: expected {file_size}, got {actual_size}"
            )

        logger.info(
            "Successfully downloaded %s (%s)",
            filename,
            format_file_size(actual_size),
        )
        return destination

    async def _download_simple(self, url: str, destination: Path) -> None:
        """Download file without progress tracking."""
        async with httpx.AsyncClient(timeout=self.timeout) as client:
            async with client.stream(
                "GET", url, follow_redirects=True
            ) as response:
                response.raise_for_status()

                with destination.open("wb") as f:
                    async for chunk in response.aiter_bytes(chunk_size=8192):
                        f.write(chunk)

    async def _download_with_progress(
        self, url: str, destination: Path, file_size: int
    ) -> None:
        """Download file with progress bar."""
        with Progress(
            SpinnerColumn(),
            TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
            BarColumn(bar_width=None),
            "[progress.percentage]{task.percentage:>3.1f}%",
            "•",
            DownloadColumn(),
            "•",
            TimeElapsedColumn(),
        ) as progress:
            filename = destination.name
            task = progress.add_task(
                "download", filename=filename, total=file_size
            )

            async with httpx.AsyncClient(timeout=self.timeout) as client:
                async with client.stream(
                    "GET", url, follow_redirects=True
                ) as response:
                    response.raise_for_status()

                    with destination.open("wb") as f:
                        downloaded = 0
                        async for chunk in response.aiter_bytes(
                            chunk_size=8192
                        ):
                            f.write(chunk)
                            downloaded += len(chunk)
                            progress.update(task, completed=downloaded)

    async def download_datasets(
        self,
        datasets: list[str],
        destination_dir: Path | None = None,
        force_download: bool = False,
    ) -> dict[str, Path]:
        """Download multiple datasets.

        Args:
            datasets: List of dataset names to download
            destination_dir: Optional destination directory
            force_download: If True, download even if files exist

        Returns:
            Dictionary mapping dataset names to downloaded file paths
        """
        if destination_dir is None:
            destination_dir = self.cache_dir
        else:
            destination_dir = Path(destination_dir)

        # Validate dataset names
        invalid_datasets = set(datasets) - set(SUPPORTED_DATASETS)
        if invalid_datasets:
            raise ValidationError(
                f"Invalid datasets: {', '.join(invalid_datasets)}"
            )

        # Create download tasks
        download_tasks = {}
        for dataset in datasets:
            file_info = DEPMAP_FILES[dataset]
            filename = file_info.filename
            destination = destination_dir / filename

            download_tasks[dataset] = self.download_file(
                filename=filename,
                destination=destination,
                force_download=force_download,
                show_progress=True,
            )

        # Execute downloads concurrently
        logger.info("Starting download of %s datasets...", len(datasets))

        try:
            results = await asyncio.gather(
                *download_tasks.values(), return_exceptions=True
            )

            # Process results
            downloaded_files: dict[str, Path] = {}
            errors = []

            for dataset, result in zip(datasets, results, strict=False):
                if isinstance(result, BaseException):
                    errors.append(f"{dataset}: {result}")
                else:
                    downloaded_files[dataset] = result

            if errors:
                error_msg = "Some downloads failed:\n" + "\n".join(errors)
                logger.error(error_msg)
                raise DownloadError(error_msg)

            logger.info(
                "Successfully downloaded %s datasets",
                len(downloaded_files),
            )
            return downloaded_files

        except DownloadError:
            raise
        except (httpx.HTTPError, OSError) as e:
            logger.error("Download batch failed: %s", e)
            raise

    def get_cached_files(self) -> dict[str, Path]:
        """Get list of already cached files.

        Returns:
            Dictionary mapping dataset names to cached file paths
        """
        cached = {}

        for dataset_name, file_info in DEPMAP_FILES.items():
            cached_path = self.cache_dir / file_info.filename
            if cached_path.exists():
                cached[dataset_name] = cached_path

        return cached

    def clear_cache(self, datasets: list[str] | None = None) -> None:
        """Clear cached files.

        Args:
            datasets: Optional list of specific datasets to clear. If None, clears all.
        """
        if datasets is None:
            # Clear entire cache directory
            for file_path in self.cache_dir.glob("*.csv"):
                try:
                    file_path.unlink()
                    logger.info("Removed cached file: %s", file_path.name)
                except OSError as e:
                    logger.error("Failed to remove %s: %s", file_path, e)
        else:
            # Clear specific datasets
            for dataset in datasets:
                if dataset in DEPMAP_FILES:
                    file_info = DEPMAP_FILES[dataset]
                    cached_path = self.cache_dir / file_info.filename
                    if cached_path.exists():
                        try:
                            cached_path.unlink()
                            logger.info(
                                "Removed cached file: %s",
                                file_info.filename,
                            )
                        except OSError as e:
                            logger.error(
                                "Failed to remove %s: %s", cached_path, e
                            )


# Convenience functions for sync usage
def download_file_sync(filename: str, destination: Path | None = None) -> Path:
    """Synchronous wrapper for downloading a single file."""
    client = DepMapClient()
    return asyncio.run(client.download_file(filename, destination))


def download_datasets_sync(
    datasets: list[str], destination_dir: Path | None = None
) -> dict[str, Path]:
    """Synchronous wrapper for downloading multiple datasets."""
    client = DepMapClient()
    return asyncio.run(client.download_datasets(datasets, destination_dir))
