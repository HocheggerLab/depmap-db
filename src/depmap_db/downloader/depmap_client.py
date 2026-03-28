"""DepMap API client for downloading data files."""

import asyncio
import csv
import hashlib
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Any

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
from ..utils.constants import (
    DEPMAP_FILES,
    SUPPORTED_DATASETS,
    get_dataset_release_label,
)
from ..utils.helpers import format_file_size, is_valid_depmap_file
from .exceptions import DownloadError, NetworkError, ValidationError

logger = get_logger(__name__)


@dataclass(frozen=True)
class ManifestFile:
    """A single file entry from the DepMap download manifest."""

    release: str
    release_date: str
    filename: str
    url: str
    md5_hash: str


class DepMapClient:
    """Client for downloading DepMap data files."""

    def __init__(self) -> None:
        """Initialize the DepMap client."""
        self.settings = get_settings()
        self.manifest_url = self.settings.depmap.base_url
        self.datasets_url = self.settings.depmap.datasets_url
        self.cache_dir = self.settings.depmap.cache_dir
        self.timeout = httpx.Timeout(self.settings.depmap.timeout_seconds)
        self.release_label = self.settings.depmap.release_label
        self.verify_checksums = self.settings.depmap.verify_checksums
        self._manifest_cache: list[ManifestFile] | None = None

        # Ensure cache directory exists
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    async def _fetch_manifest(self) -> list[ManifestFile]:
        """Fetch and parse the DepMap download manifest CSV."""
        if self._manifest_cache is not None:
            return self._manifest_cache

        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(self.manifest_url)
                response.raise_for_status()
        except httpx.RequestError as e:
            logger.error("Network error fetching DepMap manifest: %s", e)
            raise NetworkError(f"Failed to fetch DepMap manifest: {e}") from e
        except httpx.HTTPStatusError as e:
            logger.error("HTTP error fetching DepMap manifest: %s", e)
            raise DownloadError(f"Failed to fetch DepMap manifest: {e}") from e

        manifest = [
            ManifestFile(
                release=row["release"],
                release_date=row["release_date"],
                filename=row["filename"],
                url=row["url"],
                md5_hash=row.get("md5_hash", ""),
            )
            for row in csv.DictReader(io.StringIO(response.text))
            if row.get("release")
            and row.get("release_date")
            and row.get("filename")
            and row.get("url")
        ]
        self._manifest_cache = manifest
        return manifest

    async def get_dataset_catalogue(self) -> list[dict[str, Any]]:
        """Fetch the DepMap dataset catalogue JSON."""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.get(self.datasets_url)
                response.raise_for_status()
        except httpx.RequestError as e:
            logger.error("Network error fetching DepMap dataset catalogue: %s", e)
            raise NetworkError(
                f"Failed to fetch DepMap dataset catalogue: {e}"
            ) from e
        except httpx.HTTPStatusError as e:
            logger.error("HTTP error fetching DepMap dataset catalogue: %s", e)
            raise DownloadError(
                f"Failed to fetch DepMap dataset catalogue: {e}"
            ) from e

        data = response.json()
        if not isinstance(data, list):
            raise ValidationError("Unexpected dataset catalogue response shape")
        return [item for item in data if isinstance(item, dict)]

    async def resolve_manifest_file(
        self, filename: str, release_label: str | None = None
    ) -> ManifestFile:
        """Resolve a filename to its manifest entry and signed download URL."""
        manifest = await self._fetch_manifest()
        matches = [entry for entry in manifest if entry.filename == filename]
        if not matches:
            raise DownloadError(f"File not found in DepMap manifest: {filename}")

        target_release = release_label or self.release_label
        release_matches = [
            entry for entry in matches if entry.release == target_release
        ]
        if release_matches:
            return release_matches[0]

        # Fallback to the newest manifest entry for this filename.
        logger.warning(
            "Configured release '%s' not found for %s; using newest manifest entry",
            target_release,
            filename,
        )
        return max(matches, key=lambda entry: (entry.release_date, entry.release))

    async def get_download_sources(
        self, datasets: list[str]
    ) -> dict[str, ManifestFile]:
        """Resolve dataset names to manifest-backed download entries."""
        invalid_datasets = set(datasets) - set(SUPPORTED_DATASETS)
        if invalid_datasets:
            raise ValidationError(
                f"Invalid datasets: {', '.join(sorted(invalid_datasets))}"
            )

        sources: dict[str, ManifestFile] = {}
        for dataset_name in datasets:
            filename = DEPMAP_FILES[dataset_name].filename
            release_label = get_dataset_release_label(
                dataset_name, self.release_label
            )
            sources[dataset_name] = await self.resolve_manifest_file(
                filename, release_label=release_label
            )
        return sources

    async def get_available_files(self) -> dict[str, dict[str, object]]:
        """Get available configured datasets for the selected release."""
        sources = await self.get_download_sources(list(DEPMAP_FILES.keys()))
        available: dict[str, dict[str, object]] = {}

        for dataset_name, file_info in DEPMAP_FILES.items():
            source = sources[dataset_name]
            available[dataset_name] = {
                "filename": file_info.filename,
                "description": file_info.description,
                "table_name": file_info.table_name,
                "priority": file_info.priority,
                "url": source.url,
                "md5_hash": source.md5_hash,
                "release": source.release,
                "release_date": source.release_date,
            }

        return available

    async def _get_remote_file_size(self, url: str) -> int | None:
        """Try to determine remote file size from the signed download URL."""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.head(url, follow_redirects=True)
        except httpx.RequestError:
            return None

        if response.status_code != 200:
            return None

        content_length = response.headers.get("content-length")
        return int(content_length) if content_length else None

    async def check_file_exists(
        self, filename: str
    ) -> tuple[bool, int | None]:
        """Check if a file is present in the DepMap manifest for the release."""
        try:
            manifest_file = await self.resolve_manifest_file(filename)
        except DownloadError:
            return False, None

        file_size = await self._get_remote_file_size(manifest_file.url)
        return True, file_size

    def _verify_md5_checksum(self, file_path: Path, expected_md5: str) -> None:
        """Validate a downloaded file against the manifest md5 checksum."""
        digest = hashlib.md5(usedforsecurity=False)
        with file_path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(8192), b""):
                digest.update(chunk)

        actual_md5 = digest.hexdigest()
        if actual_md5 != expected_md5:
            raise ValidationError(
                f"MD5 mismatch for {file_path.name}: expected {expected_md5}, got {actual_md5}"
            )

    async def download_file(
        self,
        filename: str,
        destination: Path | None = None,
        force_download: bool = False,
        show_progress: bool = True,
    ) -> Path:
        """Download a single file from the DepMap portal."""
        if destination is None:
            destination = self.cache_dir / filename
        else:
            destination = Path(destination)

        if destination.exists() and not force_download:
            logger.info("File already exists: %s", destination)
            return destination

        if not is_valid_depmap_file(filename):
            raise ValidationError(f"Invalid DepMap filename: {filename}")

        manifest_file = await self.resolve_manifest_file(filename)
        file_size = await self._get_remote_file_size(manifest_file.url)

        destination.parent.mkdir(parents=True, exist_ok=True)

        try:
            if show_progress and file_size:
                await self._download_with_progress(
                    manifest_file.url, destination, file_size
                )
            else:
                await self._download_simple(manifest_file.url, destination)
        except (httpx.HTTPError, OSError) as e:
            if destination.exists():
                destination.unlink()
            raise DownloadError(f"Failed to download {filename}: {e}") from e

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

        if self.verify_checksums and manifest_file.md5_hash:
            try:
                self._verify_md5_checksum(destination, manifest_file.md5_hash)
            except ValidationError:
                destination.unlink()
                raise

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
        """Download multiple datasets."""
        if destination_dir is None:
            destination_dir = self.cache_dir
        else:
            destination_dir = Path(destination_dir)

        invalid_datasets = set(datasets) - set(SUPPORTED_DATASETS)
        if invalid_datasets:
            raise ValidationError(
                f"Invalid datasets: {', '.join(sorted(invalid_datasets))}"
            )

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

        logger.info("Starting download of %s datasets...", len(datasets))

        try:
            results = await asyncio.gather(
                *download_tasks.values(), return_exceptions=True
            )

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
        """Get list of already cached files."""
        cached = {}

        for dataset_name, file_info in DEPMAP_FILES.items():
            cached_path = self.cache_dir / file_info.filename
            if cached_path.exists():
                cached[dataset_name] = cached_path

        return cached

    def clear_cache(self, datasets: list[str] | None = None) -> None:
        """Clear cached files."""
        if datasets is None:
            for file_path in self.cache_dir.glob("*.csv"):
                try:
                    file_path.unlink()
                    logger.info("Removed cached file: %s", file_path.name)
                except OSError as e:
                    logger.error("Failed to remove %s: %s", file_path, e)
        else:
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


def get_download_sources_sync(datasets: list[str]) -> dict[str, ManifestFile]:
    """Synchronous wrapper for manifest-backed dataset resolution."""
    client = DepMapClient()
    return asyncio.run(client.get_download_sources(datasets))


def get_dataset_catalogue_sync() -> list[dict[str, Any]]:
    """Synchronous wrapper for the dataset catalogue JSON endpoint."""
    client = DepMapClient()
    return asyncio.run(client.get_dataset_catalogue())
