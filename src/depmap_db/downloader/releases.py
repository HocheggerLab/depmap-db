"""Release tracking and refresh planning for DepMap datasets."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from ..config import get_settings
from ..utils.constants import DEPMAP_FILES, get_dataset_release_label
from .file_manager import FileManager


@dataclass
class ReleaseSnapshot:
    """A lightweight description of the configured release state."""

    release_label: str
    release_url: str
    datasets: dict[str, str]
    dataset_releases: dict[str, str]
    generated_at: str

    @classmethod
    def current(cls, datasets: list[str]) -> ReleaseSnapshot:
        """Create a snapshot from current configuration and dataset mapping."""
        settings = get_settings()
        dataset_releases = {
            dataset_name: get_dataset_release_label(
                dataset_name, settings.depmap.release_label
            )
            for dataset_name in datasets
        }
        distinct_release_labels = sorted(set(dataset_releases.values()))
        summary_release_label = (
            distinct_release_labels[0]
            if len(distinct_release_labels) == 1
            else "mixed"
        )
        return cls(
            release_label=summary_release_label,
            release_url=settings.depmap.release_url,
            datasets={
                dataset_name: DEPMAP_FILES[dataset_name].filename
                for dataset_name in datasets
            },
            dataset_releases=dataset_releases,
            generated_at=datetime.now(UTC).isoformat(),
        )


@dataclass
class RefreshPlan:
    """The actionable plan for a refresh run."""

    snapshot: ReleaseSnapshot
    datasets: list[str]
    datasets_to_download: list[str]
    cached_datasets: list[str]
    reason: str


class ReleaseTracker:
    """Persist the last applied release snapshot."""

    def __init__(self, tracking_file: Path | None = None) -> None:
        settings = get_settings()
        self.tracking_file = tracking_file or settings.depmap.release_tracking_file
        self.tracking_file.parent.mkdir(parents=True, exist_ok=True)

    def load(self) -> ReleaseSnapshot | None:
        """Load the last applied release snapshot if present."""
        if not self.tracking_file.exists():
            return None

        data = json.loads(self.tracking_file.read_text())
        if "dataset_releases" not in data:
            legacy_release_label = data.get("release_label") or "unknown"
            data["dataset_releases"] = dict.fromkeys(data.get("datasets", {}), legacy_release_label)
        return ReleaseSnapshot(**data)

    def save(self, snapshot: ReleaseSnapshot) -> None:
        """Persist a release snapshot."""
        self.tracking_file.write_text(json.dumps(asdict(snapshot), indent=2))


class RefreshPlanner:
    """Plan refresh runs against cached files and configured release state."""

    def __init__(
        self,
        file_manager: FileManager | None = None,
        tracker: ReleaseTracker | None = None,
    ) -> None:
        self.file_manager = file_manager or FileManager()
        self.tracker = tracker or ReleaseTracker()

    def build_plan(self, datasets: list[str]) -> RefreshPlan:
        """Build a refresh plan for the requested datasets."""
        snapshot = ReleaseSnapshot.current(datasets)
        previous = self.tracker.load()

        cached_datasets: list[str] = []
        datasets_to_download: list[str] = []

        release_changed = (
            previous is None
            or previous.dataset_releases != snapshot.dataset_releases
        )
        filenames_changed = previous is None or previous.datasets != snapshot.datasets

        for dataset_name in datasets:
            if release_changed or filenames_changed:
                datasets_to_download.append(dataset_name)
                continue

            if self.file_manager.is_file_cached(dataset_name, validate_integrity=False):
                cached_datasets.append(dataset_name)
            else:
                datasets_to_download.append(dataset_name)

        if previous is None:
            reason = "no previous release has been applied"
        elif release_changed:
            reason = "configured dataset release labels changed"
        elif filenames_changed:
            reason = "dataset-to-filename mapping changed"
        elif datasets_to_download:
            reason = "some requested datasets are missing from the local cache"
        else:
            reason = "all requested datasets are already cached for this release"

        return RefreshPlan(
            snapshot=snapshot,
            datasets=datasets,
            datasets_to_download=datasets_to_download,
            cached_datasets=cached_datasets,
            reason=reason,
        )

    def mark_applied(self, snapshot: ReleaseSnapshot) -> None:
        """Record that a release snapshot has been applied."""
        self.tracker.save(snapshot)


def refresh_plan_to_dict(plan: RefreshPlan) -> dict[str, Any]:
    """Convert a refresh plan to a serializable dictionary."""
    return {
        "release_label": plan.snapshot.release_label,
        "release_url": plan.snapshot.release_url,
        "datasets": plan.datasets,
        "dataset_releases": plan.snapshot.dataset_releases,
        "datasets_to_download": plan.datasets_to_download,
        "cached_datasets": plan.cached_datasets,
        "reason": plan.reason,
    }
