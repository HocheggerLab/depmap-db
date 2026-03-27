"""Exceptions for the download system."""


class DownloadError(Exception):
    """Base exception for download-related errors."""


class ValidationError(DownloadError):
    """Exception raised when file validation fails."""


class NetworkError(DownloadError):
    """Exception raised for network-related issues."""


class FileSystemError(DownloadError):
    """Exception raised for filesystem-related issues."""


class RateLimitError(DownloadError):
    """Exception raised when rate limits are exceeded."""
