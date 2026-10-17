"""Installed Jarvis-HEP-Sampler release required by this Jarvis-HEP.

The plug-in CLI, registry, and scan identity all use this package. The historical
release floor is 1.0.3; this unreleased follow-up also requires host contract 2.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from packaging.version import InvalidVersion, Version

SDK_DISTRIBUTION = "Jarvis-HEP-Sampler"
SDK_SPECIFIER = ">=1.0.3"
_MIN_VERSION = Version("1.0.3")


class SamplerSDKError(RuntimeError):
    """The SDK is missing, too old, or lacks the required host interface."""


def _release(value: str) -> Version:
    try:
        return Version(value)
    except InvalidVersion as exc:
        raise SamplerSDKError(
            f"{SDK_DISTRIBUTION} version {value!r} is not a release number. "
            f'Install it with: python3 -m pip install "{SDK_DISTRIBUTION}{SDK_SPECIFIER}"'
        ) from exc


def require_sampler_sdk() -> str:
    """Return the installed distribution version when it meets the SDK minimum."""
    try:
        installed = version(SDK_DISTRIBUTION)
    except PackageNotFoundError as exc:
        raise SamplerSDKError(
            f"{SDK_DISTRIBUTION}{SDK_SPECIFIER} is not installed. "
            f'Install it with: python3 -m pip install "{SDK_DISTRIBUTION}{SDK_SPECIFIER}"'
        ) from exc
    if _release(installed) < _MIN_VERSION:
        raise SamplerSDKError(
            f"{SDK_DISTRIBUTION} {installed} does not satisfy {SDK_SPECIFIER}. "
            f'Install it with: python3 -m pip install "{SDK_DISTRIBUTION}{SDK_SPECIFIER}"'
        )
    import jarvis_hep_sampler

    if getattr(jarvis_hep_sampler, "HOST_CONTRACT_VERSION", 0) < 2:
        raise SamplerSDKError(
            f"{SDK_DISTRIBUTION} {installed} lacks sampler host contract 2. "
            "Install a matching SDK build from the updated Jarvis-HEP-Sampler source; "
            "the published 1.0.3 release predates this interface."
        )
    return installed


__all__ = [
    "SDK_DISTRIBUTION",
    "SDK_SPECIFIER",
    "SamplerSDKError",
    "require_sampler_sdk",
]
