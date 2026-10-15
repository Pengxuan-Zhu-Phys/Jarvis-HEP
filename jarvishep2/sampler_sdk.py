"""Installed Jarvis-HEP-Sampler release required by this Jarvis-HEP.

The plug-in CLI, registry, and scan identity all use this package. Any release
newer than 1.0.0 is accepted.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

SDK_DISTRIBUTION = "Jarvis-HEP-Sampler"
SDK_SPECIFIER = ">1"
_V1 = (1, 0, 0)


class SamplerSDKError(RuntimeError):
    """The installed sampler SDK is missing or not newer than 1.0.0."""


def _release(value: str) -> tuple[int, int, int]:
    numbers: list[int] = []
    for part in str(value).split("."):
        digits = ""
        for char in part:
            if not char.isdigit():
                break
            digits += char
        if not digits:
            break
        numbers.append(int(digits))
    if not numbers:
        raise SamplerSDKError(
            f"{SDK_DISTRIBUTION} version {value!r} is not a release number. "
            f'Install it with: python3 -m pip install "{SDK_DISTRIBUTION}{SDK_SPECIFIER}"'
        )
    return tuple((numbers + [0, 0])[:3])


def require_sampler_sdk() -> str:
    """Return the installed distribution version when it is newer than 1.0.0."""
    try:
        installed = version(SDK_DISTRIBUTION)
    except PackageNotFoundError as exc:
        raise SamplerSDKError(
            f"{SDK_DISTRIBUTION}{SDK_SPECIFIER} is not installed. "
            f'Install it with: python3 -m pip install "{SDK_DISTRIBUTION}{SDK_SPECIFIER}"'
        ) from exc
    if _release(installed) <= _V1:
        raise SamplerSDKError(
            f"{SDK_DISTRIBUTION} {installed} does not satisfy {SDK_SPECIFIER}. "
            f'Install it with: python3 -m pip install "{SDK_DISTRIBUTION}{SDK_SPECIFIER}"'
        )
    return installed


__all__ = [
    "SDK_DISTRIBUTION",
    "SDK_SPECIFIER",
    "SamplerSDKError",
    "require_sampler_sdk",
]
