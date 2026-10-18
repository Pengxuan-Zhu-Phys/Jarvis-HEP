"""Download the public PDF manual without changing the local file on failure."""

from __future__ import annotations

from contextlib import suppress
import hashlib
from http.client import HTTPException
import math
import os
from pathlib import Path
import re
import tempfile
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

REPOSITORY_URL = "https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP"
MANUAL_FILENAME = "Jarvis-HEP-Manual.pdf"
MAX_MANUAL_BYTES = 256 * 1024 * 1024


class ManualDownloadError(Exception):
    """A user-facing failure with a stable exit status and recovery hint."""

    def __init__(self, message: str, *, hint: str, exit_code: int = 1) -> None:
        super().__init__(message)
        self.hint = hint
        self.exit_code = exit_code


def download_manual(
    *, output: str | Path = MANUAL_FILENAME, release: str = "latest",
    force: bool = False, timeout: float = 30.0,
) -> dict:
    """Fetch one named release asset, validate it, then atomically publish it."""
    if not math.isfinite(timeout) or timeout <= 0:
        raise ManualDownloadError(
            "Timeout must be a finite positive number of seconds.",
            hint="Use Jarvis download-manual --timeout 30.", exit_code=2,
        )
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}", release):
        raise ManualDownloadError(
            f"Invalid release: {release!r}.",
            hint="Use --release latest or an explicit tag such as --release v2.1.0.",
            exit_code=2,
        )
    if re.fullmatch(r"\d+\.\d+\.\d+", release):
        release = "v" + release
    source = (f"{REPOSITORY_URL}/releases/latest/download/{MANUAL_FILENAME}"
              if release == "latest" else
              f"{REPOSITORY_URL}/releases/download/{release}/{MANUAL_FILENAME}")
    try:
        target = Path(output).expanduser().absolute()
        if target.is_dir():
            target = target / MANUAL_FILENAME
        if target.is_symlink() or (target.exists() and not target.is_file()):
            raise ManualDownloadError(
                f"Destination is not a regular file: {target}",
                hint="Choose a regular PDF file path with --output.", exit_code=2,
            )
        if target.exists() and not force:
            raise ManualDownloadError(
                f"Destination already exists: {target}",
                hint="Choose another --output path or use --force to replace this file.",
                exit_code=2,
            )
        if not target.parent.is_dir():
            raise ManualDownloadError(
                f"Destination directory does not exist: {target.parent}",
                hint="Create the directory first or choose an existing directory with --output.",
                exit_code=2,
            )
    except (OSError, ValueError, RuntimeError) as exc:
        raise ManualDownloadError(
            f"Invalid destination: {exc}",
            hint="Choose an accessible file or existing directory with --output.",
            exit_code=2,
        ) from exc

    temporary: Path | None = None
    try:
        request = Request(source, headers={
            "User-Agent": "Jarvis-HEP-manual-download",
            "Accept": "application/pdf, application/octet-stream",
        })
        with urlopen(request, timeout=timeout) as response:
            raw_length = response.headers.get("Content-Length")
            try:
                expected = int(raw_length) if raw_length is not None else None
            except ValueError:
                expected = None
            if expected is not None and (expected < 0 or expected > MAX_MANUAL_BYTES):
                raise ValueError("The reported PDF size is invalid or exceeds 256 MiB.")
            with tempfile.NamedTemporaryFile(
                mode="wb", dir=target.parent, prefix=f".{target.name}.",
                suffix=".part", delete=False,
            ) as stream:
                temporary = Path(stream.name)
                size = 0
                header = b""
                tail = b""
                digest = hashlib.sha256()
                while chunk := response.read(64 * 1024):
                    size += len(chunk)
                    if size > MAX_MANUAL_BYTES:
                        raise ValueError("The PDF exceeds the 256 MiB download limit.")
                    header = (header + chunk)[:8]
                    tail = (tail + chunk)[-1024:]
                    if len(header) >= 5 and not header.startswith(b"%PDF-"):
                        raise ValueError("The downloaded content is not a PDF.")
                    stream.write(chunk)
                    digest.update(chunk)
                if expected is not None and size != expected:
                    raise ValueError("The PDF download is incomplete (size mismatch).")
                if not header.startswith(b"%PDF-") or b"%%EOF" not in tail:
                    raise ValueError("The PDF is empty, incomplete or invalid.")
                stream.flush()
                os.fsync(stream.fileno())
        if force:
            os.replace(temporary, target)
        else:
            # Atomic no-clobber publication, including another writer winning a race.
            os.link(temporary, target)
        return {
            "ok": True, "command": "download-manual", "release": release,
            "source": source, "path": str(target), "bytes": size,
            "sha256": digest.hexdigest(),
        }
    except HTTPError as exc:
        exc.close()
        if exc.code == 404:
            message = (f"No PDF manual is available for release {release!r}. "
                       f"The release or its {MANUAL_FILENAME} attachment may not be published yet.")
        else:
            message = f"Manual download failed: GitHub returned HTTP {exc.code}."
        raise ManualDownloadError(
            message, hint=f"Check {REPOSITORY_URL}/releases or retry later.",
        ) from exc
    except FileExistsError as exc:
        raise ManualDownloadError(
            f"Destination appeared during download: {target}",
            hint="Choose another --output path or use --force to replace this file.",
            exit_code=2,
        ) from exc
    except (OSError, URLError, HTTPException, ValueError) as exc:
        raise ManualDownloadError(
            f"Cannot download manual to {target}: {exc}",
            hint="Check the network and destination permissions, then retry Jarvis download-manual.",
        ) from exc
    finally:
        if temporary is not None:
            with suppress(OSError):
                temporary.unlink(missing_ok=True)
