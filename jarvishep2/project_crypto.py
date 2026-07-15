#!/usr/bin/env python3
"""OpenSSL-compatible AES-256-CBC (PBKDF2) helpers for restricted project packs.

Scheme id: ``openssl-aes-256-cbc`` — same on-disk format as::

    openssl enc -aes-256-cbc -pbkdf2 -salt -in plain.tar.gz -out pack.jenc -pass pass:KEY

Uses the system ``openssl`` binary so no extra PyPI crypto package is required.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from typing import BinaryIO

# OpenSSL salted ciphertext magic (Salted__)
OPENSSL_SALTED_MAGIC = b"Salted__"
ENCRYPTION_SCHEME_OPENSSL = "openssl-aes-256-cbc"
ENCRYPTION_SCHEME_NONE = "none"
PROJECT_FETCH_KEY_ENV = "JARVIS_PROJECT_FETCH_KEY"


class ProjectCryptoError(RuntimeError):
    pass


def resolve_fetch_key(explicit: str | None = None) -> str | None:
    """Return explicit key, else ``JARVIS_PROJECT_FETCH_KEY``, else None."""
    if explicit is not None and str(explicit).strip():
        return str(explicit).strip()
    env = os.environ.get(PROJECT_FETCH_KEY_ENV, "")
    if env.strip():
        return env.strip()
    return None


def is_openssl_encrypted_file(path: str) -> bool:
    try:
        with open(path, "rb") as handle:
            return handle.read(8) == OPENSSL_SALTED_MAGIC
    except OSError:
        return False


def _require_openssl() -> str:
    path = shutil.which("openssl")
    if not path:
        raise ProjectCryptoError(
            "openssl is required to encrypt/decrypt restricted project archives; "
            "install OpenSSL and ensure it is on PATH"
        )
    return path


def encrypt_file(
    plain_path: str,
    encrypted_path: str,
    *,
    key: str,
) -> str:
    """Encrypt ``plain_path`` → ``encrypted_path`` (OpenSSL salted AES-256-CBC)."""
    passphrase = str(key or "").strip()
    if not passphrase:
        raise ProjectCryptoError("encryption key must not be empty")
    openssl = _require_openssl()
    plain = os.path.abspath(plain_path)
    out = os.path.abspath(encrypted_path)
    parent = os.path.dirname(out)
    if parent:
        os.makedirs(parent, exist_ok=True)
    try:
        subprocess.run(
            [
                openssl,
                "enc",
                "-aes-256-cbc",
                "-pbkdf2",
                "-salt",
                "-in",
                plain,
                "-out",
                out,
                "-pass",
                f"pass:{passphrase}",
            ],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or b"").decode("utf-8", errors="replace").strip()
        raise ProjectCryptoError(f"openssl encrypt failed: {detail or exc}") from exc
    return out


def decrypt_file(
    encrypted_path: str,
    plain_path: str,
    *,
    key: str,
) -> str:
    """Decrypt OpenSSL salted AES-256-CBC blob to ``plain_path``."""
    passphrase = str(key or "").strip()
    if not passphrase:
        raise ProjectCryptoError("decryption key must not be empty")
    openssl = _require_openssl()
    enc = os.path.abspath(encrypted_path)
    out = os.path.abspath(plain_path)
    parent = os.path.dirname(out)
    if parent:
        os.makedirs(parent, exist_ok=True)
    try:
        subprocess.run(
            [
                openssl,
                "enc",
                "-aes-256-cbc",
                "-pbkdf2",
                "-d",
                "-in",
                enc,
                "-out",
                out,
                "-pass",
                f"pass:{passphrase}",
            ],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or b"").decode("utf-8", errors="replace").strip()
        raise ProjectCryptoError(
            "openssl decrypt failed (wrong key or corrupt archive)"
            + (f": {detail}" if detail else "")
        ) from exc
    return out


def maybe_decrypt_archive(
    archive_path: str,
    *,
    requires_key: bool,
    encryption_scheme: str,
    key: str | None,
) -> str:
    """Return path to a plain archive (decrypting in place via temp file if needed)."""
    scheme = str(encryption_scheme or ENCRYPTION_SCHEME_NONE).strip().lower()
    looks_encrypted = is_openssl_encrypted_file(archive_path)
    need_decrypt = requires_key or scheme not in {"", ENCRYPTION_SCHEME_NONE, "none"} or looks_encrypted

    if not need_decrypt and not looks_encrypted:
        return archive_path

    if scheme in {"", ENCRYPTION_SCHEME_NONE, "none"} and not looks_encrypted:
        return archive_path

    if scheme not in {
        "",
        ENCRYPTION_SCHEME_NONE,
        "none",
        ENCRYPTION_SCHEME_OPENSSL,
        "openssl",
        "openssl-aes-256-cbc-pbkdf2",
    } and looks_encrypted:
        # Unknown scheme but OpenSSL magic — try openssl path.
        pass
    elif scheme not in {
        ENCRYPTION_SCHEME_OPENSSL,
        "openssl",
        "openssl-aes-256-cbc-pbkdf2",
        "",
        ENCRYPTION_SCHEME_NONE,
        "none",
    } and not looks_encrypted:
        raise ProjectCryptoError(f"unsupported encryption scheme: {encryption_scheme}")

    resolved = resolve_fetch_key(key)
    if not resolved:
        raise ProjectCryptoError(
            "this project archive is encrypted; provide --key or set "
            f"{PROJECT_FETCH_KEY_ENV}"
        )

    fd, plain_path = tempfile.mkstemp(suffix=".tar.gz", prefix="jarvis-decrypt-")
    os.close(fd)
    try:
        decrypt_file(archive_path, plain_path, key=resolved)
        return plain_path
    except Exception:
        try:
            os.unlink(plain_path)
        except OSError:
            pass
        raise


__all__ = [
    "ENCRYPTION_SCHEME_NONE",
    "ENCRYPTION_SCHEME_OPENSSL",
    "OPENSSL_SALTED_MAGIC",
    "PROJECT_FETCH_KEY_ENV",
    "ProjectCryptoError",
    "decrypt_file",
    "encrypt_file",
    "is_openssl_encrypted_file",
    "maybe_decrypt_archive",
    "resolve_fetch_key",
]
