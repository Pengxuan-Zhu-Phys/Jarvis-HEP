#!/usr/bin/env python3
"""Crypto helpers for restricted project packs (used only via Jarvis2 project CLI).

End users should **not** run openssl by hand. Supported UX:

- ``Jarvis2 project fetch NAME --key …``  (decrypt + unpack)
- ``Jarvis2 project pack PATH --encrypt --key …``  (pack then encrypt)
- ``Jarvis2 project encrypt FILE --key …``  (encrypt an existing tarball)

On-disk format matches OpenSSL ``enc -aes-256-cbc -pbkdf2 -salt`` so archives
stay interoperable. Backend order: system ``openssl`` if present, else the
optional ``cryptography`` package (same file format).
"""

from __future__ import annotations

import os
import secrets
import shutil
import subprocess
import tempfile

# OpenSSL salted ciphertext magic (Salted__)
OPENSSL_SALTED_MAGIC = b"Salted__"
ENCRYPTION_SCHEME_OPENSSL = "openssl-aes-256-cbc"
ENCRYPTION_SCHEME_NONE = "none"
PROJECT_FETCH_KEY_ENV = "JARVIS_PROJECT_FETCH_KEY"
# OpenSSL 1.1.1+ / 3.x default for ``enc -pbkdf2`` without -iter
_OPENSSL_PBKDF2_ITERS = 10000


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


def _openssl_bin() -> str | None:
    return shutil.which("openssl")


def _crypto_backend_hint() -> str:
    return (
        "Restricted project crypto needs either `openssl` on PATH "
        "(usual on macOS/Linux) or `pip install cryptography`. "
        "You still only use: Jarvis2 project fetch|pack|encrypt — not raw openssl."
    )


def _encrypt_with_openssl(plain: str, out: str, passphrase: str) -> None:
    openssl = _openssl_bin()
    if not openssl:
        raise ProjectCryptoError("openssl not found")
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
        raise ProjectCryptoError(f"encrypt failed: {detail or exc}") from exc


def _decrypt_with_openssl(enc: str, out: str, passphrase: str) -> None:
    openssl = _openssl_bin()
    if not openssl:
        raise ProjectCryptoError("openssl not found")
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
            "decrypt failed (wrong key or corrupt archive)"
            + (f": {detail}" if detail else "")
        ) from exc


def _encrypt_with_cryptography(plain: str, out: str, passphrase: str) -> None:
    try:
        from cryptography.hazmat.backends import default_backend
        from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
        from cryptography.hazmat.primitives import hashes, padding
        from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
    except ImportError as exc:
        raise ProjectCryptoError("cryptography package not installed") from exc

    salt = secrets.token_bytes(8)
    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=48,  # 32 key + 16 iv
        salt=salt,
        iterations=_OPENSSL_PBKDF2_ITERS,
        backend=default_backend(),
    )
    key_iv = kdf.derive(passphrase.encode("utf-8"))
    key, iv = key_iv[:32], key_iv[32:48]
    with open(plain, "rb") as handle:
        data = handle.read()
    padder = padding.PKCS7(128).padder()
    padded = padder.update(data) + padder.finalize()
    encryptor = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend()).encryptor()
    ciphertext = encryptor.update(padded) + encryptor.finalize()
    with open(out, "wb") as handle:
        handle.write(OPENSSL_SALTED_MAGIC + salt + ciphertext)


def _decrypt_with_cryptography(enc: str, out: str, passphrase: str) -> None:
    try:
        from cryptography.hazmat.backends import default_backend
        from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
        from cryptography.hazmat.primitives import hashes, padding
        from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
    except ImportError as exc:
        raise ProjectCryptoError("cryptography package not installed") from exc

    with open(enc, "rb") as handle:
        blob = handle.read()
    if len(blob) < 16 or blob[:8] != OPENSSL_SALTED_MAGIC:
        raise ProjectCryptoError("not an OpenSSL salted encrypted archive")
    salt = blob[8:16]
    ciphertext = blob[16:]
    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=48,
        salt=salt,
        iterations=_OPENSSL_PBKDF2_ITERS,
        backend=default_backend(),
    )
    key_iv = kdf.derive(passphrase.encode("utf-8"))
    key, iv = key_iv[:32], key_iv[32:48]
    decryptor = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend()).decryptor()
    try:
        padded = decryptor.update(ciphertext) + decryptor.finalize()
        unpadder = padding.PKCS7(128).unpadder()
        data = unpadder.update(padded) + unpadder.finalize()
    except Exception as exc:
        raise ProjectCryptoError("decrypt failed (wrong key or corrupt archive)") from exc
    with open(out, "wb") as handle:
        handle.write(data)


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
    plain = os.path.abspath(plain_path)
    out = os.path.abspath(encrypted_path)
    parent = os.path.dirname(out)
    if parent:
        os.makedirs(parent, exist_ok=True)

    errors: list[str] = []
    if _openssl_bin():
        try:
            _encrypt_with_openssl(plain, out, passphrase)
            return out
        except ProjectCryptoError as exc:
            errors.append(str(exc))
    try:
        _encrypt_with_cryptography(plain, out, passphrase)
        return out
    except ProjectCryptoError as exc:
        errors.append(str(exc))
    raise ProjectCryptoError(
        "Could not encrypt project archive. " + _crypto_backend_hint()
        + (f" Details: {'; '.join(errors)}" if errors else "")
    )


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
    enc = os.path.abspath(encrypted_path)
    out = os.path.abspath(plain_path)
    parent = os.path.dirname(out)
    if parent:
        os.makedirs(parent, exist_ok=True)

    errors: list[str] = []
    if _openssl_bin():
        try:
            _decrypt_with_openssl(enc, out, passphrase)
            return out
        except ProjectCryptoError as exc:
            errors.append(str(exc))
    try:
        _decrypt_with_cryptography(enc, out, passphrase)
        return out
    except ProjectCryptoError as exc:
        errors.append(str(exc))
    raise ProjectCryptoError(
        "Could not decrypt project archive (wrong key, or missing crypto backend). "
        "Use: Jarvis2 project fetch NAME --key YOUR_KEY. "
        + _crypto_backend_hint()
        + (f" Details: {'; '.join(errors)}" if errors else "")
    )


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
            "this project is restricted. Fetch with a key:\n"
            "  Jarvis2 project fetch NAME --key YOUR_KEY\n"
            f"or set {PROJECT_FETCH_KEY_ENV}."
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
