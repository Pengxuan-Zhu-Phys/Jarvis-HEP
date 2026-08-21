"""First-run Redis installation guidance for the Jarvis CLI.

Jarvis V2 can manage a local Redis process when the ``redis-server`` binary is
available, but it must not install operating-system packages on a user's
behalf.  This module therefore performs a small PATH check and, once per user
installation, prints the package-manager command appropriate for the detected
operating system when Redis is missing.
"""

from __future__ import annotations

import os
import platform
import shutil
import threading
from pathlib import Path
from typing import Mapping


_NOTICE_STATE_ENV = "JARVIS_REDIS_NOTICE_STATE"
_DEFAULT_NOTICE_STATE = Path.home() / ".jarvis" / "redis-install-check-v1"
_REDIS_SERVER_BINARIES = ("redis-server", "redis6-server", "valkey-server")
_checked_state_paths: set[str] = set()
_state_lock = threading.Lock()


def find_redis_server() -> str | None:
    """Return a Redis-compatible server executable path, or ``None`` if absent."""

    for executable in _REDIS_SERVER_BINARIES:
        path = shutil.which(executable)
        if path:
            return path
    return None


def _read_os_release(path: Path = Path("/etc/os-release")) -> dict[str, str]:
    """Read Linux distribution metadata without importing a platform package."""

    values: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return values
    for line in lines:
        key, separator, value = line.partition("=")
        if not separator or not key or key.startswith("#"):
            continue
        value = value.strip().strip('"').strip("'")
        values[key.strip()] = value
    return values


def _linux_install_command(
    *,
    distro: Mapping[str, str] | None = None,
    which=shutil.which,
) -> str:
    """Return an install command for a Linux distribution/package manager."""

    metadata = dict(_read_os_release() if distro is None else distro)
    distro_ids = {
        item.strip().lower()
        for item in (
            metadata.get("ID", ""),
            *str(metadata.get("ID_LIKE", "")).replace(",", " ").split(),
        )
        if item.strip()
    }
    if distro_ids & {
        "ubuntu",
        "debian",
        "linuxmint",
        "lmde",
        "pop",
        "elementary",
        "zorin",
        "kali",
        "raspbian",
        "parrot",
        "deepin",
        "mx",
        "devuan",
        "neon",
        "pureos",
        "trisquel",
    }:
        return "sudo apt install -y redis-server"
    if distro_ids & {"amazon", "amzn", "amazonlinux"}:
        if str(metadata.get("VERSION_ID", "")).strip().startswith("2"):
            return "sudo amazon-linux-extras install redis6 -y"
        # Current Amazon Linux guidance favors Valkey; it speaks the Redis
        # protocol and provides a Redis-compatible server.
        return "sudo dnf install -y valkey"
    if distro_ids & {"fedora"}:
        # Fedora's maintained package is now Valkey, with Redis compatibility.
        return "sudo dnf install -y valkey"
    if distro_ids & {
        "rhel",
        "redhat",
        "centos",
        "rocky",
        "almalinux",
        "ol",
        "oracle",
        "scientific",
        "springdale",
        "openmandriva",
    }:
        if which("dnf"):
            return "sudo dnf install -y redis"
        return "sudo yum install -y redis"
    if distro_ids & {"arch", "manjaro", "artix"}:
        # Arch moved the official package from Redis to Valkey.
        return "sudo pacman -S valkey"
    if distro_ids & {"opensuse", "opensuse-leap", "opensuse-tumbleweed", "sles", "suse"}:
        return "sudo zypper install -y redis"
    if distro_ids & {"alpine"}:
        return "sudo apk add redis"
    if distro_ids & {"gentoo"}:
        return "sudo emerge --ask dev-db/redis"
    if distro_ids & {"void"}:
        return "sudo xbps-install -S redis"
    if distro_ids & {"nixos"}:
        return "nix profile install nixpkgs#redis"
    if distro_ids & {"solus"}:
        return "sudo eopkg install redis"
    if distro_ids & {"mageia"}:
        return "sudo urpmi redis"

    # For derivatives with incomplete /etc/os-release metadata, prefer the
    # package manager that is actually present before falling back to apt.
    if which("apt") or which("apt-get"):
        return "sudo apt install -y redis-server"
    if which("dnf"):
        return "sudo dnf install -y redis"
    if which("dnf5"):
        return "sudo dnf5 install -y redis"
    if which("microdnf"):
        return "sudo microdnf install -y redis"
    if which("yum"):
        return "sudo yum install -y redis"
    if which("pacman"):
        return "sudo pacman -S redis"
    if which("zypper"):
        return "sudo zypper install -y redis"
    if which("apk"):
        return "sudo apk add redis"
    if which("emerge"):
        return "sudo emerge --ask dev-db/redis"
    if which("xbps-install"):
        return "sudo xbps-install -S redis"
    if which("eopkg"):
        return "sudo eopkg install redis"
    if which("urpmi"):
        return "sudo urpmi redis"
    if which("nix"):
        return "nix profile install nixpkgs#redis"
    if which("guix"):
        return "guix install redis"
    if which("brew"):
        return "brew install redis"
    if which("pkgin"):
        return "sudo pkgin install redis"
    if which("tdnf"):
        return "sudo tdnf install redis"
    return "Install Redis with your Linux distribution's package manager (redis-server)."


def redis_install_command(
    *,
    system: str | None = None,
    distro: Mapping[str, str] | None = None,
    which=shutil.which,
) -> str:
    """Return a non-destructive Redis installation command for the host OS."""

    normalized = str(system or platform.system()).strip().lower()
    if normalized in {"darwin", "mac", "macos", "osx"}:
        # Homebrew's formula is ``redis``; it supplies the redis-server binary.
        return "brew install redis"
    if normalized == "linux":
        return _linux_install_command(distro=distro, which=which)
    if normalized in {"windows", "win32"}:
        return "winget install --id Redis.Redis -e"
    if normalized in {"freebsd"}:
        return "sudo pkg install redis"
    if normalized in {"openbsd"}:
        return "doas pkg_add redis"
    if normalized in {"dragonfly"}:
        return "sudo pkg install redis"
    return "Install Redis using your operating system's package manager (redis-server)."


def notice_state_path() -> Path:
    """Return the user-level marker used to make the notice a first-run hint."""

    override = os.environ.get(_NOTICE_STATE_ENV, "").strip()
    return Path(override).expanduser() if override else _DEFAULT_NOTICE_STATE


def _claim_first_run(state_path: Path) -> bool:
    """Atomically claim a state path, returning whether this process won it."""

    resolved = str(state_path.expanduser().resolve())
    with _state_lock:
        if resolved in _checked_state_paths:
            return False
        _checked_state_paths.add(resolved)

    try:
        state_path.parent.mkdir(parents=True, exist_ok=True)
        flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
        descriptor = os.open(str(state_path), flags, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write("Jarvis Redis installation check completed.\n")
        return True
    except FileExistsError:
        return False
    except OSError:
        # A read-only home directory should not make an otherwise valid Jarvis
        # command fail. The in-process set still suppresses repeated notices.
        return True


def emit_first_run_redis_notice(*, stream=None) -> bool:
    """Print a one-time Redis install hint after a Jarvis command completes.

    Returns ``True`` when the first-run check was claimed by this invocation,
    including the case where Redis is already installed and no message is
    needed.  Diagnostics go to stderr by default so JSON/stdout command output
    remains machine-readable.
    """

    if not _claim_first_run(notice_state_path()):
        return False
    if find_redis_server():
        return True

    output = stream
    if output is None:
        import sys

        output = sys.stderr
    system = platform.system()
    command = redis_install_command(system=system)
    label = {
        "Darwin": "macOS",
        "Linux": "Linux",
        "Windows": "Windows",
    }.get(system, system or "this operating system")
    print("[Jarvis] Redis server was not found on PATH.", file=output)
    print(f"[Jarvis] Install Redis for {label} with:", file=output)
    print(f"  {command}", file=output)
    return True


__all__ = [
    "emit_first_run_redis_notice",
    "find_redis_server",
    "notice_state_path",
    "redis_install_command",
]
