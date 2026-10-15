"""Redis environment inspection and first-run guidance for the Jarvis CLI.

Jarvis V2 can manage a local Redis process when the ``redis-server`` binary is
available, but it must not install operating-system packages on a user's
behalf. The first-run notice only checks PATH; an explicit environment query can
also run the discovered executable's bounded ``--version`` command. Neither
operation connects to a server, installs packages, or starts a service.
"""

from __future__ import annotations

import os
import platform
import re
import shlex
import shutil
import subprocess
import threading
from pathlib import Path
from typing import Any, Callable, Mapping


_NOTICE_STATE_ENV = "JARVIS_REDIS_NOTICE_STATE"
_DEFAULT_NOTICE_STATE = Path.home() / ".jarvis" / "redis-install-check-v1"
_REDIS_SERVER_BINARIES = ("redis-server", "redis6-server", "valkey-server")
_VERSION_TIMEOUT_SECONDS = 2.0
_VERSION_PATTERN = re.compile(r"\b(Redis|Valkey) server v=(\d+\.\d+(?:\.\d+)?)", re.IGNORECASE)
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
    which=None,
) -> str:
    """Return an install command for a Linux distribution/package manager."""

    which = shutil.which if which is None else which
    metadata = dict(_read_os_release() if distro is None else distro)
    distro_ids = {
        item.strip().lower()
        for item in (
            metadata.get("ID", ""),
            *str(metadata.get("ID_LIKE", "")).replace(",", " ").split(),
        )
        if item.strip()
    }
    # ID is authoritative for a known distribution. In particular Rocky's
    # ID_LIKE includes Fedora, but its package choices follow RHEL's releases.
    primary_id = str(metadata.get("ID", "")).strip().lower()
    known_ids = {
        "ubuntu", "debian", "linuxmint", "lmde", "pop", "elementary",
        "zorin", "kali", "raspbian", "parrot", "deepin", "mx", "devuan",
        "neon", "pureos", "trisquel", "amazon", "amzn", "amazonlinux",
        "fedora", "rhel", "redhat", "centos", "rocky", "almalinux", "ol",
        "oracle", "scientific", "springdale", "openmandriva", "arch",
        "manjaro", "artix", "opensuse", "opensuse-leap", "opensuse-tumbleweed",
        "sles", "suse", "alpine", "gentoo", "void", "nixos", "solus", "mageia",
    }
    for candidate_id in (primary_id, *str(metadata.get("ID_LIKE", "")).replace(",", " ").split()):
        if candidate_id.lower() in known_ids:
            distro_ids = {candidate_id.lower()}
            break
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
        if str(metadata.get("VERSION_ID", "")).strip().split(".")[0] == "2":
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
        version = str(metadata.get("VERSION_ID", "")).strip().split(".")[0]
        if version == "8" and distro_ids & {"rhel", "redhat", "centos", "rocky", "almalinux", "ol", "oracle"}:
            return "sudo dnf module install -y redis:6"
        package = "valkey" if version.isdigit() and int(version) >= 10 else "redis"
        if which("dnf"):
            return f"sudo dnf install -y {package}"
        return f"sudo yum install -y {package}"
    if distro_ids & {"arch", "manjaro", "artix"}:
        # Arch moved the official package from Redis to Valkey.
        return "sudo pacman -S valkey"
    if distro_ids & {"opensuse", "opensuse-leap", "opensuse-tumbleweed", "sles", "suse"}:
        return "sudo zypper install -y redis"
    if distro_ids & {"alpine"}:
        version = str(metadata.get("VERSION_ID", "")).strip()
        pieces = version.split(".")
        legacy = len(pieces) >= 2 and all(piece.isdigit() for piece in pieces[:2]) and (
            int(pieces[0]), int(pieces[1])
        ) < (3, 20)
        package = "redis" if legacy or version == "edge" else "valkey"
        return f"sudo apk add {package}"
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

    # For derivatives with incomplete /etc/os-release metadata, prefer an
    # available package manager. With no evidence, leave a textual instruction.
    if which("apt"):
        return "sudo apt install -y redis-server"
    if which("apt-get"):
        return "sudo apt-get install -y redis-server"
    if which("dnf"):
        return "sudo dnf install -y redis"
    if which("dnf5"):
        return "sudo dnf5 install -y redis"
    if which("microdnf"):
        return "sudo microdnf install -y redis"
    if which("yum"):
        return "sudo yum install -y redis"
    if which("pacman"):
        return "sudo pacman -S valkey"
    if which("zypper"):
        return "sudo zypper install -y redis"
    if which("apk"):
        return "sudo apk add valkey"
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
    which=None,
) -> str:
    """Return a non-destructive Redis installation command for the host OS."""

    normalized = str(system or platform.system()).strip().lower()
    if normalized in {"darwin", "mac", "macos", "osx"}:
        # Homebrew's formula is ``redis``; it supplies the redis-server binary.
        return "brew install redis"
    if normalized == "linux":
        return _linux_install_command(distro=distro, which=which)
    if normalized in {"windows", "win32"}:
        return "wsl --install -d Ubuntu"
    if normalized in {"freebsd"}:
        return "sudo pkg install redis"
    if normalized in {"openbsd"}:
        return "doas pkg_add redis"
    if normalized in {"dragonfly"}:
        return "sudo pkg install redis"
    return "Install Redis using your operating system's package manager (redis-server)."


def _probe_server_version(path: str) -> dict[str, Any]:
    """Inspect an executable without starting Redis or contacting any endpoint."""

    result: dict[str, Any] = {
        "version": None,
        "version_text": None,
        "version_status": "error",
        "version_error": None,
        "supports_blmove": None,
        "usable": False,
    }
    try:
        completed = subprocess.run(
            [path, "--version"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
            timeout=_VERSION_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired:
        result.update(
            version_status="timeout",
            version_error=f"Version query timed out after {_VERSION_TIMEOUT_SECONDS:g} seconds.",
        )
        return result
    except OSError as exc:
        result["version_error"] = f"{type(exc).__name__}: {exc}"
        return result

    output = (str(completed.stdout or "") + "\n" + str(completed.stderr or "")).strip()
    result["version_text"] = output[:2048] or None
    if completed.returncode:
        result["version_error"] = f"Version query exited with status {completed.returncode}."
        return result
    match = _VERSION_PATTERN.search(output)
    if match is None:
        result["version_error"] = "Executable returned an unrecognized Redis/Valkey version."
        return result
    version = match.group(2)
    version_parts = tuple(int(part) for part in version.split(".")[:2])
    compatible = match.group(1).lower() == "valkey" or version_parts >= (6, 2)
    result.update(
        version=version,
        version_status="ok",
        supports_blmove=compatible,
        usable=compatible,
    )
    return result


def _legacy_apt_release(metadata: Mapping[str, str]) -> tuple[str, str | None] | None:
    """Identify explicit releases whose distribution Redis predates BLMOVE.

    Only the named Ubuntu/Debian releases have a verified vendor APT recipe.
    Do not apply a derivative's own VERSION_ID to its parent's repositories.
    """

    distro_id = str(metadata.get("ID", "")).strip().lower()
    version = str(metadata.get("VERSION_ID", "")).strip()
    parts = version.split(".")
    if not version or not all(part.isdigit() for part in parts):
        return None
    if distro_id == "ubuntu" and len(parts) >= 2:
        release = (int(parts[0]), int(parts[1]))
        if release <= (22, 4):
            codename = {(20, 4): "focal", (22, 4): "jammy"}.get(release)
            return f"Ubuntu {version}", codename
    elif distro_id == "debian" and int(parts[0]) <= 11:
        return f"Debian {version}", "bullseye" if int(parts[0]) == 11 else None
    return None


def _vendor_apt_commands(codename: str) -> list[str]:
    """Print the official signed APT recipe for a known release, without running it."""

    # Redis's documented recipe, with an initial refresh/TLS prerequisites and
    # the already-known codename in place of $(lsb_release -cs):
    # https://redis.io/docs/latest/operate/oss_and_stack/install/install-stack/apt/
    return [
        "sudo apt-get update",
        "sudo apt-get install -y ca-certificates curl gpg",
        "curl -fsSL https://packages.redis.io/gpg | sudo gpg --dearmor "
        "-o /usr/share/keyrings/redis-archive-keyring.gpg",
        "sudo chmod 644 /usr/share/keyrings/redis-archive-keyring.gpg",
        'echo "deb [signed-by=/usr/share/keyrings/redis-archive-keyring.gpg] '
        f'https://packages.redis.io/deb {codename} main" '
        "| sudo tee /etc/apt/sources.list.d/redis.list",
        "sudo apt-get update",
        "sudo apt-get install -y redis",
    ]


def build_redis_environment_report(
    *,
    system: str | None = None,
    distro: Mapping[str, str] | None = None,
    which: Callable[[str], str | None] | None = None,
    probe_version: bool = True,
) -> dict[str, Any]:
    """Return fresh local Redis detection and copyable installation guidance.

    ``ok`` describes successful report creation, not scan readiness. The version
    probe only checks the executable; it does not prove a running endpoint has
    the commands, configuration, or permissions a scan needs. Passing
    ``probe_version=False`` performs PATH/OS inspection only.
    """

    which = shutil.which if which is None else which
    host_system = str(system or platform.system()).strip()
    normalized = host_system.lower()
    is_linux = normalized == "linux"
    is_macos = normalized in {"darwin", "mac", "macos", "osx"}
    is_windows = normalized in {"windows", "win32"}
    metadata = dict((_read_os_release() if distro is None else distro) if is_linux else {})
    label = "macOS" if is_macos else "Windows" if is_windows else host_system or "Unknown OS"
    distro_id = str(metadata.get("ID", "")).strip() or None
    distro_name = str(metadata.get("PRETTY_NAME") or metadata.get("NAME") or "").strip() or None
    is_wsl = is_linux and (
        bool(os.environ.get("WSL_DISTRO_NAME"))
        or "microsoft" in platform.release().lower()
    )
    path = None
    executable = None
    for candidate in _REDIS_SERVER_BINARIES:
        candidate_path = which(candidate)
        if candidate_path:
            path, executable = str(candidate_path), candidate
            break
    server: dict[str, Any] = {
        "found": path is not None,
        "path": path,
        "executable": executable,
        "version": None,
        "version_text": None,
        "version_status": "not_checked" if path else "not_found",
        "version_error": None,
        "minimum_version": "6.2",
        "supports_blmove": None,
        "usable": False,
    }
    if path and probe_version:
        server.update(_probe_server_version(path))

    install_command = redis_install_command(system=host_system, distro=metadata, which=which)
    install_commands = [] if install_command.startswith("Install Redis ") else [install_command]
    install_notes = [
        "Package availability and versions depend on your OS release and enabled repositories. "
        "After installation, rerun Jarvis redis check."
    ]
    documentation = ["https://redis.io/docs/latest/operate/oss_and_stack/install/"]
    if not install_commands:
        install_notes.insert(0, install_command)
    if is_macos:
        documentation = ["https://formulae.brew.sh/formula/redis"]
        if not which("brew"):
            install_notes.append("Homebrew was not found on PATH; install Homebrew first: https://brew.sh/.")
    elif is_windows:
        install_commands.extend(["sudo apt update", "sudo apt install -y redis-server"])
        install_notes.append(
            "Run the first command in Administrator PowerShell. Restart if requested, open Ubuntu, "
            "then run the remaining commands in Ubuntu. Install and run Jarvis-HEP inside the same WSL "
            "distribution; installing Redis in WSL does not make it available to native Windows Python."
        )
        documentation = ["https://learn.microsoft.com/en-us/windows/wsl/install"]
    elif is_linux:
        if " apt " in install_command:
            install_commands.insert(0, "sudo apt update")
        elif " apt-get " in install_command:
            install_commands.insert(0, "sudo apt-get update")
        elif " apk " in install_command:
            install_commands.insert(0, "sudo apk update")
        if "valkey" in install_command:
            documentation = ["https://valkey.io/topics/installation/"]
            install_notes.append("Valkey supplies a Redis-compatible server accepted by Jarvis.")
        elif install_command == "sudo dnf module install -y redis:6":
            documentation = [
                "https://docs.redhat.com/en/documentation/red_hat_enterprise_linux/8/html/8.4_release_notes/enhancement_dynamic-programming-languages-web-and-database-servers",
                "https://docs.redhat.com/en/documentation/red_hat_enterprise_linux/8/html/8.7_release_notes/new-features",
            ]
            install_notes.insert(
                0,
                "RHEL 8 and compatible distributions default to an older Redis stream. Select redis:6 "
                "with updated AppStream repositories (RHEL 8.7 rebased it to Redis 6.2.7). Repositories "
                "pinned to older releases may still supply Redis 6.0; update the enabled repositories "
                "first, then rerun Jarvis redis check. If redis:5 is already enabled or installed, "
                "follow your distribution's stream-switching instructions before changing streams.",
            )
        legacy_release = _legacy_apt_release(metadata)
        if legacy_release is not None:
            release_label, codename = legacy_release
            documentation = ["https://redis.io/docs/latest/operate/oss_and_stack/install/install-stack/apt/"]
            if codename:
                install_commands = _vendor_apt_commands(codename)
                install_command = install_commands[-1]
                install_notes.insert(
                    0,
                    f"{release_label}'s default Redis package is below Jarvis's Redis >= 6.2 requirement. "
                    "The commands add Redis's official signed APT repository before installing Redis. "
                    "They change the system's APT key/source files only if you choose to run them.",
                )
            else:
                install_commands = []
                install_command = (
                    "Upgrade to a supported OS release or follow Redis's official installation guide "
                    "to obtain Redis >= 6.2 or compatible Valkey."
                )
                install_notes.insert(
                    0,
                    f"{release_label}'s default Redis package may be too old. No release-specific "
                    "vendor repository recipe has been verified here; an ordinary redis-server "
                    "package install is not a sufficient remedy. " + install_command,
                )

    expected_binary = "valkey-server" if "valkey" in install_command else (
        "redis6-server" if "redis6" in install_command else "redis-server"
    )
    # The optional manual example is separate from Jarvis's scan-owned Redis.
    # Disable persistence so experimenting in a working directory writes no DB.
    start_binary = expected_binary if is_windows else path or expected_binary
    start_command = (
        f"{shlex.quote(start_binary)} --bind 127.0.0.1 --port 6380 "
        "--save '' --appendonly no"
    )
    start_notes = [
        "Jarvis starts and stops its own Redis process for each local scan; a global Redis service "
        "does not need to be started first.",
        "For an optional manual test, run the command in a terminal and stop it with Ctrl-C. "
        "It binds to localhost, uses port 6380, and disables persistence. Choose another free port "
        "if 6380 is in use; Jarvis uses the port selected by the scan configuration.",
    ]
    if is_windows:
        start_notes.append("Run this Linux command inside WSL after installing Redis there.")
    notes = [
        "Detection checks PATH, not the OS package database: a missing executable may already be "
        "installed outside PATH.",
        "Jarvis requires Redis >= 6.2 (BLMOVE) or a compatible Valkey server. "
        "A version check does not connect to or verify a running Redis endpoint.",
        "These commands are instructions only; this query does not install packages or start services.",
    ]
    if is_windows:
        notes.append(
            "Native Windows is not supported by Jarvis's managed Redis launcher, which uses bash "
            "and Unix process signals. Run Jarvis-HEP and Redis inside WSL."
        )
    return {
        "ok": True,
        "kind": "redis-environment",
        "system": {
            "name": host_system,
            "label": label,
            "distro_id": distro_id,
            "distro_name": distro_name,
            "wsl": is_wsl,
        },
        "runtime_supported": False if is_windows else True if is_linux or is_macos else None,
        "server": server,
        "install": {
            "command": install_command,
            "commands": install_commands,
            "notes": install_notes,
            "documentation": documentation,
        },
        "start": {
            "required_for_jarvis": False,
            "command": start_command,
            "commands": [start_command],
            "notes": start_notes,
        },
        "notes": notes,
    }


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


def emit_first_run_redis_notice(*, stream=None, show_hint: bool = True) -> bool:
    """Print a one-time Redis install hint after a Jarvis command completes.

    Returns ``True`` when the first-run check was claimed by this invocation,
    including the case where Redis is already installed and no message is
    needed.  Diagnostics go to stderr by default so JSON/stdout command output
    remains machine-readable.

    ``show_hint=False`` claims the marker without repeating guidance already
    displayed by a successful explicit Redis environment query.
    """

    if not _claim_first_run(notice_state_path()):
        return False
    if not show_hint:
        return True
    if find_redis_server():
        return True

    output = stream
    if output is None:
        import sys

        output = sys.stderr
    report = build_redis_environment_report(probe_version=False)
    label = report["system"]["distro_name"] or report["system"]["label"]
    print("[Jarvis] Redis server was not found on PATH.", file=output)
    print("[Jarvis] It may be installed outside PATH; Jarvis needs Redis >= 6.2 or compatible Valkey.", file=output)
    print(f"[Jarvis] Install Redis for {label} with:", file=output)
    for command in report["install"]["commands"] or [report["install"]["command"]]:
        print(f"  {command}", file=output)
    if report["runtime_supported"] is False:
        print(f"[Jarvis] {report['install']['notes'][-1]}", file=output)
    elif report["system"]["label"] == "macOS" and not shutil.which("brew"):
        print("[Jarvis] Homebrew is required for this command: https://brew.sh/", file=output)
    print("[Jarvis] Jarvis starts/stops Redis for each local scan; no service startup is needed.", file=output)
    print("[Jarvis] Redis details: Jarvis redis | install help: Jarvis redis install | startup help: Jarvis redis start", file=output)
    return True


__all__ = [
    "build_redis_environment_report",
    "emit_first_run_redis_notice",
    "find_redis_server",
    "notice_state_path",
    "redis_install_command",
]
