"""PDF download boundary: offline HTTP fixtures and atomic destination writes."""

import hashlib
import io
import json
from urllib.error import HTTPError, URLError
from unittest.mock import Mock

import pytest

from jarvishep2 import manual_download as download
from jarvishep2.cli.agent import agent_contract
from jarvishep2.cli.client import dispatch, main
from jarvishep2.cli.help import build_parser

PDF = b"%PDF-1.4\n1 0 obj\n<< /Type /Catalog >>\nendobj\n%%EOF\n"


class Response(io.BytesIO):
    def __init__(self, body=PDF, length=None):
        super().__init__(body)
        self.headers = {} if length is None else {"Content-Length": str(length)}


@pytest.fixture
def http(monkeypatch):
    request = Mock(side_effect=lambda *a, **kw: Response(length=len(PDF)))
    monkeypatch.setattr(download, "urlopen", request)
    return request


@pytest.mark.parametrize("release,route,tag", [
    ("latest", "latest/download", "latest"),
    ("2.1.0", "download/v2.1.0", "v2.1.0"),
    ("v2.2.0-rc.1", "download/v2.2.0-rc.1", "v2.2.0-rc.1"),
])
def test_download_public_release_to_directory(tmp_path, http, release, route, tag):
    report = download.download_manual(output=tmp_path, release=release, timeout=12)
    target = tmp_path / download.MANUAL_FILENAME
    assert target.read_bytes() == PDF
    assert report["path"] == str(target)
    assert report["release"] == tag
    assert report["bytes"] == len(PDF)
    assert report["sha256"] == hashlib.sha256(PDF).hexdigest()
    assert report["source"] == f"{download.REPOSITORY_URL}/releases/{route}/{download.MANUAL_FILENAME}"
    request = http.call_args.args[0]
    assert request.full_url == report["source"] and request.get_method() == "GET"
    assert http.call_args.kwargs == {"timeout": 12}
    assert sorted(tmp_path.iterdir()) == [target]


def test_existing_file_requires_force_before_network(tmp_path, http):
    target = tmp_path / "manual.pdf"
    target.write_bytes(b"old manual")
    with pytest.raises(download.ManualDownloadError) as failure:
        download.download_manual(output=target)
    assert failure.value.exit_code == 2
    http.assert_not_called()
    assert target.read_bytes() == b"old manual"
    download.download_manual(output=target, force=True)
    assert target.read_bytes() == PDF
    assert list(tmp_path.iterdir()) == [target]


@pytest.mark.parametrize("body,length", [
    (b"<html>login required</html>", None), (b"", None),
    (b"%PDF-1.4\ntruncated", None), (PDF, len(PDF) + 1),
    (PDF, -1), (PDF, download.MAX_MANUAL_BYTES + 1),
])
def test_invalid_response_preserves_original_and_cleans_temporary(tmp_path, http, body, length):
    target = tmp_path / "manual.pdf"
    target.write_bytes(b"original")
    http.side_effect = lambda *a, **kw: Response(body, length)
    with pytest.raises(download.ManualDownloadError) as failure:
        download.download_manual(output=target, force=True)
    assert failure.value.exit_code == 1
    assert target.read_bytes() == b"original"
    assert list(tmp_path.iterdir()) == [target]


def test_stream_size_limit_without_content_length(tmp_path, http, monkeypatch):
    monkeypatch.setattr(download, "MAX_MANUAL_BYTES", len(PDF) - 1)
    http.side_effect = lambda *a, **kw: Response()
    with pytest.raises(download.ManualDownloadError, match="exceeds"):
        download.download_manual(output=tmp_path)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("failure", [TimeoutError("timed out"), KeyboardInterrupt()])
def test_interrupted_stream_keeps_original_and_cleans_partial(tmp_path, http, failure):
    target = tmp_path / "manual.pdf"
    target.write_bytes(b"original")
    response = Response()
    response.read = Mock(side_effect=[b"%PDF-1.4\npartial", failure])
    http.side_effect = lambda *a, **kw: response
    expected = KeyboardInterrupt if isinstance(failure, KeyboardInterrupt) else download.ManualDownloadError
    with pytest.raises(expected):
        download.download_manual(output=target, force=True)
    assert target.read_bytes() == b"original"
    assert list(tmp_path.iterdir()) == [target]
    assert response.closed


def test_another_writer_wins_without_overwrite(tmp_path, http):
    target = tmp_path / "manual.pdf"

    def fetch(*a, **kw):
        target.write_bytes(b"other writer")
        return Response()

    http.side_effect = fetch
    with pytest.raises(download.ManualDownloadError, match="appeared") as failure:
        download.download_manual(output=target)
    assert failure.value.exit_code == 2
    assert target.read_bytes() == b"other writer"
    assert list(tmp_path.iterdir()) == [target]


@pytest.mark.parametrize("kwargs", [
    {"release": "../main"}, {"release": "https://example.com"},
    {"timeout": 0}, {"timeout": -1}, {"timeout": float("nan")}, {"timeout": float("inf")},
])
def test_invalid_arguments_do_not_contact_network(tmp_path, http, kwargs):
    with pytest.raises(download.ManualDownloadError) as failure:
        download.download_manual(output=tmp_path, **kwargs)
    assert failure.value.exit_code == 2
    http.assert_not_called()


def test_missing_parent_and_symlink_are_rejected(tmp_path, http):
    target = tmp_path / "existing.pdf"
    target.write_bytes(PDF)
    link = tmp_path / "link.pdf"
    link.symlink_to(target)
    for path in (tmp_path / "absent" / "manual.pdf", link):
        with pytest.raises(download.ManualDownloadError) as failure:
            download.download_manual(output=path, force=True)
        assert failure.value.exit_code == 2
    http.assert_not_called()
    assert target.read_bytes() == PDF and link.is_symlink()


@pytest.mark.parametrize("failure,message", [
    (HTTPError("https://github.com", 404, "Not Found", {}, None), "may not be published yet"),
    (HTTPError("https://github.com", 403, "Forbidden", {}, None), "HTTP 403"),
    (URLError("network unavailable"), "network unavailable"),
])
def test_network_errors_are_actionable(tmp_path, http, failure, message):
    http.side_effect = failure
    with pytest.raises(download.ManualDownloadError, match=message) as result:
        download.download_manual(output=tmp_path)
    assert result.value.exit_code == 1
    assert not list(tmp_path.iterdir())


@pytest.fixture
def no_redis(monkeypatch):
    notice = Mock(side_effect=AssertionError("Download must not run Redis onboarding"))
    monkeypatch.setattr("jarvishep2.redis_notice.emit_first_run_redis_notice", notice)
    yield notice
    notice.assert_not_called()


@pytest.mark.parametrize("prefix", [True, False])
def test_cli_json_and_global_placement(tmp_path, http, capsys, no_redis, prefix):
    argv = ["download-manual", "--output", str(tmp_path)]
    argv = ["--json", *argv] if prefix else [*argv, "--json"]
    assert main(argv) == 0
    captured = capsys.readouterr()
    report = json.loads(captured.out)
    assert report["ok"] and report["command"] == "download-manual"
    assert captured.err == ""


@pytest.mark.parametrize("argv", [
    ["download-manual", "--help"], ["--help", "download-manual"],
    ["--json", "download-manual", "-h"],
])
def test_help_does_not_download(argv, http, no_redis, capsys):
    with pytest.raises(SystemExit) as status:
        main(argv)
    assert status.value.code == 0
    assert "--release" in capsys.readouterr().out
    http.assert_not_called()


@pytest.mark.parametrize("argv", [
    ["--json", "--unknown"], ["--json=bad"], ["--json", "--output"],
    ["--json", "--timeout", "nonsense"], ["--json", "--release", "../bad"],
])
def test_cli_usage_errors_are_single_json(argv, http, capsys, no_redis):
    assert main(["download-manual", *argv]) == 2
    captured = capsys.readouterr()
    report = json.loads(captured.out)
    assert not report["ok"] and report["error"]["type"] == "usage_error"
    assert captured.err == ""
    http.assert_not_called()


@pytest.mark.parametrize("failure,status,kind", [
    (URLError("offline"), 1, "download_error"),
    (KeyboardInterrupt(), 130, "interrupted"),
])
def test_cli_failures_report_json(tmp_path, http, capsys, no_redis, failure, status, kind):
    http.side_effect = failure
    assert main(["download-manual", "--json", "-o", str(tmp_path)]) == status
    captured = capsys.readouterr()
    assert json.loads(captured.out)["error"]["type"] == kind
    assert captured.err == ""
    assert not list(tmp_path.iterdir())


def test_human_success_and_errors_use_correct_streams(tmp_path, http, capsys, no_redis, monkeypatch):
    monkeypatch.setenv("NO_COLOR", "1")
    assert main(["download-manual", "-o", str(tmp_path)]) == 0
    captured = capsys.readouterr()
    assert "Downloaded manual:" in captured.out and captured.err == ""
    assert main(["download-manual", "-o", str(tmp_path)]) == 2
    captured = capsys.readouterr()
    assert captured.out == "" and "already exists" in captured.err
    assert "Traceback" not in captured.err and "\x1b[" not in captured.err


def test_root_parser_dispatch_and_agent_catalog(tmp_path, http, capsys):
    parser = build_parser()
    assert "download-manual" in parser.format_help()
    assert dispatch(parser.parse_args(["download-manual", "-o", str(tmp_path), "--json"])) == 0
    assert json.loads(capsys.readouterr().out)["ok"]
    row = next(r for r in agent_contract()["command_catalog"] if r["command"] == "download-manual")
    assert row["effect"] == "network_write" and row["arguments_complete"]
    flags = {flag for arg in row["arguments"] for flag in arg["flags"]}
    assert {"--release", "--output", "--force", "--timeout", "--json"} <= flags


def test_publish_failure_preserves_original(tmp_path, http, monkeypatch):
    target = tmp_path / "manual.pdf"
    target.write_bytes(b"original")
    monkeypatch.setattr(download.os, "replace", Mock(side_effect=PermissionError("read-only destination")))
    with pytest.raises(download.ManualDownloadError, match="read-only destination"):
        download.download_manual(output=target, force=True)
    assert target.read_bytes() == b"original"
    assert list(tmp_path.iterdir()) == [target]


def test_unresolvable_home_is_a_usage_error(http, monkeypatch):
    monkeypatch.setattr(download.Path, "expanduser", Mock(side_effect=RuntimeError("home unavailable")))
    with pytest.raises(download.ManualDownloadError, match="Invalid destination") as failure:
        download.download_manual(output="~/manual.pdf")
    assert failure.value.exit_code == 2
    http.assert_not_called()
