"""Tests for ``fedcat cromwell``: config generation, server start/status/stop and the
discover-time check that a local Cromwell server is up.

The ``wb`` CLI and the ``cromwell`` command are replaced by small scripts, so the
process handling runs for real without Cromwell installed.

Run with::

    pytest regenie_workflows/federated_analysis/fedcat/tests/test_cromwell_server.py -v
"""

from __future__ import annotations

import asyncio
import json
import socket
import stat
import sys
from pathlib import Path

import pytest

from fedcat.core.structures import WorkflowBackend
from fedcat.platforms.base import PlatformError
from fedcat.platforms.gcs_cromwell import cromwell
from fedcat.platforms.gcs_cromwell import cromwell_server as cs
from fedcat.platforms.gcs_cromwell.platform import GcsCromwellCfg, GcsCromwellPlatform

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="uses POSIX scripts and signals")

GENERATED = 'include required(classpath("application"))\nbackend { default = "GCPBATCH" }\n'

FAKE_WB = f"""#!{sys.executable}
import sys
from pathlib import Path
args = dict(a.split("=", 1) for a in sys.argv[1:] if a.startswith("--"))
if "--workspace-bucket-name" not in args:
    sys.exit("missing bucket")
Path(args["--dir"], "cromwell.conf").write_text({GENERATED!r} + "# bucket " + args["--workspace-bucket-name"] + "\\n")
"""

# Serves the two endpoints fedcat reads, on the port from -Dwebservice.port, and
# records its arguments so tests can check what it was launched with.
FAKE_CROMWELL = f"""#!{sys.executable}
import json, sys
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
Path(sys.argv[0]).with_suffix(".args").write_text(json.dumps(sys.argv[1:]))
opts = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("-D"))
default = Path(sys.argv[0]).with_suffix(".backend").read_text().strip()
class H(BaseHTTPRequestHandler):
    def do_GET(self):
        body = {{"/engine/v1/status": {{}},
                 "/api/workflows/v1/backends": {{"supportedBackends": [default], "defaultBackend": default}}}}
        self.send_response(200); self.send_header("Content-Type", "application/json"); self.end_headers()
        self.wfile.write(json.dumps(body.get(self.path, {{}})).encode())
    def log_message(self, *a): pass
HTTPServer(("127.0.0.1", int(opts["webservice.port"])), H).serve_forever()
"""


def _script(path: Path, text: str) -> Path:
    """Write an executable script."""
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def _free_port() -> int:
    """Return a free local port."""
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(autouse=True)
def _fast_polling(monkeypatch):
    """Make the server start and stop polls fast."""
    monkeypatch.setattr(cs, "START_POLL_SECONDS", 0.05)
    monkeypatch.setattr(cs, "STOP_POLL_SECONDS", 0.05)


@pytest.fixture
def fake_wb(tmp_path, monkeypatch):
    """Put a fake `wb` on PATH."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _script(bin_dir / "wb", FAKE_WB)
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path(sys.executable).parent}:/usr/bin:/bin")
    return bin_dir


@pytest.fixture
def server_dir(tmp_path):
    """A server directory with a generated `cromwell.conf`."""
    directory = tmp_path / "cromwell"
    directory.mkdir()
    (directory / cs.CONF_NAME).write_text(GENERATED)
    return directory


@pytest.fixture
def fake_cromwell(tmp_path):
    """Factory of fake `cromwell` commands reporting a given default backend."""

    def make(default_backend: str = "GCPBATCH") -> Path:
        script = _script(tmp_path / "cromwell_fake", FAKE_CROMWELL)
        script.with_suffix(".backend").write_text(default_backend)
        return script

    return make


# ── setup ─────────────────────────────────────────────────────────────────────


def test_setup_appends_fedcat_settings_to_the_generated_config(tmp_path, fake_wb):
    conf = cs.write_config(tmp_path / "cromwell", "my-bucket", 300)
    text = conf.read_text()
    assert text.startswith(GENERATED)
    assert "# bucket my-bucket" in text
    assert 'interface = "127.0.0.1"' in text
    assert "enabled = false" in text
    # Cromwell's pinned default helper image was deleted; alpine lacks the python
    # path Batch sets, so the Debian variant is the one that works on Cromwell <=92
    assert 'cloud-sdk-image-url = "gcr.io/google.com/cloudsdktool/google-cloud-cli:slim"' in text
    # a metadata read for a large scatter takes longer than Cromwell's 40s default
    assert "request-timeout = 300s" in text


def test_setup_keeps_machine_paths_out_of_the_config(tmp_path, fake_wb):
    text = cs.write_config(tmp_path / "cromwell", "b", 10).read_text()
    assert str(tmp_path) not in text
    assert "jdbc:hsqldb:file" not in text


def test_setup_refuses_to_overwrite_without_force(tmp_path, fake_wb):
    cs.write_config(tmp_path, "b", 10)
    with pytest.raises(cs.CromwellServerError, match="--force"):
        cs.write_config(tmp_path, "b", 10)
    cs.write_config(tmp_path, "b", 20, force=True)


def test_setup_reports_wb_failures(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _script(bin_dir / "wb", "#!/bin/sh\necho 'not logged in' >&2\nexit 1\n")
    monkeypatch.setenv("PATH", f"{bin_dir}:/usr/bin:/bin")
    with pytest.raises(cs.CromwellServerError, match="not logged in"):
        cs.write_config(tmp_path / "c", "b", 10)


# ── start / status / stop ─────────────────────────────────────────────────────


def test_start_status_stop(server_dir, fake_cromwell):
    port = _free_port()
    script = fake_cromwell()
    try:
        status = cs.start_server(server_dir, port=port, memory="2g", cromwell=str(script), timeout=30)
        healthy, lines = cs.describe(status, port)
        assert healthy, lines
        args = json.loads(script.with_suffix(".args").read_text())
        assert args[-1] == "server"
        assert "-Xmx2g" in args
        assert f"-Dconfig.file={(server_dir / cs.CONF_NAME).resolve()}" in args
        assert f"-Dwebservice.port={port}" in args
        db = next(a for a in args if a.startswith("-Ddatabase.db.url="))
        assert f"jdbc:hsqldb:file:{server_dir.resolve() / 'db' / 'cromwell'};" in db
        assert cs.describe(cs.server_status(port), port)[0]
        # a second start finds the running server rather than launching another
        cs.start_server(server_dir, port=port, cromwell="/nonexistent", timeout=5)
    finally:
        assert cs.stop_server(server_dir, timeout=30)
    assert cs.server_status(port) is None
    assert not cs.stop_server(server_dir)


def test_start_needs_a_config(tmp_path):
    with pytest.raises(cs.CromwellServerError, match="fedcat cromwell setup"):
        cs.start_server(tmp_path, port=_free_port())


def test_start_reports_a_crash_with_the_log(server_dir, tmp_path):
    script = _script(tmp_path / "crash", "#!/bin/sh\necho 'Unable to access jarfile' \nexit 1\n")
    with pytest.raises(cs.CromwellServerError, match="Unable to access jarfile"):
        cs.start_server(server_dir, port=_free_port(), cromwell=str(script), timeout=30)


def test_a_local_default_backend_is_not_accepted(server_dir, fake_cromwell):
    port = _free_port()
    try:
        status = cs.start_server(server_dir, port=port, cromwell=str(fake_cromwell("Local")), timeout=30)
        healthy, lines = cs.describe(status, port)
        assert not healthy and "Local" in lines[-1]
    finally:
        cs.stop_server(server_dir, timeout=30)


def test_no_server_is_reported_with_the_start_command():
    healthy, lines = cs.describe(None)
    assert not healthy and "fedcat cromwell start" in lines[0]


# ── discover-time check ───────────────────────────────────────────────────────


@pytest.fixture
def fresh_client(monkeypatch):
    """Reset the shared Cromwell client."""
    monkeypatch.setattr(cromwell, "_CROMWELL_CLIENT", None)
    monkeypatch.setattr(cromwell, "_CROMWELL_CLIENT_CONFIGURED", False)


def _platform() -> GcsCromwellPlatform:
    """Return a GCS platform on the `local_cromwell` backend."""
    return GcsCromwellPlatform(GcsCromwellCfg(bucket="b", backend=WorkflowBackend.LOCAL_CROMWELL))


def test_connect_fails_fast_without_a_server(fresh_client, monkeypatch):
    monkeypatch.setenv("CROMWELL_URL", f"http://localhost:{_free_port()}")
    with pytest.raises(PlatformError, match="fedcat cromwell start"):
        asyncio.run(_platform().connect(need_submit=True))


def test_connect_accepts_a_running_server(fresh_client, monkeypatch, server_dir, fake_cromwell):
    port = _free_port()
    monkeypatch.setenv("CROMWELL_URL", cs.server_url(port))
    try:
        cs.start_server(server_dir, port=port, cromwell=str(fake_cromwell()), timeout=30)
        asyncio.run(_platform().connect(need_submit=True))
    finally:
        cs.stop_server(server_dir, timeout=30)


def test_connect_rejects_a_local_default_backend(fresh_client, monkeypatch, server_dir, fake_cromwell):
    port = _free_port()
    monkeypatch.setenv("CROMWELL_URL", cs.server_url(port))
    try:
        cs.start_server(server_dir, port=port, cromwell=str(fake_cromwell("Local")), timeout=30)
        with pytest.raises(PlatformError, match="'Local'"):
            asyncio.run(_platform().connect(need_submit=True))
    finally:
        cs.stop_server(server_dir, timeout=30)


def test_dry_run_does_not_need_a_server(fresh_client, monkeypatch):
    monkeypatch.setenv("CROMWELL_URL", f"http://localhost:{_free_port()}")
    asyncio.run(_platform().connect(need_submit=False))


# ── local submission sources ──────────────────────────────────────────────────


def test_local_wdl_and_options_prefer_working_directory(monkeypatch, tmp_path):
    from loguru import logger

    from fedcat.platforms.gcs_cromwell import workflow_setup as ws

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(ws, "_ANNOUNCED_LOCAL", set())
    stem = "meta_analysis_discovery_aou"

    # Nothing local: the repo checkout's WDL and fedcat's default options.
    repo_wdl = ws.resolve_local_wdl(stem)
    assert repo_wdl is not None and repo_wdl.parent.name == "wdls"
    assert json.loads(ws.local_workflow_options()) == ws.DEFAULT_WORKFLOW_OPTIONS
    assert ws.resolve_local_wdl("no_such_workflow") is None

    (tmp_path / f"{stem}.wdl").write_text("version 1.0\n")
    (tmp_path / ws.LOCAL_OPTIONS_NAME).write_text('{"read_from_cache": true}')
    messages: list[str] = []
    sink = logger.add(lambda m: messages.append(m.record["message"]), level="INFO")
    try:
        assert ws.resolve_local_wdl(stem) == Path(f"{stem}.wdl")
        assert ws.resolve_local_wdl(stem) == Path(f"{stem}.wdl")
        assert json.loads(ws.local_workflow_options()) == {"read_from_cache": True}
    finally:
        logger.remove(sink)
    # Announced once each, not on every submission.
    assert sum(f"{stem}.wdl" in m for m in messages) == 1
    assert sum(ws.LOCAL_OPTIONS_NAME in m for m in messages) == 1
