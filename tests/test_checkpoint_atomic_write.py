"""StateSaver file durability and preservation of the previous recovery point."""

import os
import pickle
import stat

import pytest

from jarvishep2.sampling import runtime_checkpoint as checkpoint


def test_checkpoint_fsync_precedes_atomic_replace_and_syncs_directory(
    tmp_path, monkeypatch
):
    path = tmp_path / "state.pkl"
    events = []
    real_fsync, real_replace = os.fsync, os.replace

    def fsync(fd):
        events.append("directory" if stat.S_ISDIR(os.fstat(fd).st_mode) else "file")
        real_fsync(fd)

    def replace(source, destination):
        assert events == ["file"]
        assert not path.exists()
        events.append("replace")
        real_replace(source, destination)

    monkeypatch.setattr(checkpoint.os, "fsync", fsync)
    monkeypatch.setattr(checkpoint.os, "replace", replace)
    checkpoint.atomic_pickle_dump(str(path), {"cursor": 12})
    assert pickle.loads(path.read_bytes()) == {"cursor": 12}
    assert events == (
        ["file", "replace", "directory"] if os.name == "posix" else ["file", "replace"]
    )
    assert list(tmp_path.glob("*.tmp")) == []


def test_failed_state_save_preserves_previous_file_and_cleans_temp(subtests, tmp_path_factory):
    def fail(*args, **kwargs):
        raise OSError("injected save failure")

    stages = {
        "serialize": (checkpoint.pickle, "dump"),
        "fsync": (checkpoint.os, "fsync"),
        "replace": (checkpoint.os, "replace"),
    }
    for stage, (owner, name) in stages.items():
        with subtests.test(stage=stage), pytest.MonkeyPatch.context() as monkeypatch:
            tmp_path = tmp_path_factory.mktemp("case")
            path = tmp_path / "state.pkl"
            checkpoint.atomic_pickle_dump(str(path), {"cursor": 12})
            before = path.read_bytes()
            monkeypatch.setattr(owner, name, fail)
            with pytest.raises(OSError, match="injected save failure"):
                checkpoint.atomic_pickle_dump(str(path), {"cursor": 24})
            assert path.read_bytes() == before
            assert list(tmp_path.iterdir()) == [path]
