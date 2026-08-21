#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Jarvis-HEP V2 UUID channel helpers for the vendored dynesty package.

V1 contract (reused here):
  * ``prior_transform(u)`` may return ``np.append(coords, uuid_str)``
  * Dynesty tracks ``live_uid`` / results ``samples_uid``
  * ``loglikelihood`` receives the **full** vector including the trailing uuid
    so the NestedLikelihoodBridge can stamp Sample.uuid.

Stock dynesty is unchanged when prior_transform returns a plain float vector
of length ``ndim``.
"""

from __future__ import annotations

import numpy as np

# Package-level flag (tests may toggle).
JARVIS_UUID_CHANNEL = True


def looks_uuid_augmented(row) -> bool:
    """Return True if *row* looks like [...coords, uuid_str].

    Only **string-like** trailing elements count. Non-numeric objects such as
    ``numpy.random.SeedSequence`` (dynesty bound-bootstrap args) must **not**
    be treated as uuids — otherwise RedisEvaluationPool routes bootstrap work
    through the logL Redis path and crashes.
    """
    arr = np.asarray(row, dtype=object).reshape(-1)
    if arr.size < 2:
        return False
    last = arr[-1]
    if isinstance(last, (str, bytes, np.str_)):
        text = last.decode("utf-8", errors="replace") if isinstance(last, bytes) else str(last)
        # Empty trailing string is not a uuid channel marker.
        return bool(text.strip())
    # uuid.UUID instances (rare)
    if type(last).__name__ == "UUID":
        return True
    return False


def split_vid_row(row):
    """Split one prior_transform output into (v_float, uid_str|None)."""
    arr = np.asarray(row, dtype=object).reshape(-1)
    if looks_uuid_augmented(arr):
        v = np.asarray(arr[:-1], dtype=float)
        uid = str(arr[-1])
        return v, uid
    return np.asarray(arr, dtype=float), None


def split_vid_batch(vid):
    """Split a batch of prior_transform outputs.

    Returns
    -------
    v : ndarray (n, ndim) float
    uid : ndarray (n,) U36 or None if no uuid channel detected
    full : list of object rows suitable for loglikelihood (with uuid if present)
    """
    rows = [np.asarray(r, dtype=object).reshape(-1) for r in vid]
    if not rows:
        return np.zeros((0, 0)), None, []
    if not looks_uuid_augmented(rows[0]):
        v = np.asarray([np.asarray(r, dtype=float) for r in rows], dtype=float)
        return v, None, list(v)
    v = np.asarray([np.asarray(r[:-1], dtype=float) for r in rows], dtype=float)
    uid = np.asarray([str(r[-1]) for r in rows], dtype="U36")
    full = [np.asarray(r, dtype=object) for r in rows]
    return v, uid, full


def join_v_uid(v, uid):
    """Recombine physical/unit coords with uuid for loglikelihood bridge."""
    if uid is None:
        return np.asarray(v, dtype=float)
    out = np.empty(np.asarray(v).size + 1, dtype=object)
    out[:-1] = np.asarray(v, dtype=float)
    out[-1] = str(uid)
    return out
