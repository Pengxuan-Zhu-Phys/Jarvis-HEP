#!/usr/bin/env python3
import os
import threading
from typing import Callable

class BucketAllocator:
    def __init__(
        self,
        base_path: str,
        limit: int = 200,
        width: int = 6,
        start_bucket: int = 1,
        on_bucket_sealed: Callable[[str], None] | None = None,
    ):
        self.base_path = base_path
        self.limit = int(limit)
        self.width = int(width)
        self.bucket = int(start_bucket)
        self.count = 0
        self._lock = threading.Lock()
        self._on_bucket_sealed = on_bucket_sealed
        # Deferred archive callback state:
        # - bucket is "sealed" when allocator starts using next bucket id
        # - callback is emitted only when sealed bucket has no in-flight samples
        self._sealed_buckets = set()
        self._sealed_dispatched = set()
        self._inflight_by_bucket = {}

        os.makedirs(self.base_path, exist_ok=True)

    def _bucket_dir(self, bucket_id: int) -> str:
        return os.path.join(self.base_path, f"{int(bucket_id):0{self.width}d}")

    def _parse_bucket_id(self, bucket_dir: str | int | None) -> int | None:
        if bucket_dir is None:
            return None
        if isinstance(bucket_dir, int):
            return int(bucket_dir)
        name = os.path.basename(str(bucket_dir).rstrip(os.sep))
        if not name.isdigit():
            return None
        try:
            return int(name)
        except Exception:
            return None

    def _should_dispatch_locked(self, bucket_id: int) -> bool:
        if self._on_bucket_sealed is None:
            return False
        if int(bucket_id) not in self._sealed_buckets:
            return False
        if int(bucket_id) in self._sealed_dispatched:
            return False
        return int(self._inflight_by_bucket.get(int(bucket_id), 0)) == 0

    def _dispatch_if_ready(self, bucket_id: int) -> bool:
        if self._on_bucket_sealed is None:
            return False
        sealed_dir = self._bucket_dir(bucket_id)
        try:
            self._on_bucket_sealed(sealed_dir)
        except Exception:
            return False
        with self._lock:
            self._sealed_dispatched.add(int(bucket_id))
        return True

    def next_bucket_dir(self) -> str:
        """Return bucket dir like <base_path>/000001, thread-safe."""
        sealed_bucket_id = None
        should_dispatch = False
        with self._lock:
            if self.count >= self.limit:
                sealed_bucket_id = self.bucket
                self._sealed_buckets.add(int(sealed_bucket_id))
                should_dispatch = self._should_dispatch_locked(int(sealed_bucket_id))
                self.bucket += 1
                self.count = 0
            self.count += 1
            bucket_id = self.bucket
            self._inflight_by_bucket[bucket_id] = int(self._inflight_by_bucket.get(bucket_id, 0)) + 1

        bucket_dir = self._bucket_dir(bucket_id)
        os.makedirs(bucket_dir, exist_ok=True)

        if sealed_bucket_id is not None and should_dispatch:
            self._dispatch_if_ready(int(sealed_bucket_id))

        return bucket_dir

    def mark_sample_started(self, bucket_dir: str | int) -> None:
        """Idempotent helper for explicit sample lifecycle integrations."""
        bucket_id = self._parse_bucket_id(bucket_dir)
        if bucket_id is None:
            return
        with self._lock:
            self._inflight_by_bucket[bucket_id] = int(self._inflight_by_bucket.get(bucket_id, 0)) + 1

    def mark_sample_finished(self, bucket_dir: str | int) -> None:
        bucket_id = self._parse_bucket_id(bucket_dir)
        if bucket_id is None:
            return
        should_dispatch = False
        with self._lock:
            cur = int(self._inflight_by_bucket.get(bucket_id, 0))
            if cur <= 1:
                self._inflight_by_bucket.pop(bucket_id, None)
            else:
                self._inflight_by_bucket[bucket_id] = cur - 1
            should_dispatch = self._should_dispatch_locked(bucket_id)
        if should_dispatch:
            self._dispatch_if_ready(bucket_id)

    def seal_current_bucket(self) -> None:
        """Seal current bucket and dispatch callback once in-flight reaches zero."""
        should_dispatch = False
        with self._lock:
            bucket_id = int(self.bucket)
            self._sealed_buckets.add(bucket_id)
            should_dispatch = self._should_dispatch_locked(bucket_id)
        if should_dispatch:
            self._dispatch_if_ready(bucket_id)

    def current_bucket_dir(self) -> str:
        with self._lock:
            bucket_id = self.bucket
        return self._bucket_dir(bucket_id)

    def get_state(self) -> dict:
        with self._lock:
            return {
                "base_path": self.base_path,
                "bucket": self.bucket,
                "count": self.count,
                "limit": self.limit,
                "width": self.width,
            }

    def set_state(self, state: dict):
        with self._lock:
            # base_path 可以不改；如你希望可恢复也能改，就保留这行
            self.base_path = state.get("base_path", self.base_path)
            self.bucket = int(state.get("bucket", self.bucket))
            self.count  = int(state.get("count", self.count))
            self.limit  = int(state.get("limit", self.limit))
            self.width  = int(state.get("width", self.width))
        os.makedirs(self.base_path, exist_ok=True)

    def check_and_update(self) -> None:
        """Inspect existing bucket directories under base_path and update (bucket, count).

        This is used when resuming a run where some buckets already exist on disk.
        Strategy:
          - Find the largest numeric bucket directory name under base_path.
          - Set `self.bucket` to that bucket index.
          - Set `self.count` to the number of immediate subdirectories in that bucket.
          - If the bucket is already full (count >= limit), advance to next bucket and reset count.

        Note: This scan happens only once at allocator creation/resume, not per-sample.
        """
        os.makedirs(self.base_path, exist_ok=True)

        try:
            entries = os.listdir(self.base_path)
        except Exception:
            return

        # bucket dirs are like "000001" (digits only)
        bucket_nums = []
        for d in entries:
            p = os.path.join(self.base_path, d)
            if d.isdigit() and os.path.isdir(p):
                try:
                    bucket_nums.append(int(d))
                except Exception:
                    pass

        if not bucket_nums:
            # nothing on disk yet
            return

        max_bucket = max(bucket_nums)
        max_bucket_dir = os.path.join(self.base_path, f"{max_bucket:0{self.width}d}")

        try:
            # count immediate children (uuid subdirs) in the max bucket
            count = 0
            for name in os.listdir(max_bucket_dir):
                if os.path.isdir(os.path.join(max_bucket_dir, name)):
                    count += 1
        except Exception:
            count = 0

        with self._lock:
            self.bucket = int(max_bucket)
            self.count = int(count)

            # If bucket is already full, advance to next bucket and reset
            if self.count >= self.limit:
                self.bucket += 1
                self.count = 0
