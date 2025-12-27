#!/usr/bin/env python3
import os
import threading

class BucketAllocator:
    def __init__(self, base_path: str, limit: int = 200, width: int = 6, start_bucket: int = 1):
        self.base_path = base_path
        self.limit = int(limit)
        self.width = int(width)
        self.bucket = int(start_bucket)
        self.count = 0
        self._lock = threading.Lock()

        os.makedirs(self.base_path, exist_ok=True)

    def next_bucket_dir(self) -> str:
        """Return bucket dir like <base_path>/000001, thread-safe."""
        with self._lock:
            if self.count >= self.limit:
                self.bucket += 1
                self.count = 0
            self.count += 1
            bucket_id = self.bucket

        bucket_dir = os.path.join(self.base_path, f"{bucket_id:0{self.width}d}")
        os.makedirs(bucket_dir, exist_ok=True)
        return bucket_dir

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