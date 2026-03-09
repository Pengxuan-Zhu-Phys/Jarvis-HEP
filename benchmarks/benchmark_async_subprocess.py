#!/usr/bin/env python3
from __future__ import annotations

import argparse
import asyncio
import json
import os
import resource
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass, asdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from jarvishep.async_subprocess import AsyncSubprocessScheduler, SubprocessJob, SubprocessRuntimeConfig


@dataclass
class BenchmarkResult:
    runner: str
    workload: str
    tasks: int
    concurrency: int
    duration_sec: float
    throughput_tps: float
    ok: int
    failed: int
    timed_out: int
    stdout_bytes: int
    stderr_bytes: int
    peak_running: int | None
    rss_mb: float
    run_dir: str | None


def rss_mb() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if rss > 10_000_000:
        return rss / (1024.0 * 1024.0)
    return rss / 1024.0


def build_command(python_bin: str, workload: str) -> list[str]:
    if workload == "small":
        return [python_bin, "-c", "import sys; sys.stdout.write('ok\\n')"]
    if workload == "large":
        return [
            python_bin,
            "-c",
            (
                "import sys;"
                "chunk='x'*65536;"
                "[sys.stdout.write(chunk) for _ in range(16)];"
                "sys.stdout.flush();"
                "[sys.stderr.write(chunk) for _ in range(16)];"
                "sys.stderr.flush()"
            ),
        ]
    raise ValueError(f"Unknown workload: {workload}")


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


async def _run_naive_async(tasks: int, concurrency: int, cmd: list[str]) -> tuple[int, int, int, int]:
    sem = asyncio.Semaphore(concurrency)

    async def _one() -> tuple[int, int, int]:
        async with sem:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            out, err = await proc.communicate()
            return int(proc.returncode), len(out), len(err)

    task_list = [asyncio.create_task(_one()) for _ in range(tasks)]
    rows = await asyncio.gather(*task_list)

    ok = sum(1 for rc, _, _ in rows if rc == 0)
    failed = tasks - ok
    out_bytes = sum(o for _, o, _ in rows)
    err_bytes = sum(e for _, _, e in rows)
    return ok, failed, out_bytes, err_bytes


def run_naive(tasks: int, concurrency: int, cmd: list[str], workload: str) -> BenchmarkResult:
    start = time.perf_counter()
    ok, failed, out_bytes, err_bytes = asyncio.run(_run_naive_async(tasks, concurrency, cmd))
    dur = max(1e-9, time.perf_counter() - start)
    return BenchmarkResult(
        runner="naive",
        workload=workload,
        tasks=tasks,
        concurrency=concurrency,
        duration_sec=dur,
        throughput_tps=tasks / dur,
        ok=ok,
        failed=failed,
        timed_out=0,
        stdout_bytes=out_bytes,
        stderr_bytes=err_bytes,
        peak_running=None,
        rss_mb=rss_mb(),
        run_dir=None,
    )


def run_scheduler(tasks: int, concurrency: int, cmd: list[str], workload: str, run_dir: Path) -> BenchmarkResult:
    scheduler = AsyncSubprocessScheduler(
        config=SubprocessRuntimeConfig(
            max_concurrency=concurrency,
            max_pending=max(128, concurrency * 8),
            log_policy="file",
            progress_interval_sec=999.0,
            diagnostics_enabled=False,
        ),
        logger=_NoopLogger(),
        status_path=str(run_dir / "status.jsonl"),
    )

    futures = []
    start = time.perf_counter()
    try:
        for idx in range(tasks):
            futures.append(
                scheduler.submit(
                    SubprocessJob(
                        cmd=cmd,
                        shell=False,
                        task_id=f"{workload}-{idx:06}",
                        log_dir=str(run_dir / "tasks" / f"{idx:06}"),
                    )
                )
            )

        results = [f.result(timeout=300.0) for f in futures]
    finally:
        scheduler.shutdown(wait=True, timeout=120.0)

    dur = max(1e-9, time.perf_counter() - start)
    ok = sum(1 for r in results if r.ok)
    failed = len(results) - ok
    timed_out = sum(1 for r in results if r.timed_out)
    out_bytes = sum(int(r.stdout_bytes) for r in results)
    err_bytes = sum(int(r.stderr_bytes) for r in results)
    snap = scheduler.snapshot()

    return BenchmarkResult(
        runner="scheduler",
        workload=workload,
        tasks=tasks,
        concurrency=concurrency,
        duration_sec=dur,
        throughput_tps=tasks / dur,
        ok=ok,
        failed=failed,
        timed_out=timed_out,
        stdout_bytes=out_bytes,
        stderr_bytes=err_bytes,
        peak_running=int(snap.get("peak_running", 0)),
        rss_mb=rss_mb(),
        run_dir=str(run_dir),
    )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Benchmark Jarvis async subprocess scheduler vs naive asyncio subprocess orchestration.")
    ap.add_argument("--tasks", type=int, default=10000, help="Number of subprocess jobs per scenario")
    ap.add_argument("--concurrency", type=str, default="8,16,32", help="Comma-separated concurrency values")
    ap.add_argument("--workloads", type=str, default="small,large", help="Comma-separated workloads: small,large")
    ap.add_argument("--runner", type=str, default="both", choices=["scheduler", "naive", "both"])
    ap.add_argument("--python", type=str, default=sys.executable, help="Python executable used for synthetic subprocess tasks")
    ap.add_argument("--run-root", type=str, default="", help="Optional output root for scheduler logs")
    ap.add_argument("--cleanup", action="store_true", help="Remove scheduler run directories after benchmarking")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    conc_list = [int(x.strip()) for x in str(args.concurrency).split(",") if x.strip()]
    workloads = [x.strip() for x in str(args.workloads).split(",") if x.strip()]

    run_root = Path(args.run_root).expanduser().resolve() if args.run_root else Path(tempfile.mkdtemp(prefix="jarvis-bench-"))
    run_root.mkdir(parents=True, exist_ok=True)

    results: list[BenchmarkResult] = []
    for workload in workloads:
        cmd = build_command(args.python, workload)
        for conc in conc_list:
            if args.runner in {"naive", "both"}:
                results.append(run_naive(args.tasks, conc, cmd, workload))
            if args.runner in {"scheduler", "both"}:
                scenario_dir = run_root / f"scheduler_{workload}_c{conc}"
                scenario_dir.mkdir(parents=True, exist_ok=True)
                results.append(run_scheduler(args.tasks, conc, cmd, workload, scenario_dir))

    payload = {
        "python": args.python,
        "tasks": args.tasks,
        "concurrency": conc_list,
        "workloads": workloads,
        "results": [asdict(r) for r in results],
    }

    print(json.dumps(payload, indent=2, ensure_ascii=False))

    if args.cleanup:
        shutil.rmtree(run_root, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
