"""Logging setup, the rate-limited state write queue, and limits on concurrent CPU work."""

from __future__ import annotations

import asyncio
import contextlib
import contextvars
import datetime
import inspect
import os
import sys
import time
from collections.abc import AsyncIterator, Awaitable, Callable
from pathlib import Path
from typing import Any

from loguru import logger

# ── CPU / thread pool sizing ────────────────────────────────────────────────
# CPUs fedcat uses, at most 16
NCPU = min(os.cpu_count() or 1, 16)
POLARS_THREADS = max(1, NCPU - 1)
os.environ.setdefault("POLARS_MAX_THREADS", str(POLARS_THREADS))

# ── Module-level state ──────────────────────────────────────────────────────
_LOGGING_CONFIGURED = False
_CONFIGURED_CONTEXTS: set[str] = set()
# Run a log line belongs to, which selects that run's log files
RUN_CONTEXT_VAR: contextvars.ContextVar[str | None] = contextvars.ContextVar("run_context", default=None)

# What each log line is about, shown as its [covar_stratum stage phenotype stepN] prefix.
# Set where that scope starts; tasks started inside it inherit the values.
COVAR_STRATUM_VAR: contextvars.ContextVar[str | None] = contextvars.ContextVar("log_covar_stratum", default=None)
STAGE_VAR: contextvars.ContextVar[str | None] = contextvars.ContextVar("log_stage", default=None)
PHENO_VAR: contextvars.ContextVar[str | None] = contextvars.ContextVar("log_pheno", default=None)
# The running ConditioningStep, read when a line is logged so the step number stays current.
COND_STEP_VAR: contextvars.ContextVar[Any] = contextvars.ContextVar("log_cond_step", default=None)


def _log_label() -> str:
    """Return the `[covar_stratum stage phenotype stepN] ` prefix of the current context."""
    step = COND_STEP_VAR.get()
    parts = [
        COVAR_STRATUM_VAR.get(),
        STAGE_VAR.get(),
        PHENO_VAR.get(),
        f"step{step.step_number}" if step is not None else None,
    ]
    label = " ".join(p for p in parts if p)
    return f"[{label}] " if label else ""


@contextlib.contextmanager
def log_context(var: contextvars.ContextVar, value: Any):
    """Set a log context variable for the `with` block only, e.g. for one phenotype of a loop.

    Args:
        var: Context variable, e.g. `PHENO_VAR`.
        value: Value for the block.
    """
    token = var.set(value)
    try:
        yield
    finally:
        var.reset(token)


def _patch_record(record) -> None:
    """Add the run context and line prefix to a log record."""
    record["extra"].update(context=RUN_CONTEXT_VAR.get(), label=_log_label())


_FILE_FORMAT = "{time:YYYY-MM-DD HH:mm:ss.SSS} | {level: <8} | {name}:{function}:{line} - {extra[label]}{message}"

# ── Async semaphores ─────────────────────────────────────────────────────────
_THREAD_SEMAPHORE = asyncio.BoundedSemaphore(NCPU)
_PERMITS_RESERVATION_LOCK = asyncio.Lock()


# ── Permits context-manager ──────────────────────────────────────────────────


class Permits:
    """Async context manager holding CPU slots from a shared pool of `NCPU`, to limit concurrent CPU work.

    ```python
    async with Permits(4):
        await asyncio.to_thread(heavy_cpu_work)
    ```
    """

    def __init__(self, permits: int | None = None) -> None:
        """Set the number of slots to hold.

        Args:
            permits: Slots to hold, capped at `NCPU`; None for all of them.

        Raises:
            ValueError: If `permits` is below 1.
        """
        if permits is None:
            permits = NCPU
        elif permits < 1:
            raise ValueError("permits must be >= 1")
        elif permits > NCPU:
            permits = NCPU
        self.permits = permits
        self._acquired = 0

    async def __aenter__(self) -> None:
        """Wait for the slots, taking them all before another `Permits` can start taking any."""
        async with _PERMITS_RESERVATION_LOCK:
            try:
                for _ in range(self.permits):
                    await _THREAD_SEMAPHORE.acquire()
                    self._acquired += 1
                logger.trace(f"Permits acquired: {self._acquired}/{self.permits}")
            except Exception:
                for _ in range(self._acquired):
                    _THREAD_SEMAPHORE.release()
                self._acquired = 0
                raise

    async def __aexit__(self, exc_type, exc, tb) -> None:
        """Release the slots."""
        for _ in range(self._acquired):
            _THREAD_SEMAPHORE.release()


# ── Periodic progress reporter ────────────────────────────────────────────────


@contextlib.asynccontextmanager
async def periodic_progress(
    message_fn: Callable[[], str | Awaitable[str]],
    interval: float = 300.0,
) -> AsyncIterator[None]:
    """Log a progress line at INFO every `interval` seconds while the block runs.

    The first line comes after `interval` seconds, so short blocks are silent. Errors
    in `message_fn` are logged at DEBUG and do not affect the block.

    ```python
    async with periodic_progress(lambda: f"{done}/{total} done"):
        await asyncio.gather(*tasks)
    ```

    Args:
        message_fn: Returns the line to log (sync or async); an empty line is skipped.
        interval: Seconds between lines.
    """

    async def _reporter() -> None:
        """Log `message_fn`'s line every `interval` seconds until cancelled."""
        while True:
            await asyncio.sleep(interval)
            try:
                msg = message_fn()
                if inspect.isawaitable(msg):
                    msg = await msg
            except Exception as exc:  # noqa: BLE001 — never let reporting break the run
                logger.debug(f"progress reporter error: {exc}")
                continue
            if msg:
                logger.info(msg)

    task = asyncio.create_task(_reporter(), name="progress-reporter")
    try:
        yield
    finally:
        task.cancel()
        with contextlib.suppress(asyncio.CancelledError):
            await task


# ── Logging configuration ────────────────────────────────────────────────────


def setup_run_logging(
    run_context: str | None = None,
    base_log_path=None,
    log_warnings_separately: bool = True,
    logging_level: str = "INFO",
    covar_stratum: str | None = None,
) -> None:
    """Configure the console log, and the log files of a run.

    The console handler is installed on the first call; a run's log files are added
    the first time its `run_context` is seen. Log files record at least DEBUG,
    whatever the console level.

    Args:
        run_context: Run identifier; log lines logged in its context go to its files.
        base_log_path: Directory for the run's log files; None for no files.
        log_warnings_separately: Also write warnings to a `*.warnings.log` file.
        logging_level: Console level: `"TRACE"`, `"DEBUG"` or `"INFO"`. The caller maps
            `--logging-level NOFILE` to `"TRACE"` with no `base_log_path`.
        covar_stratum: Covariate stratum shown in each line's prefix, for a
            covariate-stratified run.
    """
    global _LOGGING_CONFIGURED

    if not _LOGGING_CONFIGURED:
        console_format = (
            "<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | "
            "<cyan>{extra[label]}</cyan><level>{message}</level>"
        )
        logger.configure(
            handlers=[
                {
                    "sink": sys.stderr,
                    "format": console_format,
                    "level": logging_level,
                    "enqueue": True,
                }
            ],
            patcher=_patch_record,
        )
        _LOGGING_CONFIGURED = True

    if covar_stratum is not None:
        COVAR_STRATUM_VAR.set(covar_stratum)

    if run_context:
        RUN_CONTEXT_VAR.set(run_context)

        if base_log_path and run_context not in _CONFIGURED_CONTEXTS:
            from pathlib import Path as _Path

            base_log_path = _Path(base_log_path)
            ts = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            log_filter = lambda record: record["extra"].get("context") == run_context
            base_log_path.mkdir(parents=True, exist_ok=True)

            file_level = "TRACE" if logging_level == "TRACE" else "DEBUG"
            info_log = base_log_path / f"{run_context}_{ts}.log"
            logger.add(sink=info_log, filter=log_filter, enqueue=True, level=file_level, format=_FILE_FORMAT)
            _CONFIGURED_CONTEXTS.add(run_context)
            logger.info(f"Logging to file ({file_level}+): {info_log}")

            if log_warnings_separately:
                warn_log = base_log_path / f"{run_context}_{ts}.warnings.log"
                logger.add(sink=warn_log, filter=log_filter, level="WARNING", enqueue=True, format=_FILE_FORMAT)
                logger.info(f"Warnings for {run_context} also saved in: {warn_log}")


# ── StateWriteQueue ──────────────────────────────────────────────────────────


WriteFn = Callable[[], Awaitable[None]]


class StateWriteQueue:
    """Rate-limited queue of state writes (and uploads), keeping only the latest per key.

    A key's writes run at most once per `min_interval`, one at a time and in the order
    started, so a slow upload never lands after a newer write. Different keys write
    concurrently.

    ```python
    STATE_WRITE_QUEUE.start()
    STATE_WRITE_QUEUE.enqueue(key, write_fn, max_delay=10)
    await STATE_WRITE_QUEUE.flush(key)
    ```
    """

    def __init__(self, min_interval: float = 120.0) -> None:
        """Create an unstarted queue.

        Args:
            min_interval: Least seconds between writes of one key, unless a write sets
                a shorter `max_delay`.
        """
        self._min_interval = min_interval
        self._pending: dict[str, tuple[WriteFn, float | None]] = {}
        self._last_write: dict[str, float] = {}
        self._key_locks: dict[str, asyncio.Lock] = {}
        self._inflight: dict[str, set[asyncio.Task]] = {}
        self._wake: asyncio.Event | None = None
        self._worker_task: asyncio.Task | None = None

    def start(self) -> None:
        """Start the background worker; must be called in a running event loop."""
        if self._worker_task is None or self._worker_task.done():
            # Recreated so the queue can be restarted in a new event loop
            self._wake = asyncio.Event()
            self._key_locks.clear()
            self._worker_task = asyncio.get_running_loop().create_task(self._worker(), name="StateWriteQueue")

    async def stop(self) -> None:
        """Stop the worker, then run every pending write and wait for all writes to finish."""
        if self._worker_task and not self._worker_task.done():
            self._worker_task.cancel()
            with contextlib.suppress(asyncio.CancelledError):
                await self._worker_task
        # enqueue() refuses until restarted, as nothing would run the write
        self._wake = None
        await self.flush()

    def enqueue(self, key: str, write_fn: WriteFn, max_delay: float | None = None) -> None:
        """Queue a write, replacing any write pending for the same key.

        Args:
            key: What is written, e.g. one run's state file.
            write_fn: Coroutine function performing the write.
            max_delay: Most seconds before the write runs: None to wait up to
                `min_interval` since the key's last write, 0 for as soon as possible.
                An earlier deadline already pending for the key is kept.

        Raises:
            RuntimeError: If the queue is not started.
        """
        if self._wake is None:
            raise RuntimeError("StateWriteQueue.start() must be called before enqueue()")
        deadline = time.monotonic() + max_delay if max_delay is not None else None
        if (existing := self._pending.get(key)) is not None and existing[1] is not None:
            deadline = existing[1] if deadline is None else min(existing[1], deadline)
        self._pending[key] = (write_fn, deadline)
        self._wake.set()

    def _launch(self, key: str, write_fn: WriteFn) -> asyncio.Task:
        """Start a write as a task, tracked until it finishes."""
        self._last_write[key] = time.monotonic()
        task = asyncio.get_running_loop().create_task(self._run(key, write_fn), name=f"state-write:{key}")
        tasks = self._inflight.setdefault(key, set())
        tasks.add(task)
        task.add_done_callback(tasks.discard)
        return task

    async def _run(self, key: str, write_fn: WriteFn) -> None:
        """Run a write after the key's earlier writes, logging any error."""
        async with self._key_locks.setdefault(key, asyncio.Lock()):
            try:
                await write_fn()
            except Exception as e:  # noqa: BLE001 — a failed state write must not stop the queue
                logger.error(f"State write failed for {key!r}: {e}")

    async def flush(self, key: str | None = None) -> None:
        """Run pending writes now, ignoring the rate limit, and wait for them and any running writes.

        Args:
            key: Key to flush; None for every key.
        """
        keys = [key] if key is not None else list(self._pending.keys() | self._inflight.keys())
        for k in keys:
            if k in self._pending:
                self._launch(k, self._pending.pop(k)[0])
        tasks = {t for k in keys for t in self._inflight.get(k, ())}
        if tasks:
            await asyncio.wait(tasks)

    async def _worker(self) -> None:
        """Start each pending write once due, sleeping until the next is due."""
        assert self._wake is not None
        while True:
            await self._wake.wait()
            self._wake.clear()
            now = time.monotonic()
            soonest = float("inf")
            for key, (write_fn, deadline) in list(self._pending.items()):
                due = self._last_write.get(key, float("-inf")) + self._min_interval
                if deadline is not None:
                    due = min(due, deadline)
                if now >= due:
                    del self._pending[key]
                    # Not awaited, so cancelling the worker does not cancel the write
                    self._launch(key, write_fn)
                else:
                    soonest = min(soonest, due - now)
            if soonest != float("inf"):
                asyncio.get_running_loop().call_later(max(soonest, 0.1), self._wake.set)


STATE_WRITE_QUEUE: StateWriteQueue = StateWriteQueue(min_interval=20.0)
"""Shared state write queue; `start()` it in the running event loop before use."""


def write_text_atomic(path: Path, text: str) -> None:
    """Write a text file through a temporary file, so an interrupted write never leaves it truncated.

    Args:
        path: File to write; its directory is created if needed.
        text: Content.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    os.replace(tmp, path)
