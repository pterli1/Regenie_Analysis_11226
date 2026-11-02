"""Saving run state, of discovery phases and conditioning steps, through the state write queue."""

from __future__ import annotations

import asyncio
from pathlib import Path

from loguru import logger

from fedcat.core.logging import STATE_WRITE_QUEUE, write_text_atomic
from fedcat.core.structures import make_conv
from fedcat.platforms.gcs_cromwell import gcs


class PersistentState:
    """Mixin saving an object's state as JSON to a local file, and optionally uploading it to GCS.

    Writes go through [`STATE_WRITE_QUEUE`][fedcat.core.logging.STATE_WRITE_QUEUE]. The
    file path and content are both taken when the write runs, so call `flush_state`
    before changing anything the path depends on (a new phase or step) to leave the old
    file with its final state.
    """

    __slots__ = ()

    def _state_file(self) -> Path:
        """Return the local state file; implemented by subclasses."""
        raise NotImplementedError

    def _state_queue_key(self) -> str:
        """Return the key the object's writes are queued under; implemented by subclasses."""
        raise NotImplementedError

    def _state_upload_destination(self, path: Path) -> str | None:
        """Return the GCS path (without `gs://`) to upload the state file to, or None to keep it local."""
        return None

    def persist_state(self, reason: str = "checkpoint", max_delay: float | None = None) -> None:
        """Queue a write of the current state; nothing is saved outside a running event loop.

        Args:
            reason: Why the state is saved, for the log.
            max_delay: Most seconds before the write runs; see
                [`StateWriteQueue.enqueue`][fedcat.core.logging.StateWriteQueue.enqueue].
        """
        try:
            asyncio.get_running_loop()
        except RuntimeError:
            return

        async def _write() -> None:
            """Write the state as it is when the write runs."""
            await self._write_state(reason)

        STATE_WRITE_QUEUE.enqueue(self._state_queue_key(), _write, max_delay=max_delay)

    async def flush_state(self) -> None:
        """Write any pending state now, and wait for it and any write already running."""
        await STATE_WRITE_QUEUE.flush(self._state_queue_key())

    async def _write_state(self, reason: str) -> None:
        """Write the state file, then upload it; a failed upload is only warned about."""
        path = self._state_file()
        write_text_atomic(path, make_conv().dumps(self, indent=4))
        logger.debug(f"Saved state to {path} (reason={reason})")
        if destination := self._state_upload_destination(path):
            try:
                await gcs.upload_file_async(path, destination=destination)
                logger.debug(f"Uploaded state to gs://{destination}")
            except Exception as exc:  # noqa: BLE001 — the local state file is written; the upload is best effort
                logger.warning(f"Failed to upload state {path.name} to gs://{destination}: {exc}")
        await logger.complete()
