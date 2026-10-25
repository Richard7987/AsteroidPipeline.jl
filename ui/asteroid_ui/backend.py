"""Client for the analysis engine (ui/backend/server.jl).

The engine is a Julia process kept open for the whole session: loading and
compiling AsteroidPipeline takes ~20 s, so it is paid once. It speaks JSON,
one message per line (see server.jl).
"""

from __future__ import annotations

import asyncio
import itertools
import json
from collections.abc import Awaitable, Callable
from typing import Any

from asteroid_ui import engine
from asteroid_ui.config import BACKEND_LOG

Event = dict[str, Any]


class BackendError(Exception):
    """The engine answered with an error, or could not be started."""


class Backend:
    def __init__(self) -> None:
        self._proc: asyncio.subprocess.Process | None = None
        self._ids = itertools.count(1)
        self._queues: dict[int, asyncio.Queue[Event]] = {}
        self._lock = asyncio.Lock()  # the engine handles one command at a time
        self.ready = asyncio.Event()
        self.failed: str | None = None
        self.on_status: Callable[[str], None] = lambda message: None

    @property
    def busy(self) -> bool:
        return self._lock.locked()

    async def start(self) -> None:
        """Find Julia, set the engine up if this is the first run (or an upgrade), and start it."""
        try:
            julia = await asyncio.to_thread(engine.julia_executable)
            if await asyncio.to_thread(engine.needs_setup):
                self.on_status("Installing the analysis engine — the first time only, it takes several minutes…")
                loop = asyncio.get_running_loop()

                def progress(line: str) -> None:
                    loop.call_soon_threadsafe(self.on_status, f"Installing the analysis engine: {line[:90]}")

                await asyncio.to_thread(engine.setup, julia, progress)
                self.on_status("Starting the analysis engine…")
        except engine.EngineError as e:
            self.failed = str(e)
            self.ready.set()
            return
        BACKEND_LOG.parent.mkdir(parents=True, exist_ok=True)
        with open(BACKEND_LOG, "a", buffering=1) as log:  # the child keeps its own copy
            self._proc = await asyncio.create_subprocess_exec(
                *engine.server_command(julia),
                stdin=asyncio.subprocess.PIPE, stdout=asyncio.subprocess.PIPE, stderr=log,
                limit=64 * 1024 * 1024,
            )
        asyncio.create_task(self._read_loop())

    async def _read_loop(self) -> None:
        assert self._proc and self._proc.stdout
        while line := await self._proc.stdout.readline():
            try:
                event = json.loads(line)
            except ValueError:
                continue
            if event.get("event") == "ready":
                self.ready.set()
                continue
            if event.get("event") == "status":
                self.on_status(event.get("message", ""))
            queue = self._queues.get(event.get("id"))
            if queue is not None:
                await queue.put(event)
        code = await self._proc.wait()
        self.failed = f"The analysis engine stopped (exit code {code}). Details in {BACKEND_LOG}"
        self.ready.set()
        for queue in self._queues.values():
            await queue.put({"event": "error", "message": self.failed})

    async def request(self, cmd: str, on_event: Callable[[Event], Awaitable[None] | None] | None = None,
                      **args: Any) -> Event:
        """Send a command and wait for its final event (`result`, `report`, `report_check`, `inspection` or `error`).

        Events along the way (progress, status) go to `on_event`.
        """
        await self.ready.wait()
        if self.failed:
            raise BackendError(self.failed)
        async with self._lock:
            assert self._proc and self._proc.stdin
            request_id = next(self._ids)
            queue: asyncio.Queue[Event] = asyncio.Queue()
            self._queues[request_id] = queue
            try:
                self._proc.stdin.write((json.dumps({"cmd": cmd, "id": request_id, **args}) + "\n").encode())
                await self._proc.stdin.drain()
                while True:
                    event = await queue.get()
                    kind = event.get("event")
                    if kind == "error":
                        raise BackendError(event.get("message", "unknown error"))
                    if kind in ("result", "report", "report_check", "inspection"):
                        return event
                    if on_event is not None:
                        maybe = on_event(event)
                        if asyncio.iscoroutine(maybe):
                            await maybe
            finally:
                del self._queues[request_id]

    async def stop(self) -> None:
        if self._proc and self._proc.returncode is None:
            try:
                assert self._proc.stdin
                self._proc.stdin.write(b'{"cmd": "quit"}\n')
                await self._proc.stdin.drain()
                await asyncio.wait_for(self._proc.wait(), 3)
            except (TimeoutError, BrokenPipeError, ConnectionResetError):
                self._proc.kill()
