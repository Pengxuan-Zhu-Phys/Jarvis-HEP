"""Atomic Worker admission and per-task archive queue delay accounting."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

DISPATCH_GATE = "hep:dispatch:gate"
DISPATCH_WAKE = "hep:dispatch:wake"
DISPATCH_WAIT = "hep:dispatch:wait"
DISPATCH_IDS = "hep:dispatch:ids"
DISPATCH_KEYS = (DISPATCH_GATE, DISPATCH_WAKE, DISPATCH_WAIT, DISPATCH_IDS)

# Redis is the shared clock for queue events. Absolute timestamps stay on Redis.
# Limit credit to the lease even if the controller never comes back to open it.
_CLOCK = """
local gate = KEYS[1]
local state = redis.call('HMGET', gate, 'owner', 'slots', 'start_us', 'until_us',
                         'credit_us', 'stopped')
local stamp = redis.call('TIME')
local now = tonumber(stamp[1]) * 1000000 + tonumber(stamp[2])
local start = tonumber(state[3]) or 0
local deadline = tonumber(state[4]) or 0
local closed = start > 0 and now < deadline
local clock = tonumber(state[5]) or 0
if start > 0 then
    clock = clock + math.max(0, math.min(now, deadline) - start)
end
local function integer(n)
    return string.format('%.0f', n)
end
"""

_WAKE = """
local function wake(list, count)
    local limit = tonumber(state[2]) or 0
    local needed = math.min(count, limit - redis.call('LLEN', list))
    for _ = 1, needed do redis.call('RPUSH', list, '1') end
end
"""

_GATE = _CLOCK + _WAKE + """
if state[1] ~= ARGV[1] or state[6] == '1' then return 0 end
local action = ARGV[2]
if action == 'close' then
    if not closed then
        redis.call('HSET', gate, 'start_us', integer(now), 'credit_us', integer(clock))
        redis.call('HINCRBY', gate, 'revision', 1)
        redis.call('DEL', KEYS[2])
    end
    redis.call('HSET', gate, 'until_us', integer(now + tonumber(ARGV[3])))
else
    redis.call('HSET', gate, 'start_us', 0, 'until_us', 0, 'credit_us', integer(clock))
    if action == 'stop' then redis.call('HSET', gate, 'stopped', 1) end
    redis.call('HINCRBY', gate, 'revision', 1)
    wake(KEYS[2], tonumber(state[2]))
end
return tonumber(redis.call('HGET', gate, 'revision'))
"""

_CLEAR = """
if redis.call('HGET', KEYS[1], 'owner') ~= ARGV[1] then return 0 end
return redis.call('DEL', unpack(KEYS))
"""

_ENQUEUE = _CLOCK + _WAKE + """
if state[1] ~= ARGV[2] or state[6] == '1' then
    return redis.error_reply('archive admission is stopped or belongs to another run')
end
local count = (#ARGV - 2) / 2
for i = 3, #ARGV, 2 do
    local uuid, payload = ARGV[i], ARGV[i + 1]
    local row = redis.call('HGET', KEYS[3], uuid)
    local credit = row and cjson.decode(row)[1] or 0
    redis.call('HSET', KEYS[3], uuid, cjson.encode({credit, clock}))
    redis.call('HSET', KEYS[4], redis.sha1hex(payload), uuid)
    redis.call(ARGV[1], KEYS[5], payload)
end
redis.call('INCRBY', KEYS[6], count)
if not closed and state[6] ~= '1' then wake(KEYS[2], count) end
return count
"""

_CLAIM = _CLOCK + """
if state[1] ~= ARGV[1] or state[6] == '1' then return {'stopped'} end
local existing = redis.call('LINDEX', KEYS[6], 0)
if existing then return {'ok', existing, 0} end
if closed then return {'paused'} end
local payload = redis.call('LMOVE', KEYS[5], KEYS[6], 'LEFT', 'LEFT')
if not payload then return {'empty'} end
redis.call('HINCRBY', KEYS[7], 'running', 1)
local digest = redis.sha1hex(payload)
local uuid = redis.call('HGET', KEYS[4], digest)
if uuid then
    redis.call('HDEL', KEYS[4], digest)
    local row = redis.call('HGET', KEYS[3], uuid)
    if row then
        local credit = cjson.decode(row)
        redis.call('HSET', KEYS[3], uuid,
                   cjson.encode({credit[1] + math.max(0, clock - credit[2])}))
    end
end
return {'ok', payload, 1}
"""

_CREDITS = _CLOCK + """
if state[1] ~= ARGV[1] then return {} end
local result = {}
for i = 2, #ARGV do
    local uuid = ARGV[i]
    local row = redis.call('HGET', KEYS[2], uuid)
    local credit = 0
    if row then
        local v = cjson.decode(row)
        credit = v[1]
        if v[2] then credit = credit + math.max(0, clock - v[2]) end
    end
    result[i - 1] = integer(credit)
end
return result
"""


class _Dispatch:
    """Private RedisQueue mixin. A run opts in before spawning its Workers."""

    _dispatch_epoch: str | None = None
    _dispatch_has_paused = False
    _dispatch_wait_reason = ""
    _dispatch_revision = 0

    def _dispatch_eval(self, source: str, numkeys: int, *args: Any) -> Any:
        self._require_client()
        client = self._ctrl()
        scripts = self.__dict__.setdefault("_dispatch_scripts", {})
        script = scripts.get(source)
        if script is None:
            script = scripts[source] = client.register_script(source)
        return script(keys=args[:numkeys], args=args[numkeys:], client=client)

    def init_archive_admission(self, epoch: str, slots: int) -> None:
        with self._ctrl().pipeline(transaction=True) as pipe:
            pipe.delete(*DISPATCH_KEYS)
            pipe.hset(DISPATCH_GATE, mapping={"owner": epoch, "slots": slots, "revision": 0})
            pipe.execute()
        self._dispatch_epoch = epoch
        self._dispatch_has_paused = False
        self._dispatch_revision = 0

    def set_archive_admission(self, action: str, *, lease_sec: float = 5.0) -> bool:
        if self._dispatch_epoch is None:
            return False
        if action == "close":
            # Set before the Redis transition so a concurrently starting wait
            # either snapshots the clock or starts before the first pause.
            self._dispatch_has_paused = True
        revision = self._dispatch_eval(
            _GATE, 2, DISPATCH_GATE, DISPATCH_WAKE, self._dispatch_epoch,
            action, int(lease_sec * 1000000),
        )
        if revision:
            self._dispatch_revision = int(revision)
        return bool(revision)

    def clear_archive_admission(self) -> None:
        if self._dispatch_epoch is not None:
            self._dispatch_eval(_CLEAR, len(DISPATCH_KEYS), *DISPATCH_KEYS, self._dispatch_epoch)
            self._dispatch_epoch = None
            self._dispatch_has_paused = False

    def _enqueue_with_admission(
        self, direction: str, tasks: Sequence[Mapping[str, Any]], encoded: Sequence[str | bytes]
    ) -> int:
        from jarvishep2.queue.redis_queue import OP_COUNT, TASK_QUEUE

        args = [direction, self._dispatch_epoch]
        for task, payload in zip(tasks, encoded):
            args.extend((str(task["uuid"]), payload))
        return int(self._dispatch_eval(
            _ENQUEUE, 6, *DISPATCH_KEYS, TASK_QUEUE, OP_COUNT.format(kind="task"), *args
        ))

    def _claim_with_admission(self, inflight_key: str) -> dict[str, Any] | None:
        from jarvishep2.queue.redis_queue import SAMPLE_STATS, TASK_QUEUE, _redis_text

        result = self._dispatch_eval(
            _CLAIM, 7, *DISPATCH_KEYS, TASK_QUEUE, inflight_key, SAMPLE_STATS,
            self._dispatch_epoch,
        )
        self._dispatch_wait_reason = _redis_text(result[0])
        return self._task_from_inflight_result(result)

    def _wait_for_admission(self) -> None:
        self._blpop(DISPATCH_WAKE, timeout=1)

    def dispatch_wait_credits(self, uuids: Sequence[str]) -> dict[str, float]:
        if not uuids or not self._dispatch_has_paused:
            return {}
        values = self._dispatch_eval(
            _CREDITS, 2, DISPATCH_GATE, DISPATCH_WAIT, self._dispatch_epoch, *uuids
        )
        return {uuid: int(value) / 1000000 for uuid, value in zip(uuids, values)}

    def _drop_dispatch_wait(self, uuids: Sequence[str]) -> None:
        if self._dispatch_epoch is not None and uuids:
            # A feedback record has already been popped. Keep delivering it
            # even if optional accounting cannot be collected at this instant.
            from jarvishep2.queue.redis_queue import _is_redis_timeout, _warn_control_timeout

            try:
                self._ctrl().hdel(DISPATCH_WAIT, *uuids)
            except Exception as exc:
                if _is_redis_timeout(exc):
                    _warn_control_timeout("drop_dispatch_wait", exc)
                else:
                    from jarvishep2.logging import get_jarvis_logger

                    get_jarvis_logger("redis").warning("dispatch cleanup failed -> %s", exc)
