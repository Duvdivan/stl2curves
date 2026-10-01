"""
Worker processes for the slow steps that split into independent pieces: the surface
passes over separate smooth regions (features.py) and the blend fills (build.py).

OpenCascade, and numpy on the small arrays this program mostly handles, keep hold of
Python's lock while they work, so threads take turns rather than run at once; separate
processes don't. Starting them takes a second or so each on Windows, so they are started
in the background as a conversion begins (start), while the mesh is repaired and its
first surfaces found. Big data they all need (the mesh) is written once to a file
(share) and read by each worker once (load).
"""
import atexit
import os
import pickle
import tempfile
import threading
# (imported before stop is registered to run at exit: exit handlers run last-registered
# first, and multiprocessing's own one waits for every worker to end by itself, which
# workers still waiting for jobs never do)
from concurrent.futures import ProcessPoolExecutor

WORKERS = max(1, min(16, (os.cpu_count() or 2) - 2))     # processes working at once

_pool = None                        # the pool, once its workers are all up
_launching = None                   # the pool while its workers are being launched
_starting = None                    # the thread launching them
_cancel = False                     # stop() was called while they were being launched
_leftovers = False
_files = []
_loaded = {}                        # (in a worker) key -> object, the few read last
memo = {}                           # (in a worker) anything worth keeping between jobs


def start():
    """Start the worker processes in the background (no-op if running or not wanted).
    They are kept for further parts and stopped when the program ends."""
    global _starting, _cancel
    if WORKERS > 1 and _pool is None and _starting is None:
        _cancel = False
        _starting = threading.Thread(target=_start, daemon=True)
        _starting.start()


def _start():
    global _pool, _launching
    try:
        pool = _launching = ProcessPoolExecutor(WORKERS)
        # (on Windows a worker is only launched when a job finds none idle, one at a
        # time, in the middle of sending work out: launch them all now, one by one,
        # minding stop() meanwhile)
        try:
            while len(pool._processes) < WORKERS and not _cancel:
                pool._spawn_process()
        except Exception:
            pass
        if not _cancel:
            pool.submit(int).result()
            _pool = pool
    except Exception:
        _pool = None                # (no worker processes to be had: the work is done here)


def get(wait=True):
    """The pool of worker processes once started, or None. If it is still starting, wait
    for it (wait True) or not (None: small jobs are quicker done in this process)."""
    global _starting
    if _starting is not None:
        if not wait and _starting.is_alive():
            return None
        _starting.join()
        _starting = None
    return _pool


def abandoned():
    """A job was given up while still running (an overrunning fill or fuse): the workers
    are restarted after this part (see finish), not left grinding into the next one."""
    global _leftovers
    _leftovers = True


def finish():
    """At the end of a part: keep the workers for the next one, unless jobs were given up."""
    global _leftovers
    if _leftovers:
        _leftovers = False
        stop()


def broken():
    """The pool failed: stop using it (the work falls back to this process)."""
    global _pool
    pool, _pool = _pool, None
    if pool is not None:
        _terminate(pool)


def stop():
    """Stop the worker processes (any still grinding on an overrunning job included, and
    any still being launched)."""
    global _pool, _launching, _starting, _cancel
    _cancel = True
    if _starting is not None:
        _starting.join()            # (at most the one launch under way)
        _starting = None
    pools = {id(p): p for p in (_pool, _launching) if p is not None}
    _pool = _launching = None
    for pool in pools.values():
        _terminate(pool)
    for path in _files:
        try:
            os.remove(path)
        except OSError:
            pass
    _files.clear()


atexit.register(stop)


def _terminate(pool):
    for p in list((getattr(pool, "_processes", None) or {}).values()):
        p.terminate()
    pool.shutdown(wait=False, cancel_futures=True)


def share(obj):
    """Write obj where every worker can read it (once each): its key for load()."""
    fd, path = tempfile.mkstemp(prefix="stl2curves_", suffix=".pkl")
    with os.fdopen(fd, "wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    _files.append(path)
    return path


def load(key):
    """In a worker: the object shared under this key (read once, then kept)."""
    if key not in _loaded:
        while len(_loaded) >= 3:    # (the mesh and what the current pass needs)
            del _loaded[next(iter(_loaded))]
        with open(key, "rb") as f:
            _loaded[key] = pickle.load(f)
    else:
        _loaded[key] = _loaded.pop(key)     # (most recently used last)
    return _loaded[key]
