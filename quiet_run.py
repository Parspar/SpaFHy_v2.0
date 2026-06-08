# -*- coding: utf-8 -*-
"""
quiet_run
=========

Thin wrapper used by ``run_pallas.py`` to call the model ``driver`` /
``parallel_driver`` while optionally silencing the very chatty per-timestep
prints emitted by the 2D groundwater solver.

``run_quiet(func, *args, quiet=True, **kwargs)`` simply calls
``func(*args, **kwargs)`` and returns its result. When ``quiet`` is True the
process' ``stdout`` is redirected to the null device for the duration of the
call, so progress spam (solver residuals, per-subset write messages, ...) is
hidden while the return value (the NetCDF path or list of paths) is preserved.

Implemented with ``os.dup2`` on the underlying file descriptor so that prints
from child processes spawned by ``multiprocessing`` (the ``--parallel`` path)
are suppressed too, not just Python-level ``print`` calls in this process.
"""

import os
import sys
from contextlib import contextmanager


@contextmanager
def _suppress_stdout():
    """Redirect fd 1 (stdout) to os.devnull, restoring it on exit."""
    # Flush any buffered Python-level output first.
    sys.stdout.flush()
    saved_fd = os.dup(1)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 1)
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved_fd, 1)
        os.close(devnull_fd)
        os.close(saved_fd)


def run_quiet(func, *args, quiet=True, **kwargs):
    """
    Call ``func(*args, **kwargs)`` and return its result.

    Args:
        func:            Callable to run (e.g. ``driver`` or ``parallel_driver``).
        *args:           Positional arguments forwarded to ``func``.
        quiet (bool):    If True, suppress stdout for the duration of the call.
        **kwargs:        Keyword arguments forwarded to ``func``.

    Returns:
        Whatever ``func`` returns (NetCDF path, or list of paths for parallel).
    """
    if not quiet:
        return func(*args, **kwargs)
    with _suppress_stdout():
        return func(*args, **kwargs)
