"""Helpers for the boilerplate around calling the Fortran ``*_c`` entry points.

Every C-callable GENESIS wrapper reports failures through a trailing
``(status, msg, msglen)`` triple and writes progress to stdout/stderr, and the
copy-based wrappers additionally need an ``s_molecule`` handed over and released
explicitly. Centralising that here keeps the analysis functions focused on their
own arguments and makes it impossible to forget the cleanup.
"""

import contextlib
import ctypes
import io
from functools import lru_cache
from typing import Tuple

from .exceptions import raise_fortran_error
from .libgenesis import LibGenesis
from .output_capture import suppress_stdout_capture_stderr

# Fallback when the library does not expose its own message buffer length.
DEFAULT_MSG_LEN = 2048


@lru_cache(maxsize=1)
def get_msg_len() -> int:
    """Return the message buffer length the loaded library expects."""
    lib = LibGenesis().lib
    fn = (getattr(lib, "get_default_msg_len", None)
          or getattr(lib, "genesis_get_default_msg_len", None))
    if fn is None:
        return DEFAULT_MSG_LEN
    fn.restype = ctypes.c_int
    try:
        n = int(fn())
    except Exception:
        return DEFAULT_MSG_LEN
    return n if n > 0 else DEFAULT_MSG_LEN


def make_msgbuf():
    """Allocate an error message buffer sized for the loaded library."""
    n = get_msg_len()
    return (ctypes.c_char * n)(), n


@contextlib.contextmanager
def fortran_status():
    """Provide the ``(status, msg, msglen)`` triple a ``*_c`` call reports through.

    Fortran output is suppressed while the call runs, and a non-zero status is
    turned into the matching GenesisFortranError subclass once it returns. An
    exception raised inside the block propagates unchanged so that the original
    failure is never masked by the status check.

    Usage:
        with fortran_status() as (status, msg, msglen):
            lib.some_analysis_c(..., ctypes.byref(status), msg,
                                ctypes.c_int(msglen))
    """
    msg, msglen = make_msgbuf()
    status = ctypes.c_int(0)

    with suppress_stdout_capture_stderr() as captured:
        yield status, msg, msglen

    if status.value != 0:
        raise_fortran_error(
            msg.value.decode("utf-8", "replace").strip(),
            code=status.value,
            stderr_output=captured.stderr,
        )


@contextlib.contextmanager
def molecule_c(molecule):
    """Hand an SMolecule to Fortran as a copy and always release it again.

    The copy is allocated by Fortran (see the zerocopy notes in CLAUDE.md), so
    it has to be freed explicitly rather than by the Python GC.
    """
    mol_c = molecule.to_SMoleculeC()
    try:
        yield mol_c
    finally:
        LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


def ctrl_to_bytes(ctrl: io.BytesIO) -> Tuple[bytes, ctypes.c_int]:
    """Return a built control file as the ``(text, length)`` pair Fortran wants."""
    data = ctrl.getvalue()
    return data, ctypes.c_int(len(data))
