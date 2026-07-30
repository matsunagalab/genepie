"""Helpers shared by more than one analysis wrapper."""

import os
import numpy as np
from typing import (
    List,
    Optional,
    Tuple,
)

from ..constants import (
    FITTING_METHOD_MAP,
    PBCC_MODE_MAP,
    TRJ_FORMAT_MAP,
    TRJ_TYPE_MAP,
)
from ..exceptions import GenesisValidationError


_TRJ_FORMAT_MAP = TRJ_FORMAT_MAP
_TRJ_TYPE_MAP = TRJ_TYPE_MAP
_FITTING_METHOD_MAP = FITTING_METHOD_MAP
_PBCC_MODE_MAP = PBCC_MODE_MAP


def _resolve_enum(value: str, mapping: dict, label: str) -> int:
    """Translate a control-file keyword into the Fortran integer it mirrors."""
    key = value.upper()
    if key not in mapping:
        raise GenesisValidationError(
            f"Invalid {label}: {value}. Valid values: {list(mapping)}"
        )
    return mapping[key]


def _prepare_lazy_trajectory(trajs, ana_period: int):
    """Validate lazy metadata and prepare common C-call inputs."""
    if ana_period <= 0:
        raise GenesisValidationError(
            f"ana_period must be positive, got {ana_period}"
        )
    if not trajs.is_lazy or not trajs.lazy_dcd_file:
        raise GenesisValidationError("A lazy DCD trajectory is required")
    if not os.path.exists(trajs.lazy_dcd_file):
        raise GenesisValidationError(
            f"DCD file not found: {trajs.lazy_dcd_file}"
        )

    selection_indices = np.ascontiguousarray(
        trajs.selection_indices, dtype=np.int32
    )
    if selection_indices.ndim != 1 or selection_indices.size != trajs.natom:
        raise GenesisValidationError(
            "Lazy trajectory selection metadata is inconsistent"
        )

    result_frames = trajs.nframe // ana_period
    if result_frames <= 0:
        raise GenesisValidationError(
            f"ana_period={ana_period} selects no trajectory frames"
        )

    return (
        trajs.lazy_dcd_file.encode("utf-8"),
        selection_indices,
        trajs.lazy_ana_period * ana_period,
        result_frames,
    )


def _pack_filenames(filenames: List[str]) -> Tuple[bytes, int, int]:
    """Pack list of filenames into fixed-width byte buffer for Fortran.

    Args:
        filenames: List of file path strings

    Returns:
        Tuple of (packed_bytes, n_files, max_filename_len)
    """
    if not filenames:
        return b'', 0, 0
    max_len = max(len(f) for f in filenames)
    # Pad each filename to max_len with null bytes
    packed = b''.join(f.encode('utf-8').ljust(max_len, b'\x00') for f in filenames)
    return packed, len(filenames), max_len


def extract_model_blocks(pdb_string):
    """
    与えられたPDB形式の文字列から、MODEL行からENDMDL行の間のスライスを順に返す。
    文字列のコピーは行わない。

    Parameters:
        pdb_string (str): PDB形式の文字列。

    Yields:
        str: 各MODELブロック（元の文字列のスライス）。
    """
    start = None  # MODEL行の開始位置
    end = None    # ENDMDL行の終了位置

    i = 0
    while i < len(pdb_string):
        # MODEL行の開始位置を探す
        if pdb_string.startswith("MODEL", i):
            start = i
            # MODEL行の終わりまで進める
            while i < len(pdb_string) and pdb_string[i] != '\n':
                i += 1
            i += 1  # 改行をスキップ

        # ENDMDL行の終了位置を探す
        elif pdb_string.startswith("ENDMDL", i):
            end = i
            while i < len(pdb_string) and pdb_string[i] != '\n':
                i += 1
            i += 1  # 改行をスキップ
            yield pdb_string[start:i]  # MODEL行からENDMDL行までのスライスを返す
            start = None
            end = None

        else:
            i += 1  # 次の文字へ進む


def _import_root() -> str:
    """Return the directory the ``genepie`` package can be imported from.

    Derived from the package itself rather than from this module's path so that
    it stays correct regardless of how deeply nested this module is.
    """
    from .. import __file__ as genepie_init
    return os.path.dirname(os.path.dirname(os.path.abspath(genepie_init)))


def run_analysis_isolated(
    func_name: str,
    timeout: Optional[float],
    task_description: str,
    **kwargs,
):
    """Run ``genesis_exe.<func_name>(**kwargs)`` in a throwaway subprocess.

    Some GENESIS solvers keep global Fortran state (iteration workspaces, file
    units, allocations) alive between calls. For the free energy estimators this
    accumulation crashes the interpreter once four or more solves run in the same
    process. Executing each solve in its own subprocess hands every call clean
    state, so any number of estimates can be computed back to back in one kernel
    session.

    The child returns the wrapper's result object unchanged (pickled), so this
    works for wrappers that return a plain ``numpy`` array as well as those that
    return a ``NamedTuple``.

    Args:
        func_name: Name of the ``genesis_exe`` function to call.
        timeout: Maximum time in seconds to wait for completion (None = no limit).
        task_description: Short label used in error/timeout messages.
        **kwargs: Arguments forwarded verbatim to the wrapper.

    Returns:
        Whatever the wrapped function returns.

    Raises:
        The same exception types the wrapped function would raise, plus
        ``TimeoutError`` if the solve exceeds ``timeout`` and ``RuntimeError`` if
        the subprocess fails to produce a decodable result.
    """
    import base64
    import pickle
    import subprocess
    import sys

    kwargs_bytes = base64.b64encode(pickle.dumps(kwargs)).decode('ascii')

    script = f'''
import sys
import pickle
import base64

try:
    from genepie import genesis_exe
except ImportError:
    sys.path.insert(0, "{_import_root()}")
    from genepie import genesis_exe

kwargs = pickle.loads(base64.b64decode("{kwargs_bytes}"))

try:
    result = genesis_exe.{func_name}(**kwargs)
    output = {{"success": True, "result": result}}
except Exception as e:
    import traceback
    output = {{
        "success": False,
        "error": str(e),
        "error_type": type(e).__name__,
        "traceback": traceback.format_exc(),
    }}

sys.stdout.buffer.write(base64.b64encode(pickle.dumps(output)))
'''

    try:
        proc = subprocess.run(
            [sys.executable, '-c', script],
            capture_output=True,
            timeout=timeout,
            env={**os.environ,
                 'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS', '1')},
        )
    except subprocess.TimeoutExpired as e:
        raise TimeoutError(
            f"{task_description} timed out after {timeout} seconds"
        ) from e

    if proc.returncode != 0:
        stderr_text = proc.stderr.decode('utf-8', errors='replace')
        raise RuntimeError(
            f"{task_description} subprocess failed with code "
            f"{proc.returncode}:\n{stderr_text}"
        )

    try:
        output = pickle.loads(base64.b64decode(proc.stdout))
    except Exception as e:
        raise RuntimeError(
            f"Failed to decode {task_description} subprocess output: {e}\n"
            f"stdout: {proc.stdout[:500]!r}\n"
            f"stderr: {proc.stderr.decode('utf-8', errors='replace')}"
        )

    if not output["success"]:
        _reraise_subprocess_error(output)

    return output["result"]


def _reraise_subprocess_error(output: dict) -> None:
    """Rebuild the child's exception in the parent process."""
    from .. import exceptions as _exc

    error_type = output.get("error_type", "")
    error_msg = output.get("error", "Unknown error")
    traceback_str = output.get("traceback", "")

    exc_cls = getattr(_exc, error_type, None)
    if isinstance(exc_cls, type) and issubclass(exc_cls, Exception):
        try:
            raise exc_cls(error_msg)
        except TypeError:
            # Exception subclass with a non-standard __init__ signature.
            raise _exc.GenesisError(f"{error_type}: {error_msg}")

    raise RuntimeError(f"{error_type}: {error_msg}\n{traceback_str}")
