# SPDX-License-Identifier: LGPL-3.0-or-later
# Minimal, robust loader for libgenesis.* without LD_LIBRARY_PATH.
from __future__ import annotations

import ctypes
import os
from pathlib import Path
from importlib.resources import files
from typing import Iterable, Optional

CANDIDATES: tuple[str, ...] = (
        "libgenesis.so", 
        "libpython_interface.so", 
        "libgenesis.dylib", 
        "genesis.dll")


def _first_existing(paths: Iterable[str]) -> Optional[str]:
    for p in paths:
        if os.path.exists(p):
            return p
    return None

def _candidate_paths_in_dir(d: Path) -> List[Path]:
    return [d / name for name in CANDIDATES]

def _find_in_package_dotlib(pkg: str) -> Optional[Path]:
    try:
        libdir = Path(files(pkg).joinpath(".lib"))
    except Exception:
        return None
    return _first_existing(_candidate_paths_in_dir(libdir))


def _find_repo_root_lib_near(mod_file: str) -> Optional[Path]:
    here = Path(mod_file).resolve()
    for parent in list(here.parents)[:8]:
        cand_dir = parent / "lib"
        if cand_dir.is_dir():
            hit = _first_existing(_candidate_paths_in_dir(cand_dir))
            if hit:
                return hit
    return None

def _find_pkg_lib_python_interface_dotlib(mod_file: str) -> Optional[Path]:
    here = Path(mod_file).resolve()
    for parent in list(here.parents)[:8]:
        cand_dir = parent / "pkg" / "lib" / "python_interface" / ".lib"
        if cand_dir.is_dir():
            hit = _first_existing(_candidate_paths_in_dir(cand_dir))
            if hit:
                return hit
    return None

def _load(path: Path) -> ctypes.CDLL:
    return ctypes.CDLL(os.fspath(path), mode=ctypes.RTLD_GLOBAL)


def load_genesis_lib() -> ctypes.CDLL:

    tried: list[str] = []

    env_dir = os.environ.get("GENEPY_LIB_DIR")
    if env_dir:
        d = Path(env_dir)
        tried.append(f"{d}/(env)")
        hit = _first_existing(_candidate_paths_in_dir(d))
        if hit:
            return _load(hit)

    # in python_interface
    p = _find_in_package_dotlib("python_interface")
    tried.append("python_interface/.lib/")
    if p:
        return _load(p)

    # in genepy/.lib
    p = _find_in_package_dotlib("genepy")
    tried.append("genepy/.lib/")
    if p:
        return _load(p)

    # in <repo>/lib
    p = _find_repo_root_lib_near(__file__)
    tried.append("<repo>/lib/")
    if p:
        return _load(p)

    # in <repo>/pkg/lib/python_interface/.lib/
    p = _find_repo_pkg_python_interface_dotlib(__file__)
    tried.append("<repo>/pkg/lib/python_interface/.lib/")
    if p:
        return _load(p)

    raise OSError(
        "genesis shared library not found.\nSearched:\n  - "
        + "\n  - ".join(tried)
        + "\nHint: set GENEPY_LIB_DIR=/absolute/path/to/dir containing one of: "
        + ", ".join(CANDIDATES)
    )

