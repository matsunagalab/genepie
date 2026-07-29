"""Helpers shared by more than one analysis wrapper."""

from typing import (
    List,
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
