"""remd_convert matches the CLI sort for coordinates."""

import math
import struct
import subprocess

import numpy as np
import pytest

from .conftest import REPO_ROOT, TREMD_PDB, TREMD_PSF
from ..exceptions import GenesisError, GenesisValidationError
from ..s_molecule import SMolecule
from .. import genesis_exe

BINARY = REPO_ROOT / "src" / "analysis" / "converter" / "remd_convert" / "remd_convert"
NREP = 4
NFRAME = 16
CRD = 50
EXC = 100
NSTEPS = NFRAME * CRD


pytestmark = pytest.mark.skipif(
    not (TREMD_PDB.is_file() and TREMD_PSF.is_file() and BINARY.is_file()),
    reason="trialanine topology or remd_convert binary is absent",
)


def _record(payload: bytes) -> bytes:
    n = len(payload)
    return struct.pack("<i", n) + payload + struct.pack("<i", n)


def _write_dcd(path, natom: int, nframe: int, replica: int) -> None:
    icntrl = [0] * 20
    icntrl[0] = nframe
    icntrl[1] = 1
    icntrl[2] = 1
    icntrl[10] = 1
    icntrl[19] = 24
    atoms = np.arange(natom, dtype=np.float32)
    with path.open("wb") as handle:
        handle.write(_record(b"CORD" + struct.pack("<20i", *icntrl)))
        handle.write(_record(struct.pack("<i", 1) + b"REMARKS TEST".ljust(80)))
        handle.write(_record(struct.pack("<i", natom)))
        for frame in range(1, nframe + 1):
            lx = 20.0 + 0.1 * frame + replica
            ly = 21.0 + 0.2 * frame
            lz = 22.0 + 0.3 * frame
            handle.write(_record(struct.pack("<6d", lx, 0.0, ly, 0.0, 0.0, lz)))
            base = np.float32(replica * 100 + frame)
            handle.write(_record((base + atoms * np.float32(0.01)).tobytes()))
            handle.write(_record((base * np.float32(0.5) + atoms).tobytes()))
            handle.write(_record((np.float32(frame) - atoms * np.float32(0.02)).tobytes()))


def _write_inputs(tmp, natom: int):
    for replica in range(1, NREP + 1):
        _write_dcd(tmp / f"rep{replica}.dcd", natom, NFRAME, replica)
        lines = []
        nrem = NSTEPS // EXC
        for k in range(nrem):
            parmid = ((replica - 1 + k) % NREP) + 1
            lines.append(f"{k * EXC:10d}{parmid:10d}")
        (tmp / f"rep{replica}.rem").write_text("\n".join(lines) + "\n")


def _cli(tmp, out, selection="all", ids="", convert_type="PARAMETER",
         with_rem=True):
    out.mkdir()
    rem = f"remfile = {tmp}/rep{{}}.rem\n" if with_rem else ""
    text = f"""[INPUT]
psffile = {TREMD_PSF}
reffile = {TREMD_PDB}
{rem}dcdfile = {tmp}/rep{{}}.dcd

[OUTPUT]
trjfile = {out}/param{{}}.dcd

[SELECTION]
group1 = {selection}

[FITTING]
fitting_method = NO

[OPTION]
convert_type = {convert_type}
num_replicas = {NREP}
convert_ids = {ids}
nsteps = {NSTEPS}
exchange_period = {EXC}
crdout_period = {CRD}
trjout_format = DCD
trjout_type = COOR+BOX
trjout_atom = 1
centering = NO
pbc_correct = NO
"""
    ctrl = tmp / "remd_convert.inp"
    ctrl.write_text(text)
    proc = subprocess.run(
        [str(BINARY), str(ctrl)], cwd=out,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert proc.returncode == 0, proc.stdout[-2000:]


def _load_dcd(path):
    with path.open("rb") as handle:
        def rec():
            raw = handle.read(4)
            if len(raw) < 4:
                return None
            n = struct.unpack("<i", raw)[0]
            data = handle.read(n)
            end = handle.read(4)
            assert len(data) == n and struct.unpack("<i", end)[0] == n
            return data

        header = rec()
        has_box = struct.unpack("<20i", header[4:84])[10] != 0
        rec()
        natom = struct.unpack("<i", rec())[0]
        frames = []
        while True:
            pos = handle.tell()
            if len(handle.read(4)) < 4:
                break
            handle.seek(pos)
            box = None
            if has_box:
                sym = struct.unpack("<6d", rec())
                box = (
                    math.sqrt(sym[0] ** 2 + sym[1] ** 2 + sym[3] ** 2),
                    math.sqrt(sym[1] ** 2 + sym[2] ** 2 + sym[4] ** 2),
                    math.sqrt(sym[3] ** 2 + sym[4] ** 2 + sym[5] ** 2),
                )
            xyz = np.stack([
                np.frombuffer(rec(), dtype="<f4"),
                np.frombuffer(rec(), dtype="<f4"),
                np.frombuffer(rec(), dtype="<f4"),
            ], axis=1).astype(np.float64)
            frames.append((box, xyz.copy()))
    return natom, frames


def _compare(out, result):
    names = sorted(out.glob("param*.dcd"))
    assert [int(p.stem.replace("param", "")) for p in names] == result.param_ids
    for path, traj in zip(names, result.trajectories):
        natom, frames = _load_dcd(path)
        assert natom == traj.coords.shape[1]
        assert len(frames) == traj.coords.shape[0]
        assert traj.coords.shape[1] == result.subset.num_atoms
        for i, (box, xyz) in enumerate(frames):
            np.testing.assert_allclose(traj.coords[i], xyz, rtol=0, atol=0)
            diag = np.array([
                traj.pbc_boxes[i, 0, 0],
                traj.pbc_boxes[i, 1, 1],
                traj.pbc_boxes[i, 2, 2],
            ])
            np.testing.assert_allclose(diag, box, rtol=0, atol=1e-10)


def test_remd_convert_matches_cli_parameter(tmp_path):
    mol = SMolecule.from_file(pdb=TREMD_PDB, psf=TREMD_PSF)
    _write_inputs(tmp_path, mol.num_atoms)
    out = tmp_path / "cli"
    _cli(tmp_path, out)
    result = genesis_exe.remd_convert(
        mol,
        dcd_files=[str(tmp_path / f"rep{i}.dcd") for i in range(1, NREP + 1)],
        rem_files=[str(tmp_path / f"rep{i}.rem") for i in range(1, NREP + 1)],
        convert_type="PARAMETER",
        nsteps=NSTEPS,
        exchange_period=EXC,
        crdout_period=CRD,
    )
    assert result.param_ids == [1, 2, 3, 4]
    _compare(out, result)


def test_remd_convert_matches_cli_selection_and_ids(tmp_path):
    mol = SMolecule.from_file(pdb=TREMD_PDB, psf=TREMD_PSF)
    _write_inputs(tmp_path, mol.num_atoms)
    out = tmp_path / "cli"
    _cli(tmp_path, out, selection="an:CA", ids="1 3")
    result = genesis_exe.remd_convert(
        mol,
        dcd_files=[str(tmp_path / f"rep{i}.dcd") for i in range(1, NREP + 1)],
        rem_files=[str(tmp_path / f"rep{i}.rem") for i in range(1, NREP + 1)],
        nsteps=NSTEPS,
        exchange_period=EXC,
        crdout_period=CRD,
        convert_ids=[1, 3],
        selection="an:CA",
    )
    assert result.param_ids == [1, 3]
    assert result.subset.num_atoms < mol.num_atoms
    _compare(out, result)


def test_remd_convert_matches_cli_replica(tmp_path):
    mol = SMolecule.from_file(pdb=TREMD_PDB, psf=TREMD_PSF)
    _write_inputs(tmp_path, mol.num_atoms)
    out = tmp_path / "cli"
    _cli(tmp_path, out, convert_type="REPLICA", with_rem=False)
    result = genesis_exe.remd_convert(
        mol,
        dcd_files=[str(tmp_path / f"rep{i}.dcd") for i in range(1, NREP + 1)],
        convert_type="REPLICA",
        nsteps=NSTEPS,
        exchange_period=EXC,
        crdout_period=CRD,
    )
    _compare(out, result)


def test_remd_convert_rejects_bad_periods():
    mol = SMolecule.from_file(pdb=TREMD_PDB, psf=TREMD_PSF)
    with pytest.raises(GenesisValidationError):
        genesis_exe.remd_convert(
            mol, dcd_files=["a.dcd"], rem_files=["a.rem"],
            nsteps=0, exchange_period=1, crdout_period=1)


def test_remd_convert_missing_dcd_is_an_error(tmp_path):
    mol = SMolecule.from_file(pdb=TREMD_PDB, psf=TREMD_PSF)
    _write_inputs(tmp_path, mol.num_atoms)
    dcds = [str(tmp_path / f"rep{i}.dcd") for i in range(1, NREP + 1)]
    rems = [str(tmp_path / f"rep{i}.rem") for i in range(1, NREP + 1)]
    dcds[0] = str(tmp_path / "missing.dcd")
    with pytest.raises(GenesisError):
        genesis_exe.remd_convert(
            mol, dcd_files=dcds, rem_files=rems,
            nsteps=NSTEPS, exchange_period=EXC, crdout_period=CRD)
