#!/usr/bin/env python
"""Package the T-REMD (trialanine / Ala3) tutorial dataset into a self-contained zip.

The full ``remd_alat`` directory from the MBAR paper is ~14 GB (mostly the raw
solvated production DCDs, 712 MB per replica). For the genepie MBAR reweighting /
resampling tutorial we only need the parameter-sorted, solute-only subset plus the
reference MBAR results. This script collects those files (~60-70 MB), writes a
self-describing README, and produces ``remd_alat_tutorial.zip`` ready to upload to
Google Drive (downloaded on demand via gdown, same pattern as the chignolin data).

Usage:
    python scripts/pack_remd_alat_tutorial.py \
        [--src /Users/yasu/gdrive/paper/paper_mbar/data/remd_alat] \
        [--out ./remd_alat_tutorial.zip]

The produced zip is NOT meant to be committed to the repository.
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib
import shutil
import tempfile
import zipfile

# Temperature ladder of the 20 T-REMD replicas (K). Index k (1-based) is the
# temperature of parameter-sorted state remd_paramID{k}.
TEMPERATURES = [
    300.00, 302.53, 305.09, 307.65, 310.24,
    312.85, 315.47, 318.12, 320.78, 323.46,
    326.16, 328.87, 331.61, 334.37, 337.14,
    339.94, 342.75, 345.59, 348.44, 351.26,
]
NREPLICA = len(TEMPERATURES)
TARGET_TEMPERATURE = 300.0
NFRAME_PER_STATE = 5000

DEFAULT_SRC = "/Users/yasu/gdrive/paper/paper_mbar/data/remd_alat"


def _readme_text() -> str:
    temp_lines = "\n".join(
        f"  {k:2d}   remd_paramID{k}   {t:.2f}"
        for k, t in enumerate(TEMPERATURES, start=1)
    )
    return f"""# T-REMD trialanine (Ala3) dataset - MBAR reweighting & resampling tutorial

This is a self-contained subset of a Temperature Replica-Exchange Molecular
Dynamics (T-REMD) simulation of **trialanine (Ala3, 42 atoms)**, prepared for the
genepie tutorial that goes from a standard MBAR free-energy analysis to
weight-based trajectory resampling.

It is derived from the input/output files of:

> Matsunaga et al., "Use of multistate Bennett acceptance ratio method for
> free-energy calculations from enhanced sampling and free energy perturbation".
> Original files: https://github.com/matsunagalab/paper_mbar (directory `remd_alat`).

## System and simulation

- Solute: trialanine (Ala3), acetyl/N-methyl capped, CHARMM36 force field, 42 atoms.
- The MD was run with explicit TIP3 water under PBC + PME, NVT (Langevin), GENESIS.
- Enhanced sampling: **Temperature REMD**, {NREPLICA} replicas.
- Production: 2,500,000 steps, timestep 0.002 ps, crdout/eneout every 500 steps
  -> **{NFRAME_PER_STATE} frames per state**.
- Trajectories/energies were sorted **by parameter (temperature)** with GENESIS
  `remd_convert` (`convert_type = PARAMETER`). So `remd_paramID{{k}}` is the
  ensemble at a fixed temperature (see the ladder below), NOT a single physical
  replica walker.
- The DCDs here contain **only the 42-atom solute** (water was dropped during
  `remd_convert`); the `.pot` energies are the **full-system** potential energy
  (this is what determines the Boltzmann/MBAR weights).

## Temperature ladder (state -> temperature)

```
  k    state          T (K)
{temp_lines}
```

Target temperature for reweighting/resampling: **{TARGET_TEMPERATURE:.1f} K**
(= state 1, remd_paramID1).

## Files

Primary tutorial inputs (flat, in the root):

- `trialanine.pdb`          - reference coordinates / topology, 42 atoms.
- `ala3.psf`                - CHARMM PSF (topology) for the 42-atom solute.
- `remd_paramID{{1..{NREPLICA}}}_trialanine.dcd`
                            - parameter-sorted solute coordinates,
                              {NFRAME_PER_STATE} frames each.
- `remd_paramID{{1..{NREPLICA}}}.pot`
                            - full-system potential energy time series,
                              {NFRAME_PER_STATE} lines of `time  PE[kcal/mol]`.
                              This is the MBAR input (input_type = EneSingle).

Reference results (in `reference/`) computed in the original paper - use these to
validate the genepie implementation:

- `weight{{1..{NREPLICA}}}.dat` - per-frame MBAR weights at {TARGET_TEMPERATURE:.0f} K,
                              {NFRAME_PER_STATE} lines of `time  weight`.
- `fene.dat`                - relative free energies of the {NREPLICA} states
                              ({NREPLICA} lines).
- `remd_paramID{{1..{NREPLICA}}}.tor` - precomputed backbone dihedrals,
                              {NFRAME_PER_STATE} lines of `index  phi  psi` (deg),
                              for Ramachandran plots.

Provenance (in `control_files/`) - original GENESIS control files:

- `production_run.inp`      - the T-REMD production run.
- `remd_convert.inp`        - parameter sorting (raw -> remd_paramID*).
- `mbar_analysis.inp`       - the MBAR analysis that produced weight*/fene.

`temperatures.txt` lists the ladder in plain text.

## How the tutorial uses it

1. Standard MBAR: feed `remd_paramID{{}}.pot` with
   `input_type = EneSingle`, `temperature = [ladder]`, `target_temperature = 300`,
   `nreplica = {NREPLICA}` -> relative free energies `fene` (compare to `fene.dat`).
2. Weights: `return_weights=True` -> per-frame weights at 300 K
   (compare to `weight{{}}.dat`).
3. Reweighted Ramachandran free-energy surface from (phi, psi) + weights.
4. Resampling with replacement using the weights -> a uniform-weight 300 K
   trajectory whose plain histogram reproduces the reweighted surface.

## Notes / caveats

- `.pot` line count (energy samples) equals DCD frame count (both {NFRAME_PER_STATE});
  they are aligned 1:1, which the resampling relies on.
- Weights/free energies were computed with `nblocks = 1`.

## License

Input files taken from GENESIS tutorials retain their original terms; the rest of
the paper_mbar material is GPL-3.0. Please cite the paper above and acknowledge the
GENESIS tutorials (https://www.r-ccs.riken.jp/labs/cbrt/) when using this data.
"""


def _copy(src: pathlib.Path, dst: pathlib.Path) -> int:
    if not src.exists():
        raise FileNotFoundError(f"missing expected source file: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst.stat().st_size


def build(src_root: pathlib.Path, out_zip: pathlib.Path) -> None:
    sort_dir = src_root / "5_analysis" / "4_sort_dcd"
    mbar_dir = src_root / "5_analysis" / "mbar"
    setup_dir = src_root / "1_setup" / "build"

    with tempfile.TemporaryDirectory() as tmp:
        stage = pathlib.Path(tmp) / "remd_alat_tutorial"
        stage.mkdir(parents=True)

        total = 0
        # Topology / reference structure.
        total += _copy(sort_dir / "trialanine.pdb", stage / "trialanine.pdb")
        total += _copy(setup_dir / "ala3.psf", stage / "ala3.psf")

        # Primary inputs: sorted solute DCDs + potential energies.
        for k in range(1, NREPLICA + 1):
            total += _copy(
                sort_dir / f"remd_paramID{k}_trialanine.dcd",
                stage / f"remd_paramID{k}_trialanine.dcd",
            )
            total += _copy(
                sort_dir / f"remd_paramID{k}.pot",
                stage / f"remd_paramID{k}.pot",
            )

        # Reference results.
        ref = stage / "reference"
        for k in range(1, NREPLICA + 1):
            total += _copy(mbar_dir / f"weight{k}.dat", ref / f"weight{k}.dat")
            total += _copy(
                sort_dir / f"remd_paramID{k}.tor", ref / f"remd_paramID{k}.tor"
            )
        total += _copy(mbar_dir / "fene.dat", ref / "fene.dat")

        # Provenance control files.
        ctrl = stage / "control_files"
        total += _copy(src_root / "4_production" / "run.inp", ctrl / "production_run.inp")
        total += _copy(
            sort_dir / "step5.remd_convert.inp", ctrl / "remd_convert.inp"
        )
        total += _copy(mbar_dir / "inp", ctrl / "mbar_analysis.inp")

        # Self-describing metadata.
        (stage / "README.md").write_text(_readme_text())
        (stage / "temperatures.txt").write_text(
            "\n".join(f"{k}\t{t:.2f}" for k, t in enumerate(TEMPERATURES, 1)) + "\n"
        )

        out_zip.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            for path in sorted(stage.rglob("*")):
                if path.is_file():
                    zf.write(path, path.relative_to(stage.parent))

    size = out_zip.stat().st_size
    sha = hashlib.sha256(out_zip.read_bytes()).hexdigest()
    print(f"staged uncompressed: {total / 1e6:.1f} MB")
    print(f"zip written:         {out_zip}  ({size / 1e6:.1f} MB)")
    print(f"sha256:              {sha}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--src", default=DEFAULT_SRC, help="source remd_alat directory")
    p.add_argument(
        "--out", default="remd_alat_tutorial.zip", help="output zip path"
    )
    args = p.parse_args()
    build(pathlib.Path(args.src).expanduser(), pathlib.Path(args.out).expanduser())


if __name__ == "__main__":
    main()
