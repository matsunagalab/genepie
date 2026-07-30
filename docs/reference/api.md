# Analysis API Reference

All analysis functions live in the `genepie.analysis` package and are re-exported
through `genepie.genesis_exe` for the historical flat API. Both import styles work:

```python
from genepie import genesis_exe, SMolecule
result = genesis_exe.rmsd_analysis(...)

# equivalent, using the split package directly
from genepie.analysis import rmsd_analysis, crd_convert
```

In the tables below:

- **Zerocopy** — results/coordinates share memory with the NumPy array (no copy).
- **Lazy** — accepts a selected, strided lazy trajectory from
  `crd_convert(..., lazy=True)` and streams frames from disk. The dedicated
  `rmsd_analysis_lazy` remains as a compatibility wrapper.

## Trajectory loading & selection

| Function | Purpose | Zerocopy | Lazy |
|----------|---------|:--------:|:----:|
| `crd_convert(mol, trj_files, ...)` | Read/convert a trajectory. Returns `(list[STrajectories], SMolecule)`. Supports `selection`, `fitting_selection`/`fitting_method`, `centering`, `pbc_correct`, `ana_period`, `rename_res`, and `lazy=True`. | ✓ | ✓ |
| `crd_convert_info(mol, trj_files, ...)` | Read only the frame counts (used internally to pre-allocate). | – | – |
| `selection(mol, "an:CA")` | Evaluate a GENESIS selection string; returns 1-indexed atom indices. | – | – |

In lazy mode, `selection` and `ana_period` define the same logical trajectory
view as eager loading: the returned trajectory and subset molecule have matching
atom counts, and an analysis-level `ana_period` applies a second stride to that
view. Lazy conversion supports one DCD file and does not perform conversion-time
fitting, centering, or PBC coordinate correction. RMSD fitting and DRMS PBC
handling remain available at analysis time.

```{admonition} crd_convert returns a tuple
:class: important
`trajs, subset_mol = crd_convert(...)`. `trajs` is a list (one entry per input
file); `subset_mol` is the molecule reduced to the selected atoms.
```

## Structural analyses

| Function | Returns (namedtuple) | Zerocopy | Lazy |
|----------|----------------------|:--------:|:----:|
| `rmsd_analysis(mol, traj, analysis_selection, [fitting_selection, fitting_method], ...)` | `RmsdAnalysisResult(rmsd)` | ✓ | ✓ |
| `rmsd_analysis_lazy(mol, dcd_file, analysis_selection, ..., has_box, max_frames)` | `RmsdLazyAnalysisResult(rmsd, dcd_nframe, dcd_natom)` | ✓ | ✓ |
| `rg_analysis(mol, traj, analysis_selection, mass_weighted, ...)` | `RgAnalysisResult(rg)` | ✓ | ✓ (via lazy traj) |
| `drms_analysis(traj, contact_list, contact_dist, ...)` | `DrmsAnalysisResult(drms)` | ✓ | ✓ (via lazy traj) |
| `trj_analysis(traj, [distance_pairs, angle_triplets, torsion_quadruplets, cdis_groups, cang_groups], ...)` | `TrjAnalysisResult(distance, angle, torsion, com_distance, com_angle)` | ✓ | – |
| `avecrd_analysis(mol, traj, selection_group, fitting_method, ...)` | `AvecrdAnalysisResult(pdb, ...)` | – | – |
| `hb_analysis(mol, traj, selection_group, output_type, ...)` | hydrogen-bond counts/records | – | – |

Selection strings and `fitting_method` values (`"NO"`, `"TR"`, `"TR+ROT"`,
`"TR+ZROT"`, `"XYTR"`, `"XYTR+ZROT"`) are **case-insensitive**.

## Dynamics

| Function | Returns | Zerocopy | Lazy |
|----------|---------|:--------:|:----:|
| `msd_analysis(mol, traj, selection_group, selection, mode, ...)` | `MsdAnalysisResult(msd, ...)` | – | – |
| `diffusion_analysis(msd_array, time_step, start_step)` | `DiffusionAnalysisResult(out_data, diffusion_coefficients)` | ✓ | – |

## Free energy

| Function | Returns | Notes |
|----------|---------|-------|
| `wham_analysis(cvfile, dimension, grids, constant, reference, ...)` | PMF array `(n_bins, n_columns)` | `cvfile` **required**; existence validated up front (raises `GenesisValidationError`). |
| `mbar_analysis(cvfile, nreplica, input_type, grids, constant, reference, ...)` | free-energy array | Same `cvfile` validation. DCD input is rejected with `GenesisFortranNotSupportedError`. |
| `pmf_analysis(cv/cvfile, weight/weightfile, grids, band_width, ...)` | `Pmf1DResult(cv, pmf, pmf_gaussian)` or `Pmf2DResult(cv1, cv2, pmf)` | Builds a PMF directly from CV samples and optional per-sample weights (histogram + Gaussian kernel, `F = -kT log P`). Accepts in-memory arrays (`cv`, `weight`) or CLI-style filename patterns (`cvfile`, `weightfile`). |

`cvfile` is a filename *pattern* where `{}` expands to the window/replica index,
e.g. `"run{}.dis"`.

`wham_analysis`, `mbar_analysis`, and `pmf_analysis` each have an
`*_isolated(..., timeout=...)` variant that runs the solve in a throwaway
subprocess, giving every call clean Fortran state so many estimates can be
computed back to back in one session.

## Clustering

| Function | Returns | Notes |
|----------|---------|-------|
| `kmeans_clustering(mol, traj, selection_group, num_clusters, ...)` | `KmeansClusteringResult(cluster_idxs, mols_from_pdb, ...)` | Writes representative-structure PDBs to the working directory. |

## MD engine (ATDYN)

| Function | Returns | Notes |
|----------|---------|-------|
| `run_atdyn_md(...)` | `AtdynMDResult(energies, final_coords, energy_labels)` | `energies` is `(n_terms, n_records)`. |
| `run_atdyn_min(...)` | `AtdynMinResult(energies, final_coords, converged, final_gradient, energy_labels)` | Steepest descent / other minimizers. |
| `run_atdyn_md_isolated(...)` | same as `run_atdyn_md` | Runs in a subprocess (crash-safe). |
| `run_atdyn_min_isolated(...)` | same as `run_atdyn_min` | Runs in a subprocess (crash-safe). |

Supported inputs: **AMBER** (`prmtopfile`/`ambcrdfile`), **GROMACS**
(`grotopfile`/`grocrdfile`), **CHARMM** (`psffile`/`pdbfile` + `parfile`/`strfile`).

```{admonition} Sequential runs
:class: note
The Fortran engine accumulates global state (FFT plans, timers) between runs.
Two or three sequential in-process runs are generally fine; for many runs, prefer
the `*_isolated` variants.
```

See the [error model](errors.md) for the exception hierarchy every function uses.
