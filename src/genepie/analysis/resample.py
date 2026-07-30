"""Weight-based trajectory resampling from MBAR.

MBAR turns multi-state (e.g. T-REMD) sampling into a single set of per-sample
weights at a chosen *target* temperature. Drawing samples with replacement
according to those weights yields a plain, uniform-weight trajectory that
approximates the target ensemble - something that can be treated like an
ordinary MD trajectory (histogrammed, deposited to a database, etc.).

This module provides:

* :func:`resample_indices` - the pure resampling step (weights -> frame indices),
  handy on its own and easy to test.
* :func:`mbar_resample_trajectory` - the full pipeline: MBAR weights ->
  resampled coordinates -> optional DCD.
"""

from collections import namedtuple
from typing import List, Optional, Sequence, Union

import numpy as np

from ..exceptions import GenesisValidationError
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from . import free_energy
from .converter import crd_convert


# Result of :func:`mbar_resample_trajectory`.
#   trajectory: STrajectories of the resampled frames (uniform weight).
#   molecule:   the (possibly subset) SMolecule matching ``trajectory``.
#   indices:    the drawn frame indices into the flattened (state, step) order.
#   weights:    the MBAR weights used, shape (n_replica, n_step).
MbarResampleResult = namedtuple(
    "MbarResampleResult", ["trajectory", "molecule", "indices", "weights"])


def resample_indices(
    weights,
    n_samples: Optional[int] = None,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Draw frame indices with replacement in proportion to ``weights``.

    Args:
        weights: MBAR weights. Either a 2-D ``(n_replica, n_step)`` array or an
            already-flattened 1-D array. 2-D input is flattened in C order, so
            index ``state * n_step + step`` refers to state ``state``'s
            ``step``-th sample - the same order
            :func:`mbar_resample_trajectory` concatenates coordinates in.
        n_samples: How many frames to draw. Defaults to the total number of
            input samples.
        seed: Seed for ``numpy.random.default_rng`` (pass an int for
            reproducibility).

    Returns:
        1-D ``int`` array of length ``n_samples`` with the chosen indices.

    Raises:
        GenesisValidationError: If the weights are empty, contain negatives, or
            do not sum to a positive finite value.
    """
    w = np.asarray(weights, dtype=np.float64).ravel()
    if w.size == 0:
        raise GenesisValidationError("resample_indices: weights are empty.")
    if np.any(w < 0) or not np.all(np.isfinite(w)):
        raise GenesisValidationError(
            "resample_indices: weights must be finite and non-negative."
        )
    total = w.sum()
    if not np.isfinite(total) or total <= 0.0:
        raise GenesisValidationError(
            f"resample_indices: weights sum to {total}, expected a positive "
            "value."
        )
    p = w / total

    if n_samples is None:
        n_samples = w.size
    if n_samples <= 0:
        raise GenesisValidationError(
            f"resample_indices: n_samples must be positive, got {n_samples}."
        )

    rng = np.random.default_rng(seed)
    return rng.choice(w.size, size=n_samples, replace=True, p=p)


def mbar_resample_trajectory(
    molecule: SMolecule,
    dcd_files: Sequence[str],
    cvfile: str,
    nreplica: int,
    temperature: Union[float, Sequence[float]],
    target_temperature: float,
    input_type: str = "EneSingle",
    dimension: int = 1,
    tolerance: float = 1.0e-8,
    self_iteration: Optional[int] = None,
    newton_iteration: Optional[int] = None,
    rest_function: Optional[Sequence[int]] = None,
    grids: Optional[Sequence] = None,
    constant: Optional[Sequence] = None,
    reference: Optional[Sequence] = None,
    is_periodic: Optional[Sequence] = None,
    box_size: Optional[Sequence] = None,
    selection_group: Optional[Sequence[str]] = None,
    selection_mole_name: Optional[Sequence[str]] = None,
    selection: str = "all",
    trj_format: str = "DCD",
    trj_type: str = "COOR+BOX",
    n_samples: Optional[int] = None,
    seed: Optional[int] = None,
    output_dcd: Optional[str] = None,
    isolated: bool = True,
) -> MbarResampleResult:
    """Resample a target-ensemble trajectory from multi-state data via MBAR.

    The MBAR weights returned by GENESIS are the *unbiased* ensemble weights at
    the requested target state, so this one function covers both canonical
    reweighting goals:

    * **Target temperature** (T-REMD): pass ``input_type="EneSingle"`` (or
      ``"REMD"``), a per-state ``temperature`` ladder and the desired
      ``target_temperature``. The weights approximate the unbiased ensemble at
      that single temperature.
    * **Bias removal** (umbrella sampling / restrained MD): pass
      ``input_type="US"`` (or ``"CV"``) together with the restraint definition
      (``rest_function``, ``grids``, ``constant``, ``reference``,
      ``is_periodic``, ``box_size``). GENESIS reweights every biased window to
      the restraint-free ensemble at ``target_temperature``; the weights are
      therefore the bias-free ensemble weights.

    The steps are:

    1. Run MBAR on ``cvfile`` to obtain per-sample weights for the target state.
    2. Load the coordinates of every state from ``dcd_files`` (optionally
       subset with ``selection``) and concatenate them in ``(state, step)``
       order, matching the weight layout.
    3. Draw ``n_samples`` frames with replacement in proportion to the weights.
    4. Optionally write the resampled frames to ``output_dcd``.

    ``dcd_files`` must be given in the same order the ``cvfile`` placeholder
    expands (state 1, 2, ... ``nreplica``), and each DCD must have exactly as
    many frames as its ``cvfile`` has samples; this 1:1 alignment is what makes
    the weights meaningful and is validated here.

    Args:
        molecule: Topology for the DCD files (full atom count).
        dcd_files: Per-state trajectory files, one per replica, in state order.
        cvfile: MBAR input pattern (e.g. ``".../remd_paramID{}.pot"`` for T-REMD
            or ``".../{}.dat"`` reaction-coordinate series for umbrella sampling).
        nreplica: Number of states/windows; must equal ``len(dcd_files)``.
        temperature: Per-state temperature(s) in K (scalar or one per state).
        target_temperature: Temperature the resampled ensemble approximates.
        input_type: MBAR input type; ``"EneSingle"``/``"REMD"`` for T-REMD,
            ``"US"``/``"CV"`` for umbrella sampling / biased MD.
        dimension: Number of reaction coordinates (umbrella case; 1 or 2).
        tolerance: MBAR convergence tolerance.
        self_iteration: MBAR self-consistent iteration count (optional).
        newton_iteration: MBAR Newton-Raphson iteration count (optional).
        rest_function: Restraint function indices used as reaction coordinates
            (umbrella case; forwarded to :func:`mbar_analysis`).
        grids: ``(min, max, num_grids)`` per dimension (umbrella case).
        constant: Restraint force constants per window (umbrella case).
        reference: Restraint centers per window (umbrella case).
        is_periodic: Whether each reaction coordinate is periodic (umbrella case).
        box_size: Period of each periodic reaction coordinate (umbrella case).
        selection_group: ``[SELECTION]`` group definitions (umbrella case).
        selection_mole_name: ``[SELECTION]`` mole-name definitions (umbrella case).
        selection: GENESIS atom selection applied when loading the DCDs
            (e.g. ``"all"`` when the DCDs already hold only the solute).
        trj_format: Trajectory format passed to :func:`crd_convert`.
        trj_type: ``"COOR+BOX"`` (default) or ``"COOR"`` passed to
            :func:`crd_convert`. Must match whether the DCDs carry box records.
        n_samples: Number of frames to draw (default: total sample count).
        seed: RNG seed for reproducible resampling.
        output_dcd: If given, the resampled trajectory is written here (mdtraj).
        isolated: Run the MBAR solve in a subprocess (recommended; avoids the
            Fortran global-state issues that crash repeated in-process solves).

    Returns:
        :class:`MbarResampleResult` ``(trajectory, molecule, indices, weights)``.

    Raises:
        GenesisValidationError: If ``dcd_files``/``nreplica`` disagree, if a
            DCD's frame count does not match its weight count, or on invalid
            weights.
    """
    dcd_files = list(dcd_files)
    if len(dcd_files) != nreplica:
        raise GenesisValidationError(
            f"mbar_resample_trajectory: got {len(dcd_files)} dcd_files but "
            f"nreplica={nreplica}. Provide one trajectory per state, in the "
            "same order the cvfile placeholder expands."
        )

    # Step 1: MBAR weights for the target state (target temperature and/or
    # bias-free). Umbrella/restraint arguments are forwarded verbatim so the US
    # and T-REMD paths share this single call.
    mbar = free_energy.mbar_analysis_isolated if isolated else free_energy.mbar_analysis
    result = mbar(
        cvfile=cvfile,
        nreplica=nreplica,
        input_type=input_type,
        dimension=dimension,
        temperature=temperature,
        target_temperature=target_temperature,
        tolerance=tolerance,
        self_iteration=self_iteration,
        newton_iteration=newton_iteration,
        rest_function=rest_function,
        grids=grids,
        constant=constant,
        reference=reference,
        is_periodic=is_periodic,
        box_size=box_size,
        selection_group=selection_group,
        selection_mole_name=selection_mole_name,
        nblocks=1,
        return_weights=True,
    )
    weights = np.asarray(result.weights, dtype=np.float64)
    n_step = weights.shape[1]

    # Step 2: load and concatenate per-state coordinates in (state, step) order.
    coords_states: List[np.ndarray] = []
    subset_mol: Optional[SMolecule] = None
    for state, dcd in enumerate(dcd_files, start=1):
        trajs, sm = crd_convert(
            molecule,
            [dcd],
            trj_format=trj_format,
            trj_type=trj_type,
            selection=selection,
        )
        traj = trajs[0]
        # Copy out of the zerocopy buffer so it survives the trajectory object.
        coords = np.array(traj.coords, dtype=np.float64)
        if coords.shape[0] != n_step:
            raise GenesisValidationError(
                f"mbar_resample_trajectory: state {state} DCD has "
                f"{coords.shape[0]} frames but its cvfile has {n_step} samples. "
                "The trajectory and energy series must be aligned 1:1."
            )
        coords_states.append(coords)
        subset_mol = sm

    coords_all = np.concatenate(coords_states, axis=0)

    # Step 3: resample with replacement according to the weights.
    indices = resample_indices(weights, n_samples=n_samples, seed=seed)
    resampled = np.ascontiguousarray(coords_all[indices])
    trajectory = STrajectories.from_numpy(resampled)

    # Step 4: optional DCD output.
    if output_dcd is not None:
        trajectory.save_dcd(output_dcd, subset_mol)

    return MbarResampleResult(
        trajectory=trajectory,
        molecule=subset_mol,
        indices=indices,
        weights=weights,
    )
