# SPDX-License-Identifier: LGPL-3.0-or-later
"""
CLI entry points for GENESIS commands.

This module provides Python entry points for bundled GENESIS binaries.
"""

import subprocess
import sys
from pathlib import Path


def get_bin_path(name: str) -> Path:
    """Get path to bundled binary.

    Args:
        name: Name of the binary (e.g., 'atdyn')

    Returns:
        Path to the binary executable
    """
    return Path(__file__).parent / "bin" / name


def _run_binary(name: str):
    """Run a bundled binary with command line arguments."""
    binary = get_bin_path(name)
    if not binary.exists():
        print(f"Error: {name} binary not found at {binary}", file=sys.stderr)
        print("The binary may not be included in this installation.", file=sys.stderr)
        print("Please build GENESIS from source or install via conda.", file=sys.stderr)
        sys.exit(1)
    sys.exit(subprocess.call([str(binary)] + sys.argv[1:]))


# MD Engine
def run_atdyn():
    """Run atdyn molecular dynamics engine."""
    _run_binary("atdyn")


# Note: spdyn requires MPI and is not included in PyPI package.
# Users who need spdyn should build GENESIS from source or use conda.


# Clustering
def run_kmeans_clustering():
    """Run k-means clustering analysis."""
    _run_binary("kmeans_clustering")


# Converters
def run_cg_convert():
    """Run CG model converter."""
    _run_binary("cg_convert")


def run_crd_convert():
    """Run coordinate/trajectory converter."""
    _run_binary("crd_convert")


def run_pcrd_convert():
    """Run parallel coordinate converter."""
    _run_binary("pcrd_convert")


def run_remd_convert():
    """Run REMD trajectory converter."""
    _run_binary("remd_convert")


def run_rst_convert():
    """Run restart file converter."""
    _run_binary("rst_convert")


def run_rst_upgrade():
    """Run restart file upgrade tool."""
    _run_binary("rst_upgrade")


# Free Energy
def run_mbar_analysis():
    """Run MBAR free energy analysis."""
    _run_binary("mbar_analysis")


def run_meanforce_analysis():
    """Run mean force analysis."""
    _run_binary("meanforce_analysis")


def run_pmf_analysis():
    """Run PMF analysis."""
    _run_binary("pmf_analysis")


def run_wham_analysis():
    """Run WHAM free energy analysis."""
    _run_binary("wham_analysis")


# Interface
def run_dssp_interface():
    """Run DSSP interface."""
    _run_binary("dssp_interface")


# Mode Analysis
def run_avecrd_analysis():
    """Run average coordinate analysis."""
    _run_binary("avecrd_analysis")


def run_eigmat_analysis():
    """Run eigenvalue matrix analysis."""
    _run_binary("eigmat_analysis")


def run_flccrd_analysis():
    """Run fluctuation coordinate analysis."""
    _run_binary("flccrd_analysis")


def run_pcavec_drawer():
    """Run PCA vector drawer."""
    _run_binary("pcavec_drawer")


def run_prjcrd_analysis():
    """Run projection coordinate analysis."""
    _run_binary("prjcrd_analysis")


# SP Analysis
def run_contact_analysis():
    """Run contact analysis."""
    _run_binary("contact_analysis")


def run_density_analysis():
    """Run density analysis."""
    _run_binary("density_analysis")


def run_hbond_analysis():
    """Run hydrogen bond analysis (sp_analysis)."""
    _run_binary("hbond_analysis")


def run_rdf_analysis():
    """Run radial distribution function analysis."""
    _run_binary("rdf_analysis")


def run_sasa_analysis():
    """Run solvent accessible surface area analysis."""
    _run_binary("sasa_analysis")


# Trajectory Analysis
def run_comcrd_analysis():
    """Run center of mass coordinate analysis."""
    _run_binary("comcrd_analysis")


def run_diffusion_analysis():
    """Run diffusion coefficient analysis."""
    _run_binary("diffusion_analysis")


def run_distmat_analysis():
    """Run distance matrix analysis."""
    _run_binary("distmat_analysis")


def run_drms_analysis():
    """Run distance RMSD analysis."""
    _run_binary("drms_analysis")


def run_energy_analysis():
    """Run energy analysis."""
    _run_binary("energy_analysis")


def run_fret_analysis():
    """Run FRET analysis."""
    _run_binary("fret_analysis")


def run_hb_analysis():
    """Run hydrogen bond analysis (trj_analysis)."""
    _run_binary("hb_analysis")


def run_lipidthick_analysis():
    """Run lipid thickness analysis."""
    _run_binary("lipidthick_analysis")


def run_msd_analysis():
    """Run mean squared displacement analysis."""
    _run_binary("msd_analysis")


def run_qval_analysis():
    """Run Q-value analysis."""
    _run_binary("qval_analysis")


def run_qval_residcg_analysis():
    """Run Q-value residue CG analysis."""
    _run_binary("qval_residcg_analysis")


def run_rg_analysis():
    """Run radius of gyration analysis."""
    _run_binary("rg_analysis")


def run_ring_analysis():
    """Run ring analysis."""
    _run_binary("ring_analysis")


def run_rmsd_analysis():
    """Run RMSD analysis."""
    _run_binary("rmsd_analysis")


def run_tilt_analysis():
    """Run tilt analysis."""
    _run_binary("tilt_analysis")


def run_trj_analysis():
    """Run trajectory analysis (distance, angle, dihedral)."""
    _run_binary("trj_analysis")


# Utilities
def run_emmap_generator():
    """Run EM map generator."""
    _run_binary("emmap_generator")


def run_morph_generator():
    """Run morph generator."""
    _run_binary("morph_generator")


def run_pathcv_analysis():
    """Run path collective variable analysis."""
    _run_binary("pathcv_analysis")


def run_qmmm_generator():
    """Run QM/MM input generator."""
    _run_binary("qmmm_generator")


def run_rpath_generator():
    """Run reaction path generator."""
    _run_binary("rpath_generator")
