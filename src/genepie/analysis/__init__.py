"""Python wrappers around the GENESIS analysis tools and the ATDYN engine.

Each module here covers one tool family. :mod:`genepie.genesis_exe` re-exports
everything so that the historical flat API keeps working.
"""

from .atdyn import (
    AtdynMDResult,
    AtdynMinResult,
    run_atdyn_md,
    run_atdyn_md_isolated,
    run_atdyn_min,
    run_atdyn_min_isolated,
)
from .avecrd import AvecrdAnalysisResult, avecrd_analysis
from .clustering import KmeansClusteringResult, kmeans_clustering
from .converter import (
    CrdConvertInfo,
    crd_convert,
    crd_convert_info,
    selection,
    selection_func,
)
from .drms import DrmsAnalysisResult, drms_analysis
from .dynamics import (
    DiffusionAnalysisResult,
    MsdAnalysisResult,
    diffusion_analysis,
    msd_analysis,
)
from .free_energy import (
    mbar_analysis,
    mbar_analysis_isolated,
    wham_analysis,
    wham_analysis_isolated,
)
from .hbond import hb_analysis
from .rg import RgAnalysisResult, rg_analysis
from .rmsd import (
    RmsdAnalysisResult,
    RmsdLazyAnalysisResult,
    rmsd_analysis,
    rmsd_analysis_lazy,
)
from .trj import TrjAnalysisResult, trj_analysis
from ._common import extract_model_blocks

__all__ = [
    # Trajectory loading and conversion
    "CrdConvertInfo",
    "crd_convert",
    "crd_convert_info",
    "selection",
    "selection_func",
    # Structural analyses
    "TrjAnalysisResult",
    "trj_analysis",
    "RgAnalysisResult",
    "rg_analysis",
    "RmsdAnalysisResult",
    "RmsdLazyAnalysisResult",
    "rmsd_analysis",
    "rmsd_analysis_lazy",
    "DrmsAnalysisResult",
    "drms_analysis",
    "AvecrdAnalysisResult",
    "avecrd_analysis",
    "hb_analysis",
    # Dynamics
    "MsdAnalysisResult",
    "msd_analysis",
    "DiffusionAnalysisResult",
    "diffusion_analysis",
    # Free energy
    "wham_analysis",
    "wham_analysis_isolated",
    "mbar_analysis",
    "mbar_analysis_isolated",
    # Clustering
    "KmeansClusteringResult",
    "kmeans_clustering",
    # ATDYN engine
    "AtdynMDResult",
    "AtdynMinResult",
    "run_atdyn_md",
    "run_atdyn_min",
    "run_atdyn_md_isolated",
    "run_atdyn_min_isolated",
    # Utilities
    "extract_model_blocks",
]
