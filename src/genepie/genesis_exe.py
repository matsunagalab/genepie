"""Backwards-compatible flat view of the analysis wrappers.

The implementations live in :mod:`genepie.analysis`, one module per tool family.
This module re-exports them so that existing ``genesis_exe.<function>`` calls
keep working unchanged.
"""

from .analysis import (  # noqa: F401
    AtdynMDResult,
    AtdynMinResult,
    AvecrdAnalysisResult,
    CrdConvertInfo,
    DiffusionAnalysisResult,
    DrmsAnalysisResult,
    KmeansClusteringResult,
    MsdAnalysisResult,
    RgAnalysisResult,
    RmsdAnalysisResult,
    RmsdLazyAnalysisResult,
    TrjAnalysisResult,
    avecrd_analysis,
    crd_convert,
    crd_convert_info,
    diffusion_analysis,
    drms_analysis,
    extract_model_blocks,
    hb_analysis,
    kmeans_clustering,
    mbar_analysis,
    msd_analysis,
    rg_analysis,
    rmsd_analysis,
    rmsd_analysis_lazy,
    run_atdyn_md,
    run_atdyn_md_isolated,
    run_atdyn_min,
    run_atdyn_min_isolated,
    selection,
    selection_func,
    trj_analysis,
    wham_analysis,
)
from .analysis import __all__ as _analysis_all

# Names that predate the split and are still imported from here by name.
from .constants import (  # noqa: F401
    FittingMethod,
    PBCCMode,
    TrjFormat,
    TrjType,
)
from ._fortran import (  # noqa: F401
    DEFAULT_MSG_LEN as _DEFAULT_MSG_LEN,
    get_msg_len,
    make_msgbuf,
)
from .libgenesis import LibGenesis  # noqa: F401

__all__ = list(_analysis_all) + [
    "LibGenesis",
    "FittingMethod",
    "PBCCMode",
    "TrjFormat",
    "TrjType",
    "get_msg_len",
    "make_msgbuf",
]
