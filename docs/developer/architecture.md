# Architecture

genepie is a **thin interface layer** over the GENESIS Fortran engine. It adds no
new science: every number comes from the same Fortran routines the GENESIS CLI
uses. Three layers are connected by standard interop mechanisms.

```{mermaid}
flowchart TD
  subgraph User["User layer"]
    U["Python / Jupyter / NumPy / PyTorch"]
  end
  subgraph Interface["Interface layer — genepie"]
    A["genepie.analysis.* (Python wrappers)"]
    L["libgenesis.py (ctypes signatures)"]
    LL["libloader.py (finds libpython_interface)"]
  end
  subgraph Engine["Engine layer — GENESIS (Fortran)"]
    C["*_c_mod.fpp (bind(C) wrappers)"]
    K["analyze_*_unified() cores"]
  end
  U --> A --> L --> C --> K
  LL -.loads.-> C
```

## The Python package layout

The interface used to be one large `genesis_exe.py`; it is now split into the
`genepie.analysis` package, one module per tool family. `genesis_exe` re-exports
everything, so `genesis_exe.rmsd_analysis` and
`from genepie.analysis import rmsd_analysis` are the same function.

| Module | Contents |
|--------|----------|
| `genepie/analysis/converter.py` | `crd_convert`, `crd_convert_info`, `selection` (+ `lazy=True`) |
| `genepie/analysis/rmsd.py` | `rmsd_analysis`, `rmsd_analysis_lazy` |
| `genepie/analysis/rg.py`, `drms.py`, `trj.py`, `avecrd.py`, `hbond.py` | structural analyses |
| `genepie/analysis/dynamics.py` | `msd_analysis`, `diffusion_analysis` |
| `genepie/analysis/free_energy.py` | `wham_analysis`, `mbar_analysis` |
| `genepie/analysis/clustering.py` | `kmeans_clustering` |
| `genepie/analysis/atdyn.py` | `run_atdyn_md/min` (+ `_isolated`) |
| `genepie/analysis/_common.py` | enum resolution, filename packing shared by wrappers |

Supporting modules (unchanged by the split): `s_molecule.py`,
`s_trajectories.py`, `libgenesis.py`, `libloader.py`, `exceptions.py`, the
`*_validators.py` files, and `output_capture.py`.

## The unified core (CLI = Python)

The key architectural idea is that the CLI tool and the Python interface call the
**same** analysis routine, `analyze_*_unified()`. Only the input and output are
abstracted, so results are identical by construction and bugs are fixed in one
place.

```{mermaid}
flowchart TD
  subgraph Sources["trj_source_mod"]
    SF["TRJ_SOURCE_FILE (CLI: read file)"]
    SM["TRJ_SOURCE_MEMORY (Python: NumPy array)"]
    SL["TRJ_SOURCE_LAZY_DCD (stream from disk)"]
  end
  SF --> CORE["analyze_*_unified()"]
  SM --> CORE
  SL --> CORE
  CORE --> Sinks
  subgraph Sinks["result_sink_mod"]
    KF["SINK_FILE (CLI: write file)"]
    KA["SINK_ARRAY (Python: zerocopy NumPy)"]
  end
```

- **`trj_source`** abstracts where frames come from: a file (CLI), a NumPy array
  in memory (Python), or on-demand reads from a DCD (lazy).
- **`result_sink`** abstracts where results go: an output file (CLI) or a
  zerocopy NumPy array (Python).

Tools built on this pattern (RMSD, RG, DRMS) share one implementation across CLI
and Python. See [Design details](design.md) for the zerocopy and lazy mechanics,
and [Adding a tool](adding-a-tool.md) for how to extend it.
