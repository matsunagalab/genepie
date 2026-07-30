# genepie

**genepie** is the Python interface to [GENESIS](https://www.r-ccs.riken.jp/labs/cbrt/)
(GENeralized-Ensemble SImulation System), a molecular dynamics package for
biomolecular systems. It exposes the GENESIS analysis tools and the ATDYN MD
engine as native Python functions that return NumPy arrays.

## Why genepie?

An MD study is not a single command — it is a workflow of *setup → simulation →
analysis*, with many analysis methods to choose from, and increasingly a final
step that feeds the results into AI/ML tooling. Doing this by stitching together
CLI tools and intermediate files is tedious and hard to reproduce.

genepie turns GENESIS into a **single, programmable Python workflow**: run
simulations, run analyses, and hand the numerical results straight to
NumPy / PyTorch / scikit-learn without leaving Python. It bridges the GENESIS
Fortran engine with researchers and the AI/ML ecosystem.

## Architecture

genepie is a thin interface layer over the existing GENESIS Fortran engine — it
does not reimplement any science. Three layers are connected by standard interop
mechanisms:

```{mermaid}
flowchart TD
  subgraph User["User layer"]
    U["Python script / Jupyter / NumPy / PyTorch"]
  end
  subgraph Interface["Interface layer"]
    G["genepie (this package)"]
  end
  subgraph Engine["Engine layer"]
    F["GENESIS (Fortran: MD &amp; analysis)"]
  end
  U -->|"ctypes (C &harr; Python)"| G
  G -->|"ISO_C_BINDING (Fortran &harr; C)"| F
```

The interface layer is built around four design ideas (see the
[Developer Guide](developer/design.md) for details):

- **Zerocopy memory** — NumPy arrays and Fortran share the same memory (no copy).
- **Lazy DCD loading** — read one frame at a time for large trajectories.
- **CLI/Python unified core** — the same Fortran routine serves both the CLI and
  Python, so results are identical.
- **Structured error handling** — a typed exception hierarchy with numeric error
  codes.

## How this site is organized

```{mermaid}
flowchart TD
  intro["Introduction (this page)"] --> install["Installation"]
  install --> qs["Getting Started (BPTI, ~10 min)"]
  qs --> tut["Tutorials"]
  qs --> ref["Reference"]
  qs --> dev["Developer Guide"]
```

- **Getting Started** — install genepie and run a complete BPTI example end to end.
- **Tutorials** — focused, runnable chapters for each capability: loading
  molecules, trajectories, lazy loading, structural analysis, dynamics, free
  energy, clustering, the MD engine, and ML integration.
- **Reference** — the analysis API tables and the error/exception model.
- **Developer Guide** — how the interface layer is built and how to add a new tool.

```{admonition} Runnable examples
:class: tip
Every notebook in **Getting Started** and **Tutorials** is executed when this site
is built, so the plots, arrays, and structures you see are real outputs produced
on the trajectories bundled with genepie's test data — no large downloads
required. The only exception is the ML integration chapter, which ships with
pre-saved outputs because it needs PyTorch and a larger dataset.
```
