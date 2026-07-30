# Design Details

Four ideas define how the interface layer connects Python to the GENESIS Fortran
engine.

## 1. Zerocopy memory management

Coordinate arrays and analysis results are **shared** between Python and Fortran
instead of copied. Python (NumPy) owns the memory; Fortran creates an alias to it
with `C_F_POINTER` and writes results directly into the NumPy array.

```
Copy-based (before):  Python -> copy -> Fortran -> copy -> Python   (2x memory)
Zerocopy   (after):   Python NumPy array  <->  Fortran C_F_POINTER  (1x, in place)
```

Array layout differs between the languages, so the wrappers transpose views
rather than copy data:

| Language | Layout | Element access |
|----------|--------|----------------|
| Python (NumPy) | `(nframe, natom, 3)`, C-order | `coords[frame, atom, 0]` |
| Fortran | `(3, natom, nframe)`, column-major | `coords_f(1, atom+1, frame+1)` |

Applied to: `crd_convert`, `rmsd_analysis`, `rg_analysis`, `drms_analysis`,
`trj_analysis`, `diffusion_analysis`.

**Safety rules:** the NumPy array must outlive the Fortran call, must not be
resized while Fortran holds the pointer, and must be contiguous
(`np.ascontiguousarray`).

```{admonition} SMolecule is copied, not zerocopy
:class: note
`SMolecule` is intentionally copy-based: it mixes string fields (need ASCII
encoding), 32-bit Fortran integers vs. NumPy int64, and row-/column-major layout.
`to_SMoleculeC()` builds the C structure and `deallocate_s_molecule_c()` frees it.
Only coordinate data and results are zerocopy.
```

## 2. Lazy DCD loading

Instead of loading a whole trajectory, lazy mode reads one frame at a time. DCD
frames are fixed-size, so any frame is reachable in O(1) via a byte offset:

```fortran
header_size = 92 + 12 + 80*ntitle + 12
frame_size  = 3 * (8 + 4*natom)            ! X, Y, Z blocks
if (has_box) frame_size = frame_size + 56  ! box block
byte_offset = header_size + (N-1)*frame_size + 1   ! frame N, 1-indexed
read(unit, pos=byte_offset) ...                    ! Fortran stream I/O
```

Peak coordinate memory is therefore ~one frame regardless of trajectory length,
and because the frames flow into the same unified core, results are identical to
the in-memory path. See the [Lazy loading tutorial](../tutorials/02b_lazy_loading.ipynb)
for the memory measurement and both APIs.

## 3. CLI / Python unified core

`trj_source_mod` and `result_sink_mod` abstract I/O so `analyze_*_unified()` is
agnostic to its caller:

```fortran
subroutine rmsd_analysis_c(...) bind(C)
  use ra_analyze_mod, only: analyze_rmsd_unified   ! the CLI's routine
  use trj_source_mod
  use result_sink_mod

  call init_source_memory(source, coords, boxes, natom, nframe, period)
  call init_sink_array(sink, result_ptr, result_size)
  call analyze_rmsd_unified(source, sink, ref_coord, mass, ...)   ! same as CLI
  call finalize_source(source)
  call finalize_sink(sink)
end subroutine
```

Source types: `TRJ_SOURCE_FILE`, `TRJ_SOURCE_MEMORY`, `TRJ_SOURCE_LAZY_DCD`.
Sink types: `SINK_FILE`, `SINK_ARRAY`.

## 4. Structured error handling

Fortran return codes map to a typed Python exception hierarchy, and Fortran's
stderr is captured onto the exception. The library-mode error guard turns what
used to be `exit(1)` on bad input into a catchable `GenesisFortranError`. Full
tables and examples are in the [error model](../reference/errors.md).
