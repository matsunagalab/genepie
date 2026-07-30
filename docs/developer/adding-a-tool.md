# Adding a New Analysis Tool

There are two paths, depending on whether the CLI analysis already supports the
`trj_source_mod` abstraction.

## Option A — Unified architecture (recommended)

Use this when a CLI analysis already exists with `trj_source_mod` support (the
RMSD / RG / DRMS pattern). You reuse the CLI's routine directly, so CLI and Python
results are identical.

1. Export `analyze_<tool>_unified()` from the CLI's `*_analyze.fpp`, taking
   primitive arguments (arrays + scalars).
2. Create `src/analysis/interface/python_interface/<tool>_c_mod.fpp` that calls
   the unified function via `init_source_memory()` + `init_sink_array()`:

   ```
   subroutine <tool>_analysis_c(...) bind(C, name="<tool>_analysis_c")
     use <tool>_analyze_mod, only: analyze_<tool>_unified
     use trj_source_mod
     use result_sink_mod
     call init_source_memory(source, coords, boxes, natom, nframe, period)
     call init_sink_array(sink, result_ptr, result_size)
     call analyze_<tool>_unified(source, sink, ...)
     call finalize_source(source)
     call finalize_sink(sink)
   end subroutine
   ```

3. Add the C signature to `src/genepie/libgenesis.py` (`argtypes` / `restype`).
4. Add a Python wrapper in `src/genepie/analysis/<tool>.py` and re-export it from
   `genepie/analysis/__init__.py` (which flows through to `genesis_exe`).
5. Add a regression test `src/genepie/tests/test_<tool>.py` and wire it into
   `all_run.sh` / the CI workflows.

**Lazy support** comes almost for free: add a `<tool>_analysis_lazy_c` that uses
`init_source_lazy_dcd(...)` and loop with `get_next_frame`. Then let the Python
wrapper dispatch on `trajs.is_lazy` (see `rmsd.py` for the template).

## Option B — Separate implementation

Use this when the unified pattern does not fit (the HB / WHAM / MBAR / k-means
pattern).

1. Create `<tool>_c_mod.fpp` with a `bind(C)` interface.
2. Create `<tool>_impl.fpp` for the analysis implementation (or reuse existing
   GENESIS code).
3. Update `src/analysis/interface/python_interface/Makefile.am` to include the new
   `.fpp` files.
4. Add the C signature to `libgenesis.py` and a wrapper to
   `genepie/analysis/<tool>.py`.
5. Add `src/genepie/tests/test_<tool>.py` and wire it into CI.

## Python wrapper conventions

- Convert the molecule once: `mol_c = molecule.to_SMoleculeC()` and free it in a
  `finally:` with `deallocate_s_molecule_c`.
- Pre-allocate result arrays in NumPy and pass `.ctypes.data_as(c_void_p)` for
  zerocopy output.
- Wrap the call in `with fortran_status() as (status, msg, msglen):` so Fortran
  errors surface as the typed exceptions in [the error model](../reference/errors.md).
- Return a small `namedtuple` so results are self-describing.
- Resolve user-facing enums with the case-insensitive helpers in `_common.py`.

## Rebuild matrix

| Change | Command |
|--------|---------|
| Python only (`.py`) | none — instant |
| Fortran (`.fpp`) | `make` |
| New Fortran files | `make clean && make` |
| `configure.ac` | `autoreconf -fi && ./configure ... && make` |

See the [CLAUDE.md](https://github.com/matsunagalab/genepie/blob/main/CLAUDE.md)
developer guide for the fully worked signatures.
