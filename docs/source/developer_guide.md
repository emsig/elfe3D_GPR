# Developer Guide

This page describes the repository structure and the development workflow for `elfe3D_GPR`.

## Repository structure

- `elfe3D_GPR/` — Fortran source tree and solver executable
- `io/` — Python model builders, input writers, execution helpers, and output readers
- `examples/` — example notebooks demonstrating typical use cases
- `tests/` — test harnesses and validation utilities
- `validation/` — reference datasets and validation scripts

## Building the project

The core procedure is:

1. compile the Fortran solver in `elfe3D_GPR/`
2. use the Python I/O layer to generate input files
3. run TetGen on the generated `.poly` file
4. execute the solver

## Documentation development

The docs live in `docs/source/` and are built with Sphinx.
The content includes user-facing guides, notebook workflows, and a scientific theory page.

## Future work

- package the Python I/O layer as an installable module under the namespace `elfe3d_gpr`
- integrate Doxygen output into RTD
- add automated Python API reference once the package layout is stable
- expand the example notebook documentation when the example workflows are finalized
