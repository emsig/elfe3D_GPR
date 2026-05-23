# Installation

This page describes the core installation steps for `elfe3D_GPR` and the documentation build.

## System prerequisites

- Fortran compiler and environment for the core solver
- `tetgen` for mesh generation
- `MUMPS` for the direct linear algebra solver
- Python 3.11+ for notebook-driven I/O and documentation building

## Build the Fortran solver

The core solver lives in `elfe3D_GPR/elfe3d_gpr` and is built from the Fortran sources in `elfe3D_GPR/`.
Follow the instructions in the repository or the original `elfe3D` documentation for compiling with your compiler and linking against `tetgen` and `MUMPS`.

## Python I/O layer

The Python I/O layer is located in the `io/` folder. It is currently a curated set of scripts and modules used by notebooks for model setup, mesh generation, solver execution, and result handling.

The Python I/O layer is not yet published as a standalone installable package. For now, use the repository root as the working directory so the modules can be imported directly.

## Documentation build

The documentation site is built with Sphinx and MyST.
Install the documentation dependencies using:

```bash
python -m pip install sphinx myst-parser myst-nb sphinx-rtd-theme
```

Then build the site locally from the repository root:

```bash
cd docs
sphinx-build -b html source _build/html
```
