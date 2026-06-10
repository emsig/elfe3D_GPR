# Installation

This page describes the core installation and build steps for `elfe3D_GPR`.

## System prerequisites

- Fortran compiler for the core solver (`gfortran`, Intel Fortran, or compatible compiler)
- `tetgen` for mesh generation (tested with `tetgen` v1.5; v1.6 may not be compatible)
- `MUMPS` for the direct linear algebra solver
- Python 3.10+ for the notebook-driven Python I/O wrapper and docs build

## Verify Python setup

Check the active Python environment and package tools:

```bash
python --version
python -m pip --version
python -m pip install --upgrade pip setuptools wheel
```

## Build the Fortran solver

The core solver source is in `elfe3D_GPR/`.
Compile the solver using the provided build support and link against `tetgen` and `MUMPS`.
If the build does not succeed with your compiler, verify the `Makefile` settings in `elfe3D_GPR/` and ensure the external libraries are available.

## Install the Python I/O wrapper

The Python I/O layer is implemented in the `io/` folder and packaged under the namespace `elfe3d_gpr`.
Install it from the repository root:

```bash
pip install -e .
```

The supported import namespace is:

```python
from elfe3d_gpr.runner import ProjectPaths, run_tetgen, run_solver
from elfe3d_gpr.inputs.survey import GPRSurvey
```

Legacy imports from `io` may still work when the repository root is on `PYTHONPATH`, but the supported package namespace is `elfe3d_gpr`.

## Build the documentation

Install the Sphinx and MyST dependencies:

```bash
python -m pip install sphinx myst-parser myst-nb sphinx-rtd-theme
```

Build the documentation from the repository root:

```bash
cd docs
sphinx-build -b html source _build/html
```

If `doxygen` is installed and available on PATH, the Sphinx build will also generate the Fortran API documentation into `docs/_build/html/doxygen`.
