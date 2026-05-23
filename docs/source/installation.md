# Installation

This page describes the core installation steps for `elfe3D_GPR` and the documentation build.

## System prerequisites

- Fortran compiler and environment for the core solver
- `tetgen` for mesh generation
- `MUMPS` for the direct linear algebra solver
- Python 3.10+ for notebook-driven I/O and documentation building

## Build the Fortran solver

The core solver lives in `elfe3D_GPR/elfe3d_gpr` and is built from the Fortran sources in `elfe3D_GPR/`.
Follow the instructions in the repository or the original `elfe3D` documentation for compiling with your compiler and linking against `tetgen` and `MUMPS`.

## Python I/O layer

The Python I/O layer is located in the `io/` folder, and it is packaged as `elfe3d-gpr` for installation.
The source files remain in `io/`, but the supported import namespace is `elfe3d_gpr`.

Install the Python package from the repository root:

```bash
pip install -e .
```

This installs the package and makes the I/O helpers available as:

```python
from elfe3d_gpr.runner import ProjectPaths, run_tetgen, run_solver
from elfe3d_gpr.inputs.survey import GPRSurvey
```

Legacy imports from `io` continue to work when the repository root is on `PYTHONPATH`, but the supported package namespace is `elfe3d_gpr`.

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
