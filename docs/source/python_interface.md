# Python Interface

The Python interface is the primary user-facing entry point for notebooks and scripted workflows.

## Main user modules

- `elfe3d_gpr.runner`
  - `ProjectPaths` — path configuration for the solver executable and working directories
  - `run_tetgen()` — launch TetGen from Python
  - `run_solver()` — run the compiled `elfe3d_gpr` executable

- `elfe3d_gpr.inputs.survey`
  - `IOConfig` — manages input and output folder locations
  - `GPRSurvey` — overall survey definition and builder

## How to use the interface

Install the package in editable mode from the repository root:

```bash
pip install -e .
```

Then import from the package namespace:

```python
from elfe3d_gpr.runner import ProjectPaths, run_tetgen, run_solver
from elfe3d_gpr.inputs.survey import GPRSurvey
```

Legacy imports from `io` continue to work when the repository root is on `PYTHONPATH`, but the supported package namespace is `elfe3d_gpr`.

To document the Python I/O layer, we maintain a curated reference rather than a full auto-generated API. This allows the docs to stay aligned with the notebooks and the current development state.
