# Python Interface

The Python interface is the primary user-facing entry point for notebooks and scripted workflows.

## Main user modules

- `io.runner`
  - `ProjectPaths` — path configuration for the solver executable and working directories
  - `run_tetgen()` — launch TetGen from Python
  - `run_solver()` — run the compiled `elfe3d_gpr` executable

- `io.inputs.survey`
  - `IOConfig` — manages input and output folder locations
  - `GPRSurvey` — overall survey definition and builder

## How to use the interface

The current library is not yet a separate installable package, so import modules directly from the repository root.
For example:

```python
from io.runner import ProjectPaths, run_tetgen, run_solver
from io.inputs.survey import GPRSurvey
```

To document the Python I/O layer, we maintain a curated reference rather than a full auto-generated API. This allows the docs to stay aligned with the notebooks and the current development state.
