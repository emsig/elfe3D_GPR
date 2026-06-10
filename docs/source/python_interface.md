# Python I/O Wrapper {#python_interface}

The Python I/O wrapper is the primary user-facing entry point for notebooks and scripted workflows.

## Main user modules

- `elfe3d_gpr.runner`
  - `ProjectPaths` — path configuration for the solver executable and working directories
  - `run_tetgen()` — launch TetGen from Python
  - `run_solver()` — run the compiled `elfe3d_gpr` executable

- `elfe3d_gpr.inputs.survey`
  - `IOConfig` — manages input and output folder locations
  - `GPRSurvey` — overall survey definition and builder

## How to use the wrapper

Install the package in editable mode from the repository root:

```bash
pip install -e .
```

Then import from the supported package namespace:

```python
from elfe3d_gpr.runner import ProjectPaths, run_tetgen, run_solver
from elfe3d_gpr.inputs.survey import GPRSurvey
```

### Recommended usage example

```python
from pathlib import Path
from elfe3d_gpr.runner import ProjectPaths, run_tetgen, run_solver
from elfe3d_gpr.inputs.survey import GPRSurvey

master_dir = Path('..') / 'elfe3D_GPR'
master_dir = master_dir.resolve()

paths = ProjectPaths(
    master_dir=master_dir,
    exec_rel='elfe3d_gpr',
    use_wsl=True,
)

survey = GPRSurvey.build(
    experiment_name='air_only',
    base_dir=master_dir,
    air_eps_r=1.0,
    air_sigma=1e-16,
    layer_thicknesses=[0.3],
    layer_eps_r=[1.0],
    layer_sigma=[1e-16],
    layer_mu_r=[1.0],
    layer_sigma_m=[0.0],
    ricker_central_f=100e6,
    antenna_position=[0.0, 0.0, 0.025],
    num_receivers_inline=48,
)

survey.generate()
run_tetgen(paths, survey.io.poly_file)
run_solver(paths, survey)
```

The first example notebook, `examples/01_homogeneous_free-space.ipynb`, shows the recommended package-based workflow.
Legacy imports from `io` may still work when the repository root is on `PYTHONPATH`, but the supported package namespace is `elfe3d_gpr`.

## Notes on page structure

This page focuses on the Python I/O wrapper usage and import conventions.
The `workflow` page explains the repository architecture and the end-to-end process from model definition to solver execution.
Keeping them separate is helpful because readers who want code examples can go directly to `python_interface.md`, while readers who want the overall process can read `workflow.md`.

## Cross-references

Use `:ref:` labels when linking between pages in MyST/Sphinx.
For example:

- See :ref:`workflow` for the high-level repository workflow.
- See :doc:`quickstart` for a minimal step-by-step run.

If you prefer a simple link, use Markdown syntax:

```markdown
[Quick Start](quickstart.md)
```

But for RTD, `:ref:` and `:doc:` produce more robust links when the target page is in the same documentation tree.
