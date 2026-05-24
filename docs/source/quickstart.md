# Quick Start

This section shows the minimal steps to run `elfe3D_GPR` using the Python notebook workflow.

## Quick start workflow

1. Build the Fortran executable for the solver.
2. Launch a Jupyter notebook or Python session from the repository root.
3. Install the Python I/O package and use `elfe3d_gpr` to define the survey, build input files, run TetGen, and execute the solver.
4. Inspect the output files in `out_<experiment_name>`.
5. Start with `examples/01_wholespace_air.ipynb` as the guided beginner notebook.

## Minimal example

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
    num_receivers_inline=5,
)

survey.generate()
run_tetgen(paths, survey.io.poly_file)
run_solver(paths, survey)
```

## Linking to other docs

- See :ref:`workflow` for the high-level repository and data flow.
- See :ref:`python_interface` for the Python package API and import conventions.

## Notes on the Python reference

The Python I/O layer is described in a curated reference page in `python_interface.md`.
This page is focused on the minimal user workflow, while `python_interface.md` is focused on the package API and supported import names.
