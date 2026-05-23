# Quick Start

This section shows the minimal steps to run `elfe3D_GPR` using the Python notebook workflow.

## Quick start workflow

1. Build the Fortran executable for the solver.
2. Launch a Jupyter notebook or Python session from the repository root.
3. Use the Python I/O layer in `io/` to define the survey, build input files, run TetGen, and execute the solver.
4. Inspect the output files in `out_<experiment_name>`.

## Minimal example

```python
from io.runner import ProjectPaths, run_tetgen, run_solver
from io.inputs.survey import GPRSurvey

paths = ProjectPaths(master_dir=r'F:\Projects\\EMGeoInversion\\elfe3D_GPR')
survey = GPRSurvey.build(
    experiment_name='air_only',
    layer_thicknesses=[],
    layer_eps_r=[],
    layer_sigma=[],
    air_eps_r=1.0,
    antenna_position=[0.0, 0.0, 0.025],
    num_receivers_inline=5,
)
survey.generate()

run_tetgen(paths, survey.io.poly_file)
run_solver(paths, survey)
```

## Notes on the Python reference

The Python I/O layer is described in a curated reference page in `api_python.md`. It is designed for notebook-driven workflows rather than a fully installable API.
