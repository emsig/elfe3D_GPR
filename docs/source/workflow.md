# Workflow

This page explains how the repository pieces fit together and how data flows from model definition to solver results.

## Repository layout

- `elfe3D_GPR/` — the Fortran source tree and solver executable
- `io/` — Python helpers for model generation, mesh preparation, solver execution, and result handling
- `examples/` — notebook-based example workflows and demonstration cases
- `tests/` — automated tests and validation scripts
- `validation/` — reference results and validation comparisons

## Typical workflow

1. Define the GPR survey in Python using `io.inputs.survey.GPRSurvey.build()`.
2. Generate the mesh input and field input files with `survey.generate()`.
3. Run TetGen on the generated `.poly` file using `io.runner.run_tetgen()`.
4. Execute the Fortran solver with `io.runner.run_solver()`.
5. Read and visualize the solver output from `out_<experiment_name>`.

## Key components

- `io.runner` contains the execution layer and host path handling.
- `io.inputs.survey` is the central model builder.
- `io.inputs.writeinputfiles` writes the solver input files.
- `io.inputs.writetetgenpoly` writes the TetGen geometry input.
