# Python API Reference

This page provides a curated reference for the Python modules that are currently used to build and run `elfe3D_GPR` simulations.

## Current strategy

The `io/` folder is not yet a complete installable Python package. For this reason, the documentation focuses on the main user-facing entry points rather than attempting a full auto-generated API.

## Key Python entry points

- `io.runner.ProjectPaths`
- `io.runner.run_tetgen`
- `io.runner.run_solver`
- `io.inputs.survey.IOConfig`
- `io.inputs.survey.GPRSurvey.build`

## Recommended reference approach

Once the Python I/O layer is packaged, it will be possible to add automatic `autodoc` pages for each module. Until then, keep this reference hand-curated and aligned with the example notebooks.

## Notes for future automation

The easiest automation path is to make `io/` importable from the repo root or to add package metadata such as `pyproject.toml`. That will allow Sphinx to import the modules and publish `autodoc`-generated documentation for the current API.
