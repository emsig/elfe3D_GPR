# Python API Reference

This page provides a curated reference for the Python modules that are currently used to build and run `elfe3D_GPR` simulations.

## Current strategy

The `io/` folder is now packaged for installation, and the supported import namespace is `elfe3d_gpr`.
This page continues to provide a curated reference for the core user-facing entry points.

## Key Python entry points

- `elfe3d_gpr.runner.ProjectPaths`
- `elfe3d_gpr.runner.run_tetgen`
- `elfe3d_gpr.runner.run_solver`
- `elfe3d_gpr.inputs.survey.IOConfig`
- `elfe3d_gpr.inputs.survey.GPRSurvey.build`

## Recommended reference approach

Once the Python I/O layer is packaged, it will be possible to add automatic `autodoc` pages for each module. Until then, keep this reference hand-curated and aligned with the example notebooks.

## Notes for future automation

The easiest automation path is to make `io/` importable from the repo root or to add package metadata such as `pyproject.toml`. That will allow Sphinx to import the modules and publish `autodoc`-generated documentation for the current API.
