# Example Notebooks

The example notebooks demonstrate the major modelling capabilities of `elfe3D_GPR`.
They provide hands-on workflows rather than full reference documentation.

## Why these examples exist

Each notebook showcases a different modelling scenario and a different aspect of the workflow:

- `01_wholespace_air.ipynb` — a guided beginner tutorial for the air-only case, verifying the basic source, mesh, solver, and postprocessing pipeline
- `02_homogeneous_earth.ipynb` — homogeneous earth layer, demonstrating a simple subsurface model
- `03_two_layered_earth.ipynb` — two-layer earth, showing layer interfaces and mesh adaptation
- `04_anomalous_earth.ipynb` — anomaly cases, illustrating heterogeneity and target response

## What the notebooks show

The notebooks are designed to highlight:

- how to define surveys with `GPRSurvey`
- how to generate solver input files from Python
- how to run TetGen and the Fortran solver
- how to inspect and interpret the resulting field outputs

These example notebooks now use the installed package namespace `elfe3d_gpr`, so they work cleanly after `pip install -e .`.

## Using the notebooks

The notebooks are stored in the `examples/` folder and may be copied or extracted into other study folders. Their purpose is to demonstrate representative workflows rather than to act as standalone reference manuals.
