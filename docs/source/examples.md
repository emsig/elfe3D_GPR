# Example Notebooks

The example notebooks demonstrate the major modelling capabilities of `elfe3D_GPR`.
They provide hands-on workflows rather than full reference documentation.

## Why these examples exist

Each notebook showcases a different modelling scenario and a different aspect of the workflow:

- `air.ipynb` — air-only domain, verifying the basic source and receiver pipeline
- `02_homogeneous_earth` — homogeneous earth layer, demonstrating a simple subsurface model
- `03_two_layered_earth` — two-layer earth, showing layer interfaces and mesh adaptation
- `04_anomalous_earth` — anomaly cases, illustrating heterogeneity and target response

## What the notebooks show

The notebooks are designed to highlight:

- how to define surveys with `GPRSurvey`
- how to generate solver input files from Python
- how to run TetGen and the Fortran solver
- how to inspect and interpret the resulting field outputs

## Using the notebooks

The notebooks are stored in the `examples/` folder and may be copied or extracted into other study folders. Their purpose is to demonstrate representative workflows rather than to act as standalone reference manuals.
