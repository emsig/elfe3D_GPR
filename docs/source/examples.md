# Example Notebooks

The example notebooks demonstrate the major modelling capabilities of `elfe3D_GPR`.
They provide hands-on workflows and show how to use the Python I/O wrapper with the core solver.

## Why these examples exist

Each notebook showcases a different modelling scenario and a different aspect of the workflow:

- `01_homogeneous_free-space.ipynb` — beginner tutorial for the homogeneous free-space case, verifying the basic source, mesh, solver, and postprocessing pipeline
- `02_homogeneous_earth.ipynb` — homogeneous earth layer, demonstrating a simple subsurface model
- `03_two_layered_earth.ipynb` — two-layer earth, showing layer interfaces and mesh adaptation
- `04_anomaly_sphere.ipynb` — single spherical anomaly, illustrating heterogeneity and target response
- `05_anomaly_multiple.ipynb` — multiple anomalies, demonstrating more complex subsurface structure

## What the notebooks show

The notebooks are designed to highlight:

- how to define surveys with the `elfe3d_gpr` Python I/O wrapper
- how to generate solver input files from Python
- how to run TetGen and the Fortran solver
- how to inspect and interpret the resulting field outputs

The notebooks now use the installed package namespace `elfe3d_gpr`, so they work cleanly after `pip install -e .`.

## Using the notebooks

Open the notebook in `examples/` with Jupyter and follow the guided workflow.
Use `01_homogeneous_free-space.ipynb` as the first notebook to verify the complete setup and solver execution.

These notebooks are intended to teach representative workflows rather than to serve as an exhaustive API reference.
