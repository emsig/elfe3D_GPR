# Overview

`elfe3D_GPR` is a 3D finite-element forward modelling framework for ground penetrating radar (GPR) simulations.
It builds on the earlier `elfe3D` code base and adds a Python I/O wrapper for model definition, mesh generation, and solver execution.

The project is intended for Earth scientists and geophysical modellers who need flexibility in subsurface geometry, source and receiver placement, and heterogeneous material distributions.

## What this project provides

- A Fortran-based core solver for the 3D total electric field forward problem
- Input file generation for structured surveys, layers, and anomalies
- TetGen mesh generation integration
- A Python I/O wrapper for simplifying model setup and execution
- Example notebooks that guide users through representative GPR workflows
- RTD-ready documentation for installation, usage, and theory

## Key references

- The original `elfe3D` solver and manual provide the foundation for the numerical and software design.
- This documentation is inspired by the more modern RTD structure found in projects such as MFEM and empymod.
- For scientific background and theoretical support, see: *[Thesis reference placeholder]*.

## How to use this documentation

1. Read `installation` to set up the environment and build the solver.
2. Follow `quickstart` for either the native Fortran run or the notebook-based workflow.
3. Use `workflow` to understand how the repository pieces fit together.
4. Use `python_interface` to learn the supported Python I/O wrapper imports.
5. Refer to `inputs_and_models` and `outputs` for model definitions and file formats.
6. Run the example notebooks in `examples` to verify the end-to-end workflow.

## Contact

For questions or support, contact p.rulff@tudelft.nl.
