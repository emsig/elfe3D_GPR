# elfe3D_GPR Theory

This page combines the theory and numerical methods behind the solver.

## Physical formulation

`elfe3D_GPR` solves the frequency-domain controlled-source electromagnetic problem in terms of the total electric field.
The physical model includes:

- dielectric permittivity
- electrical conductivity
- magnetic permeability
- conductive magnetic loss

The solver handles the coupled behaviour of the air and earth subsurface, including material contrasts and anomaly bodies.

## Numerical approach

The code uses first-order finite-element discretization on tetrahedral meshes.
Mesh generation is performed by TetGen from the geometry defined in the Python I/O layer.

## Mesh refinement and solver strategy

Adaptive mesh refinement is used to concentrate resolution where the fields vary most strongly.
The solver can refine the model iteratively based on error estimation.

## Boundary treatment

Perfectly matched layers (PMLs) are applied around the domain boundaries to absorb outgoing waves and mimic an open-region environment.

## Code architecture

The Fortran source tree implements the solver, mesh handling, and solver routines. The Python I/O layer prepares the problem setup and converts it into the solver's input format.
