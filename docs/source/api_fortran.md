# Fortran API Reference

The Fortran source tree in `elfe3D_GPR/` contains the core solver implementation, and Doxygen comments are already present for the API.

## What to document

Fortran documentation should include:

- solver entry points
- input file semantics
- module descriptions
- core numerical routines

## Current plan

At release, the Doxygen-generated API can be linked into the RTD site or integrated using a Sphinx `breathe` bridge.

## Core modules

Important modules include:

- `elfe3d_gpr` (solver executable driver)
- `mod_define_model.f90`
- `mod_solvers.f90`
- `mod_calculate_matrices.f90`
- `mod_read_mesh.f90`
- `mod_error_estimates.f90`

## Notes

This page is a placeholder for the Fortran reference. It is intentionally concise while the Doxygen HTML generation is established.
