# elfe3D GPR

`elfe3D_GPR` is a 3D edge-based **Finite Element** (FE) software for **Ground Penetrating Radar** (GPR) geophysical heterogeneous models in the frequency domain.

## Getting started:

You find the `elfe3D_GPR` manual including installation instructions at this [link](). 

## Contributions:

`elfe3D_GPR` builds on the `elfe3D` software developed by Paula Rulff ([GitHub](https://github.com/emsig/elfe3D)) that simulates diffusive-field problems in Controlled-Source Electromagnetism. 

Both `elfe3D` and `elfe3D_GPR` are programmed in Fortran 90, that use:
- [tetgen](https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) (the Fortran version 1.5 currently), for producing unstructured tetrahedral mesh, and
- [MUMPS](https://mumps-solver.org/index.php) for solving the FE system of linear equations.

The major changes between `elfe3D` and `elfe3D_GPR` are as follows:

- `elfe3D_GPR` solves the full wave equation of electromagnetism based on Maxwell's equations (whereas `elfe3D` solved the approximated equation for diffusive-field regime of very low frequencies.)
- A **Perfectly Matched Layer (PML)** is added to absorb outgoing waves from the truncated computational model domain.

These two major changes have been rigorously developed and described in the master's thesis written by Chaitanya Singh, and supervised by Paula Rulff and Evert Slob ([thesis_doc](https://resolver.tudelft.nl/uuid:b883c3d6-beb2-4842-b867-21d0c777aff7)).

Since the thesis, `elfe3D_GPR` now also includes a Python I/O module for rapid and intuitive input model generation and output processing.

## Contact: 

If you would like to discuss the theoretical background of ``elfe3D_GPR``, suggest changes to the software, or would like to contribute, 
you are welcome to contact us at: [chaitanya.singh@northumbria.ac.uk](chaitanya.singh@northumbria.ac.uk).

## Credits:

If you publish results generated with `elfe3D_GPR`, please give credit to the `elfe3D_GPR` developers by citing:

> Singh, C.D. (2025). Frequency-Domain Wideband Ground Penetrating Radar Modelling: Using Finite Elements and Perfectly Matched Layers. Master's Thesis, Delft University of Technology. [https://resolver.tudelft.nl/uuid:b883c3d6-beb2-4842-b867-21d0c777aff7](https://resolver.tudelft.nl/uuid:b883c3d6-beb2-4842-b867-21d0c777aff7).

and `elfe3D` developers by citing:

> Paula Rulff, Laura M Buntin, Thomas Kalscheuer, Efficient goal-oriented mesh refinement in 3-D finite-element modelling adapted for controlled source electromagnetic surveys, Geophysical Journal International, Volume 227, Issue 3, December 2021, Pages 1624–1645, [https://doi.org/10.1093/gji/ggab264](https://doi.org/10.1093/gji/ggab264).

Do not forget to acknowledge `MUMPS` and `TetGen` developers!