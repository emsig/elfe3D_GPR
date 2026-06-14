elfe3D_GPR Theory
=================

.. note::

    Much of the rigorous mathematical derivation of the:

    1. finite element formulation of the GPR wave equation, 
    2. resulting discrete system of equations, and 
    3. PML formulation that is implemented for the unstructured tetrahedral mesh of ``elfe3D_GPR`` 
   
    is discussed and developed in the thesis work [SIN2025]_. Please refer to it for complete details.

Following is a snippet of the theoretical background that ``elfe3D_GPR`` is based on. It will be completed soon.

Physics of Electromagnetic Waves for GPR
----------------------------------------

``elfe3D_GPR`` solves the full Maxwell's equations for wave-regime electromagnetism. This enables it to 
solve GPR problems in heterogeneous subsurface models.

First-Order Edge-Based Finite Elements
--------------------------------------

The software uses first-order edge-based FE to discretize the 3D model. Refer to [JIN2008]_ and [JIN2015]_ for 
complete analysis of the finite element method for electromagnetism, specifically for wave-regime problems. 
For ``elfe3D_GPR``, the key FE implementation details are as follows:

1. **Order** refers to the degree of interpolation. First-order means linear interpolation. This has been proven adequate for many physical problems using tetrahedral mesh with a suitable level of discretization.
2. **Edge-Based** refers to where the unknowns are formulated in the finite element mesh. A tetrahedron has 6 edges, which means per-element, there will be 6 unknowns in the system of linear equations. Multiply that by 6 for each of the electromagnetic field components :math:`E_x, E_y, E_z, H_x, H_y, H_z`, and you have 36 unknowns per element. Worth noting is that since a majority of elements are usually connected to each other, the total number of unknowns does not scale linearly to the number of elements in the model (especially including the boundary conditions).
3. The meshes ``elfe3D_GPR`` uses are **unstructured** and produced by ``tetgen``. Unstructured means they use irregularly connected elements to fit complex geometries. Unlike the regular grids of structured meshes, unstructured ones naturally favor modelling complex geometry with a high degree of accuracy.

Absorbing Boundary: Perfectly Matched Layers
--------------------------------------------

Our finite element meshes will have a finite length in all three axes. As such, **Perfectly Matched Layers (PML)** are applied 
around the domain boundaries to absorb outgoing waves. This ensures that they do not *reflect waves back from the truncated computational domain* [BER1994]_.
The specific type of PML is the **Uniaxial-PML** [BER2007]_, [PEK1995]_, [PLE2022]_. Unlike traditional U-PMLs, ``elfe3D_GPR`` implements an exact decay function that 
has been extensively studied over the years [BER2004]_, [OZG2023]_ and proven efficient for Finite Element meshes [SIN2025]_, [DIN2025]_, [FEN2019]_.