Inputs and Models
=================

This page describes the input files and model components that ``elfe3D_GPR`` works with.
It also explains how the Python I/O wrapper automates the creation of those files from
high-level model settings.

Input Files for ``elfe3D_GPR``
------------------------------

``elfe3D_GPR`` requires a few input files in the simulation input directory.
The ``examples/in`` example folder provides a current reference for these files.

``elfe3D_input.txt``
~~~~~~~~~~~~~~~~~~~~

This is the main solver configuration file written in plain text. It contains the model domain, frequency list,
receiver locations, output file names, source definition, solver parameters, and PML settings. Each keyword in the table 
below is followed by a number or a string, while for some list variables, multiple lines of numbers follow.

+----------------------------------------------------+----------------------------------------------------------------------+
| Keyword                                            | Description                                                          |
+====================================================+======================================================================+
| ``solver``                                         | Integer solver selection. Default is ``2`` for ``MUMPS``.            |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``model_size``                                     | Domain bounds. Followed by two lines:                                |
|                                                    | ``x_min y_min z_min`` and ``x_max y_max z_max``.                     |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``num_freq``                                       | Number of frequencies to simulate. Followed by lines of:             |
|                                                    | individual frequencies in Hertz.                                     |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``num_rec``                                        | Number of receivers. Followed by lines of                            |
|                                                    | Receiver coordinates in ``x y z`` format, one line per receiver.     |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``output_E_file``                                  | Path to electric field output file (no extension).                   |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``output_H_file``                                  | Path to magnetic field output file (no extension).                   |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``source_type``                                    | Integer describing source geometry and source file usage.            |
|                                                    | 6 is the most flexible option and hence the default.                 |
|                                                    | Followed by two lines of coordinates of source_start and source_end. |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``current_direction``                              | 1 - default (towards positive axis), 0 for reverse.                  |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``source_moment``                                  | Source moment. (1 is default).                                       |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``PEC_present``                                    | ``0`` or ``1`` (False or True) for perfect electric conductor (PEC)  |
|                                                    | presence in the model domain, not the boundaries of the PML.         |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``num_PEC``                                        | Number of PEC objects (0 for no PECs).                               |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``model_file_name``                                | Mesh input file stem, typically ``GPR_model_<name>.``.               |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``maxRefSteps``                                    | Maximum adaptive mesh refinement steps.                              |
|                                                    | ``0`` disables refinement (default). Needs testing in ``elfe3D_GPR``.|
+----------------------------------------------------+----------------------------------------------------------------------+
| ``maxUnknowns``                                    | Maximum number of unknowns allowed during a simulation.              |
|                                                    | Useful for refinement.                                               |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``betaRef``                                        | Error threshold for refining elements. Refer [RUL2021]_ for details. |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``accuracyTol``                                    | Accuracy tolerance. Refer [RUL2021]_ for details.                    |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``vtkRef``                                         | Set to ``1`` to write VTK files during refinement.                   |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``errorEst_method``                                | Error estimator type. Refer [RUL2021]_ for details.                  |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``refStrategy``                                    | Refinement strategy. Refer [RUL2021]_ for details.                   |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``output_fields_vtk``                              | View electric and magnetic field components in a ``.vtk`` file.      |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``PML_present``                                    | ``1`` to enable the PML. Should not be changed for GPR simulations.  |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``PML_thickness``                                  | Thickness of the PML layer added around the model domain.            |
+----------------------------------------------------+----------------------------------------------------------------------+
| ``PML_decay_type``                                 | PML decay type, e.g. ``1`` for exact reciprocal decay.               |
|                                                    | Other decays given in [SIN2025]_.                                    |
+----------------------------------------------------+----------------------------------------------------------------------+

The provided ``in/elfe3D_input.txt`` example includes a single frequency, 48 receivers,
source segment coordinates, and PML settings for a small air-only domain.

``source.txt``
~~~~~~~~~~~~~~

The source file contains the endpoints of the segmented line source.
For a straight dipole source, it contains two points:

- first line: number of source points (``2`` for a single straight segment).
- next two lines: coordinates in the format :math:`x,\ y,\ z` for the source start and end.

In the example ``source.txt``, the source is defined from :math:`(-3\times10^{-4},\ 0.0,\ 0.025)` to
:math:`(3\times10^{-4},\ 0.0,\ 0.025)`.

``regionparameters.txt``
~~~~~~~~~~~~~~~~~~~~~~~~~

This file specifies the electromagnetic properties assigned to each tetrahedral region
attribute in the mesh. The format is:

- line 1: ``# eleattr`` (number of mesh element attributes, where attributes equals regions).
- line 2: number of region entries.
- line 3: column labels ``# eleattr rho mu_r epsilon_r``
- following lines: one entry per region attribute with: attribute_num, rho, mu_r, epsilon_r.

`rho` is the electrical resistivity in Siemens/m, `mu_r` is the relative magnetic permeability, and
`epsilon_r` is the relative electric permittivity (dielectric constant) of a medium.

In the simple ``in`` example, all region attributes are assigned as
free-space parameters with very large resistivity and unit relative permittivity and permeability.

.. note::
    You might notice that there are quite a few region parameter entries in this file. It is due to the fact that each 
    PML block needs its own definitions as per ``tetgen`` meshing strategy. For a model with a single air-earth interface, 
    there are 34 PML blocks wrapping around the model domain. The air block and the earth blocks are the other two regions, 
    for a total of 36 `regions`.

``GPR_model_<name>.poly``
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``.poly`` file contains the geometry definition for ``tetgen``. It describes the domain
boundary, regions, and point markers required to build the unstructured tetrahedral mesh.

You can find the documentation on ``.poly`` file structure on `TetGen's Manual <https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual006.html>`_.

.. warning::
    If you are using ``elfe3D_GPR`` without the Python I/O module, you will have to generate 
    the input ``.poly`` file yourself, including the definition for PML blocks around the model domain. 
    This could get quite complex, especially with multi-layered media that extend beyond the model domain. 

    As such, we recommend using the Python I/O module that generates this ``.poly`` file automatically from the model domain, 
    layers, anomalies, source geometry, receivers, and PML settings. It also encapsulates some of the experimentation and 
    design choices that were made during the work on the thesis [SIN2025]_ regarding ideal mesh element sizes for source dipole, 
    different dielectric media, receiver points and surrounding tetrahedra, PML block and the layer interface modelling.

TetGen reads the ``.poly`` file and produces the mesh files required by ``elfe3D_GPR``:
``.node``, ``.ele``, ``.face``, ``.neigh``, ``.edge``, and ``.vtk``.

How the Python I/O module automates this
----------------------------------------

Many of the Python I/O module's ``GPRSurvey.build()`` arguments or less 1-on-1 match to the inputs defined in ``elfe3D_input.txt``. 
They also include details for the geophysical geometry and material model.

Key mappings from Python to Fortran input files:

- ``experiment_name``: naming of ``in_<experiment_name>`` and ``out_<experiment_name>``.
- ``base_dir``: root directory for input and output folders.
- ``x_e``, ``y_e``, ``z_e``: ``model_size`` domain extents.
- ``f_list``: ``num_freq`` and the frequency list in ``elfe3D_input.txt``.
- ``num_receivers_inline``, ``num_receivers_endfire``, ``num_receivers_oblique``: receiver coordinates in ``num_rec``.
- ``output_fields_vtk``: VTK field output control.
- ``source_type``, ``antenna_position``, ``current_direction``, ``num_segments``, ``s_f``: source geometry and ``source.txt`` endpoints.
- ``solver_type``, ``max_ref_steps``, ``max_unknowns``, ``beta_ref``, ``accuracy_tol``, ``vtk``, ``error_est_method``, ``ref_strategy``
  : solver fields in ``elfe3D_input.txt``.
- ``num_pml_layers``, ``pml_layer_thickness``, ``pml_decay_type``: PML settings in ``elfe3D_input.txt``.

The notebooks in ``examples/`` directory demonstrate this mapping well:
``GPRSurvey.build(...)`` gathers the same inputs used in the text file format,
and ``survey.generate()`` writes ``elfe3D_input.txt``, ``source.txt``, ``regionparameters.txt``,
and the TetGen ``.poly`` file.

Layers, anomalies, and PML in the Python Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Layers**:

- The air layer is always present and uses ``air_eps_r = 1.0`` and ``air_sigma = 1e-16`` variables (does not need to be explicitly given).
- Earth layers are defined by these list variables:
    - ``layer_thicknesses = []``, 
    - ``layer_eps_r = []``, 
    - ``layer_sigma = []``, 
    - ``layer_mu_r = []``, 
    - ``layer_sigma_m = []``.

- Each layer becomes a ``GeoLayer`` Python class object and is included in the ``.poly`` mesh regions. All these variables should have the same number of entries for consistency.

**Anomalies**:

Two types of anomalies are available in ``elfe3D_GPR`` (see ``io/inputs/anomalies.py``):

- ``BoxAnomaly(x: tuple[float, float], y: tuple[float, float], z: tuple[float, float], properties: tuple[float, float, float, float])`` —
  rectangular prism defined by axis extents and material properties ``(eps_r, sigma, mu_r, sigma_m)``.
- ``SphereAnomaly(center: tuple[float, float, float], radius: float, properties: tuple[float, float, float, float], subdivision_levels: int = 4)`` —
  spherical anomaly tessellated as an icosphere; ``subdivision_levels`` controls mesh detail.

You can supply multiple anomalies to ``GPRSurvey.build()`` via ``anomalies=[...]``. Anomalies are inserted into the generated
``.poly`` geometry and assigned region properties that appear in ``regionparameters.txt``.

**PML**:

- The wrapper uses ``PMLConfig`` class and writes ``PML_present``, ``PML_thickness``, and ``PML_decay_type`` variables
  into ``elfe3D_input.txt``. These inputs also help define the PML blocks for mesh generation.
- ``num_pml_layers`` and ``pml_layer_thickness`` control the mesh padding around the domain. Due to the exact-decay and unstructured meshing, it is not necessary to increase ``num_pml_layers`` to more than 1.
- The current default decay type is ``1``, corresponding to the exact reciprocal decay implementation [SIN2025]_, [DIN2025]_.

Next, we will discuss the outputs that ``elfe3D_GPR`` produces.

