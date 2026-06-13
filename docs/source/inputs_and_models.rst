Inputs and Models
=================

This page describes the input files and model components needed by `elfe3D_GPR`.
It also explains how the Python I/O wrapper automates the creation of those files from
high-level model settings.

Input Files for ``elfe3D_GPR``
------------------------------

`elfe3D_GPR` requires several input files in the simulation input directory.
The `in_air_only` example folder provides a current reference for these files.

``elfe3D_input.txt``
~~~~~~~~~~~~~~~~~~~~

This is the main solver configuration file. It contains the model domain, frequency list,
receiver locations, output file names, source definition, solver parameters, and PML settings.
The file is plain text and uses keyword-value lines followed by optional blocks.

Common keywords found in ``elfe3D_input.txt``:

+----------------------+--------------------------------------------------------------------------+
| Keyword              | Description                                                              |
+======================+==========================================================================+
| ``solver``           | Integer solver selection. Current supported option is ``2`` for MUMPS.   |
+----------------------+--------------------------------------------------------------------------+
| ``model_size``       | Domain bounds. Followed by two lines:                                    |
|                      | ``x_min y_min z_min`` and ``x_max y_max z_max``.                         |
+----------------------+--------------------------------------------------------------------------+
| ``num_freq``         | Number of frequencies to simulate.                                       |
+----------------------+--------------------------------------------------------------------------+
| frequency list       | One frequency per line after ``num_freq``.                               |
+----------------------+--------------------------------------------------------------------------+
| ``num_rec``          | Number of receivers.                                                     |
+----------------------+--------------------------------------------------------------------------+
| receiver list        | Receiver coordinates in ``x y z`` format, one line per receiver.         |
+----------------------+--------------------------------------------------------------------------+
| ``output_E_file``    | Path to electric field output file (no extension).                       |
+----------------------+--------------------------------------------------------------------------+
| ``output_H_file``    | Path to magnetic field output file (no extension).                       |
+----------------------+--------------------------------------------------------------------------+
| ``source_type``      | Integer describing source geometry and source file usage.                |
+----------------------+--------------------------------------------------------------------------+
| source endpoints     | Two lines with source start and end coordinates.                         |
+----------------------+--------------------------------------------------------------------------+
| ``current_direction``| Integer indicating the source current direction.                         |
+----------------------+--------------------------------------------------------------------------+
| ``source_moment``    | Source moment ``m = I dl``.                                              |
+----------------------+--------------------------------------------------------------------------+
| ``PEC_present``      | ``0`` or ``1`` for perfect electric conductor presence in the model.     |
+----------------------+--------------------------------------------------------------------------+
| ``num_PEC``          | Number of PEC objects (usually ``0`` for no PECs).                       |
+----------------------+--------------------------------------------------------------------------+
| ``model_file_name``  | Mesh input file stem, typically ``GPR_model_<name>.``.                   |
+----------------------+--------------------------------------------------------------------------+
| ``maxRefSteps``      | Maximum mesh refinement steps. ``0`` disables refinement.                |
+----------------------+--------------------------------------------------------------------------+
| ``maxUnknowns``      | Maximum number of unknowns allowed during refinement.                    |
+----------------------+--------------------------------------------------------------------------+
| ``betaRef``          | Error threshold for refining elements.                                   |
+----------------------+--------------------------------------------------------------------------+
| ``accuracyTol``      | Accuracy tolerance (must be less than ``1``).                            |
+----------------------+--------------------------------------------------------------------------+
| ``vtkRef``           | Set to ``1`` to write VTK files during.                                  |
+----------------------+--------------------------------------------------------------------------+
| ``errorEst_method``  | Error estimator type.                                                    |
+----------------------+--------------------------------------------------------------------------+
| ``refStrategy``      | Refinement strategy.                                                     |
+----------------------+--------------------------------------------------------------------------+
| ``output_fields_vtk``| Set to ``1`` to write field outputs in VTK format.                       |
+----------------------+--------------------------------------------------------------------------+
| ``PML_present``      | ``0`` or ``1`` to enable or disable the PML.                             |
+----------------------+--------------------------------------------------------------------------+
| ``PML_thickness``    | Thickness of the PML layer added around the domain.                      |
+----------------------+--------------------------------------------------------------------------+
| ``PML_decay_type``   | PML decay type, e.g. ``1`` for exact reciprocal decay.                   |
+----------------------+--------------------------------------------------------------------------+

The current `in_air_only/elfe3D_input.txt` example includes a single frequency, 48 receivers,
source segment coordinates, and PML settings for a small air-only domain.

``source.txt``
~~~~~~~~~~~~~~

The source file contains the endpoints of the segmented line source.
For a straight dipole source, it contains two points:

- first line: number of source points (`2` for a single straight segment).
- next two lines: `x y z` coordinates for the source start and end.

In the example `source.txt`, the source is defined from `(-0.0003, 0.0, 0.025)` to
`(0.0003, 0.0, 0.025)`.

``regionparameters.txt``
~~~~~~~~~~~~~~~~~~~~~~~~~

This file specifies the electromagnetic properties assigned to each tetrahedral region
attribute in the mesh. The format is:

- line 1: `# eleattr`
- line 2: number of region entries
- line 3: column labels `# eleattr rho mu_r epsilon_r`
- following lines: one entry per region attribute with `attribute rho mu_r epsilon_r`

`rho` is the resistivity, `mu_r` is the relative magnetic permeability, and
`epsilon_r` is the relative permittivity.

In the current `in_air_only` example, all region attributes are assigned nearly
free-space parameters with very large resistivity and unit relative permittivity and permeability.

``GPR_model_<name>.poly``
~~~~~~~~~~~~~~~~~~~~~~~~~

The `.poly` file contains the geometry definition for TetGen. It describes the domain
boundary, regions, and point markers required to build the unstructured tetrahedral mesh.

The Python I/O wrapper generates this `.poly` file automatically from the model domain,
layers, anomalies, source geometry, receivers, and PML settings.
TetGen reads the `.poly` file and produces the mesh files required by `elfe3D_GPR`:
`.node`, `.ele`, `.face`, `.neigh`, `.edge`, and `.vtk`.

How the Python I/O module automates this
----------------------------------------

The Python wrapper is implemented in `io/inputs/survey.py` and the notebook
`examples/01_homogeneous_free-space.ipynb` shows the practical workflow.
The wrapper uses `GPRSurvey.build()` to create a complete model description,
then `survey.generate()` to write all required input files.

Key mappings from Python to Fortran input files:

- `experiment_name` → naming of `in_<experiment_name>` and `out_<experiment_name>`.
- `base_dir` → root directory for input and output folders.
- `x_e`, `y_e`, `z_e` → `model_size` domain extents.
- `f_list` → `num_freq` and the frequency list in `elfe3D_input.txt`.
- `num_receivers_inline`, `num_receivers_endfire`, `num_receivers_oblique` → receiver coordinates in `num_rec`.
- `output_fields_vtk` → VTK field output control.
- `source_type`, `antenna_position`, `current_direction`, `num_segments`, `s_f` → source geometry and `source.txt` endpoints.
- `solver_type`, `max_ref_steps`, `max_unknowns`, `beta_ref`, `accuracy_tol`, `vtk`, `error_est_method`, `ref_strategy`
  → solver fields in `elfe3D_input.txt`.
- `num_pml_layers`, `pml_layer_thickness`, `pml_decay_type` → PML settings in `elfe3D_input.txt`.

The notebook `01_homogeneous_free-space.ipynb` demonstrates this mapping directly:
`GPRSurvey.build(...)` gathers the same inputs used in the legacy text file format,
and `survey.generate()` writes `elfe3D_input.txt`, `source.txt`, `regionparameters.txt`,
and the TetGen `.poly` file.

Layers, anomalies, and PML in the Python wrapper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python module keeps these concepts explicit while automating the low-level file format.

Layers:

- The air layer is always present and uses `air_eps_r` and `air_sigma`.
- Earth layers are defined by parallel lists:
  `layer_thicknesses`, `layer_eps_r`, `layer_sigma`, `layer_mu_r`, `layer_sigma_m`.
- Each layer becomes a `GeoLayer` and is included in the `.poly` mesh regions.

Anomalies:

- Preferred interface: `anomalies=[BoxAnomaly(...), SphereAnomaly(...)]`.
- Legacy single-box interface is still supported using `anomaly_x`, `anomaly_y`, `anomaly_z`,
  and `anomaly_properties`.
- Anomalies are inserted into the `.poly` geometry and assigned region properties.

PML:

- The wrapper uses `PMLConfig` and writes `PML_present`, `PML_thickness`, and `PML_decay_type`
  into `elfe3D_input.txt`.
- `num_pml_layers` and `pml_layer_thickness` control the mesh padding around the domain.
- The current default decay type is `1`, corresponding to the exact reciprocal decay implementation.

