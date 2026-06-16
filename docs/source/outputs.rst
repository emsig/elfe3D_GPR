Outputs and Postprocessing
==========================

This page explains the file outputs that ``elfe3D_GPR`` writes and the tools available for postprocessing.

Solver Outputs
--------------
The solver writes electric and magnetic field results to the configured output files. These are typically stored in:

- ``out_<experiment_name>/electric_fields*``
- ``out_<experiment_name>/magnetic_fields*``

For ``elfe3D_GPR``, the output files contain columns with frequency, receiver coordinates and the field components for the recorded receiver lines. The text output files are commonly grouped in two ways:

- per frequency, using files of names with ``*_receiver_line`` appended to the base name,
- per receiver location, using files with the base names.

Ordering is always by increasing receiver number.

The file headers in the electric and magnetic field text files are as follows:

.. code-block:: text

   frequency coordinate_x coordinate_y coordinate_z  Ex_real Ex_imaginary  Ey_real Ey_imaginary  Ez_real Ez_imaginary

.. code-block:: text

   frequency coordinate_x coordinate_y coordinate_z  Hx_real Hx_imaginary  Hy_real Hy_imaginary  Hz_real Hz_imaginary

The ``elfe3D_GPR`` Python I/O module also writes ``.vtk`` field distributions for visualization. 
These ``.vtk`` files capture the computed electromagnetic field component over the mesh volume and can be opened 
in tools like ParaView or other VTK-compatible viewers. 
In the current workflow, ``.vtk`` files are usually generated during refinement or output 
if ``output_fields_vtk`` is enabled, and they typically store field data for the first frequency only.

Python Module Visualization
----------------------------

The repository includes helper modules in ``io/outputs`` for reading, processing, and visualizing results:

- ``io.outputs.fieldreader``: loaders for solver output and analytical reference data in different formats.
- ``io.outputs.visualize``: plotting classes for receiver-line comparison and error analysis
- ``io.outputs.postprocess``: additional processing utilities for field comparisons and statistics

These modules are intended to support notebook-based analysis of solver output.

Common output tools used in the example notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The example notebooks demonstrate the most common output workflow for ``elfe3D_GPR``:

- ``AnalyticalLoader`` loads semi-analytical reference solutions from CSV files.
- ``ElfeLoader`` reads solver-generated receiver-line output, especially ``electric_fields_receiver_line.txt``.
- ``ReceiverLinePlot`` compares multiple datasets on a 2x2 grid of amplitude, phase, real, and imaginary components.
- ``ReceiverLineErrorPlot`` plots normalized errors relative to a reference dataset in a 2x2 layout.
- ``ReceiverLineCombined`` shows one numerical dataset together with an analytical reference on the top row and the corresponding errors on the bottom row.
- ``ErrorHistogramPlot`` shows amplitude and phase error distributions for one or more datasets.

These tools are used throughout the notebooks, including ``examples/01_homogeneous_free-space.ipynb``, ``examples/02_homogeneous_earth.ipynb``, ``examples/03_two_layered_earth.ipynb``, and so on.