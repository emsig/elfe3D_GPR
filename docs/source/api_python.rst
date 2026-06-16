Python I/O Module Reference
===========================

This page provides a reference for the Python modules in the ``io/`` folder that 
can be optionally used to build and run ``elfe3D_GPR`` simulations.

.. warning::

   Note that the Python I/O module is under active development due to the many possible features 
   that are useful with ``elfe3D_GPR`` modelling that will be implemented in the future.
   As such, the following page is incomplete.

The I/O Module
--------------

The ``io/`` folder is packaged for installation, and the supported import namespace is ``elfe3d_gpr_io``. 

.. note::
   Due to the dependency of ``elfe3D_GPR`` on ``MUMPS``, it is likely that converting ``elfe3D_GPR`` to a 
   Python-only API will prove difficult. This will be investigated.

Main User Modules
-----------------

These are the current major user modules for python scripting of input files for ``tetgen`` and ``elfe3D_GPR``, running them in a 
slightly ad-hoc way in a jupyter notebook, and outputs handling.

- ``elfe3d_gpr.runner``:
   - ``ProjectPaths``: path configuration for the solver executable and working directories
   - ``run_tetgen()``: launch TetGen from Python. Supports Linux and Windows-WSL installs.
   - ``run_solver()``: run the compiled ``elfe3d_gpr`` executable. Supports Linux and Windows-WSL installs.
- ``elfe3d_gpr.inputs.survey``:
   - ``GPRSurvey``: Overall GPR simulation definition class to generate input files for the Fortran executables of ``tetgen`` and ``elfe3D_GPR``.
- ``elfe3d_gpr.outputs.fieldreader``, ``elfe3d_gpr.outputs.postprocess``, ``elfe3d_gpr.outputs.visualize`` for the output handling including reference data.