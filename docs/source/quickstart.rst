Quick Start
===========

This page describes two minimal ways to start using `elfe3D_GPR`.

Option A: native Fortran run
----------------------------

The simplest path is to run the compiled Fortran solver directly with the example input files provided in the repository.
This is the default `elfe3D`-style workflow and is useful for validating the solver executable before using the Python wrapper.

Steps:

1. Build the Fortran solver in `elfe3D_GPR/`.
2. Change to the `elfe3D_GPR/` directory.
3. Run the solver executable, for example:

.. code-block:: bash

   cd elfe3D_GPR
   ./elfe3d_gpr


4. Inspect the output in the matching `out_*` folder.

Option B: notebook workflow with the Python I/O wrapper
-------------------------------------------------------

The recommended workflow for new work is to use the Python I/O wrapper and the example notebook in `examples/01_homogeneous_free-space.ipynb`.
This approach uses the installed package namespace `elfe3d_gpr_io` and demonstrates the full model build → TetGen → solver → output cycle.

1. Install or activate your Python environment from the repository root.
2. Install the Python wrapper:

.. code-block:: bash

   pip install -e .


3. Open `examples/01_homogeneous_free-space.ipynb` in JupyterLab or Jupyter Notebook.
4. Run the notebook to generate solver input files, launch TetGen, execute the solver, and inspect results.

Minimal Python example
----------------------

.. code-block:: python

   from pathlib import Path
   from elfe3d_gpr_io.runner import ProjectPaths, run_tetgen, run_solver
   from elfe3d_gpr_io.inputs.survey import GPRSurvey
   
   master_dir = Path('..') / 'elfe3D_GPR'
   master_dir = master_dir.resolve()
   
   paths = ProjectPaths(
       master_dir=master_dir,
       exec_rel='elfe3d_gpr',
       use_wsl=True,
   )
   
   survey = GPRSurvey.build(
       experiment_name='air_only',
       base_dir=master_dir,
       air_eps_r=1.0,
       air_sigma=1e-16,
       layer_thicknesses=[0.3],
       layer_eps_r=[1.0],
       layer_sigma=[1e-16],
       layer_mu_r=[1.0],
       layer_sigma_m=[0.0],
       # Note: prefer explicit frequency lists. Legacy `ricker_central_f` is deprecated.
       f_list=[100e6],
       antenna_position=[0.0, 0.0, 0.025],
       num_receivers_inline=5,
   )
   
   survey.generate()
   run_tetgen(paths, survey.io.poly_file)
   run_solver(paths, survey)


Linking to other docs
---------------------

- See :ref:`workflow` for the high-level repository and data flow.
- See :ref:`python_interface` for the Python I/O wrapper usage and supported import conventions.

Notes on the Python wrapper
---------------------------

The Python I/O wrapper is described in `python_interface`.
This page is focused on the minimal user workflow, while `python_interface` is focused on the supported package namespace and interface patterns.
