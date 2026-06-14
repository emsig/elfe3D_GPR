Quick Start
===========

This page describes the two ways to start using ``elfe3D_GPR``.

.. _fortran_quick:

Option A: Fortran-only Approach
-------------------------------

After following the Fortran installation instructions in :doc:`installation`, 
you can run your own ``elfe3D_GPR`` simulation directly. These few steps 
highlight a brief glance of what the core workflow is:

1. Going to the ``elfe3D_GPR/`` directory (the directory with the executable).
2. Creating/reviewing input files for the solver within the ``elfe3D_GPR\in`` folder (more on input files in :doc:`inputs_and_models`).
3. Run ``tetgen`` on the ``.poly`` file with preferred quality settings. This step creates the necessary unstructured finite element mesh files that ``elfe3D_GPR`` needs to solve the wave equation using the finite element method. It also includes a ``.vtk`` file that can be used for visualizing the mesh including the PML.
4. Run the ``elfe3D_GPR`` solver executable using the command:

.. code-block:: bash

   ./elfe3d_gpr

Once the simulation completes, you will have outputs in two directories:
1.  ``elfe3D_GPR\out_{name_of_experiment}`` with text files that record electric and magnetic field components per input frequency and receiver coordinates.
2.  ``elfe3D_GPR\in`` has an updated ``.vtk`` file with electric and magnetic field distributions in the entire volume. 

.. warning::

   This `.vtk` file currently only stores field data for the first frequency of simulation. 
   This is to prevent the ``.vtk`` files from growing too large. 
   A feature to allow writing fields for all frequencies of simulation would be provided in a future version.


Option B: Python I/O Module-Assisted Approach
---------------------------------------------

The easier approach is to use the Python I/O module ``elfe3d_gpr_io``. 
The example notebook in `examples/01_homogeneous_free-space.ipynb` 
gives a good starting overview of the complete simulation workflow, going from input definitions, path configurations, 
running the Fortran executable from the notebook, and visualization of outputs.

Minimal Python example
^^^^^^^^^^^^^^^^^^^^^^

Otherwise, after following installation instructions for ``elfe3d_gpr_io``, you can simply 
run the following minimal python script:

.. code-block:: python

   from pathlib import Path
   from elfe3d_gpr_io.runner import ProjectPaths, run_tetgen, run_solver
   from elfe3d_gpr_io.inputs.survey import GPRSurvey

   MASTER_PATH = (Path("..") / "elfe3D_GPR").resolve() # Base path variable, optional

   paths = ProjectPaths(
      master_dir = MASTER_PATH,   # Path where the elfe3D_GPR executable is located. Currently it is the same as repository root.
      exec_rel   = "",            # Relative path to find the executable. With the default installation it exists in the repository root, hence an empty string.
      use_wsl    = False,         # True if running the notebook from Windows with WSL installed.
   )

   # Defining elfe3D_GPR simulation inputs
   f     = 100e6     # frequency of simulation in Hertz
   wave  = 3e8 / f   # approximate length unit useful only to define simulation domain
   
   survey = GPRSurvey.build(
      experiment_name='air_only', # Name of the experiment, will also be used to create file and folder names for I/O.
      base_dir=MASTER_PATH,

      # Domain extents [in meters], along x and y (lateral) axes, and z (height) axis.
      x_e = [-wave/10, 1 + wave/10],  
      y_e = [-wave/10,     wave/10],
      z_e = [-wave/10,     wave/10],

      # Simple Homogeneous Medium Simulation: Air-only
      layer_thicknesses=[0.3],
      layer_eps_r=[1.0],
      layer_sigma=[1e-16],
      layer_mu_r=[1.0],
      layer_sigma_m=[0.0],

      # list of frequencies of simulation
      f_list=[100e6],                     # list of frequency
      antenna_position=[0.0, 0.0, 0.025], # [x, y, z], 25 mm above the surface.
      num_receivers_inline=5,             # distributes receivers automatically along a radial line from the source to the end of domain
   )
   
   survey.generate() # Creates the input files needed by the Fortran core solver.

   run_tetgen(paths, survey.io.poly_file) # Runs the tetgen meshing algorithm right from the notebook/python script.
   run_solver(paths, survey)              # As the input and mesh files are now generated, the Fortran solver is now called to run the simulation.


Running this script will result in the same outputs as the Fortran option described in the section above :ref:`fortran_quick`. 
You can now proceed with more details on the inputs of ``elfe3D_GPR`` modelling described in :doc:`inputs_and_models`.