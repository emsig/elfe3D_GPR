Python API Reference
====================

This page provides a curated reference for the Python modules that are currently used to build and run `elfe3D_GPR` simulations.

Current strategy
----------------

The `io/` folder is now packaged for installation, and the supported import namespace is `elfe3d_gpr_io`.
This page continues to provide a curated reference for the core user-facing entry points.

Key Python entry points
-----------------------

- `elfe3d_gpr_io.runner.ProjectPaths`
- `elfe3d_gpr_io.runner.run_tetgen`
- `elfe3d_gpr_io.runner.run_solver`
- `elfe3d_gpr_io.inputs.survey.IOConfig`
- `elfe3d_gpr_io.inputs.survey.GPRSurvey.build`

Recommended reference approach
------------------------------

Once the Python I/O layer is packaged, it will be possible to add automatic `autodoc` pages for each module. Until then, keep this reference hand-curated and aligned with the example notebooks.

Notes for future automation
---------------------------

The easiest automation path is to make `io/` importable from the repo root or to add package metadata such as `pyproject.toml`. That will allow Sphinx to import the modules and publish `autodoc`-generated documentation for the current API.

Python I/O Module
-----------------

The Python I/O wrapper is the primary user-facing entry point for notebooks and scripted workflows.

Main user modules
-----------------

- `elfe3d_gpr.runner`
  - `ProjectPaths` — path configuration for the solver executable and working directories
  - `run_tetgen()` — launch TetGen from Python
  - `run_solver()` — run the compiled `elfe3d_gpr` executable

- `elfe3d_gpr.inputs.survey`
  - `IOConfig` — manages input and output folder locations
  - `GPRSurvey` — overall survey definition and builder

How to use the wrapper
----------------------

Install the package in editable mode from the repository root:

.. code-block:: bash

   pip install -e .


Then import from the supported package namespace (the package is installed as `elfe3d_gpr_io`):

.. code-block:: python

   from elfe3d_gpr_io.runner import ProjectPaths, run_tetgen, run_solver
   from elfe3d_gpr_io.inputs.survey import GPRSurvey


Recommended usage example
^^^^^^^^^^^^^^^^^^^^^^^^^

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
       # Prefer `f_list=[...]` for explicit frequency control; legacy `ricker_central_f` is deprecated.
       f_list=[100e6],
       antenna_position=[0.0, 0.0, 0.025],
       num_receivers_inline=48,
   )
   
   survey.generate()
   run_tetgen(paths, survey.io.poly_file)
   run_solver(paths, survey)


The first example notebook, `examples/01_homogeneous_free-space.ipynb`, shows the recommended package-based workflow.
Legacy imports from `io` may still work when the repository root is on `PYTHONPATH`, but the supported package namespace is `elfe3d_gpr_io`.

Notes on page structure
-----------------------

This page focuses on the Python I/O wrapper usage and import conventions.
The `workflow` page explains the repository architecture and the end-to-end process from model definition to solver execution.
Keeping them separate is helpful because readers who want code examples can go directly to `python_interface`, while readers who want the overall process can read `workflow`.

Cross-references
----------------

Use `:ref:` labels when linking between pages in MyST/Sphinx.
For example:

- See :ref:`workflow` for the high-level repository workflow.
- See :doc:`quickstart` for a minimal step-by-step run.

If you prefer a simple link, use Markdown syntax:

.. code-block:: markdown

   [Quick Start](quickstart)


But for RTD, `:ref:` and `:doc:` produce more robust links when the target page is in the same documentation tree.

Workflow
--------

This page explains how the repository pieces fit together and how data flows from model definition to solver results.

Repository layout
-----------------

- `elfe3D_GPR/` — the Fortran source tree and solver executable
- `io/` — Python helpers for model generation, mesh preparation, solver execution, and result handling
- `examples/` — notebook-based example workflows and demonstration cases
- `tests/` — automated tests and validation scripts
- `validation/` — reference results and validation comparisons

Typical workflow
----------------

1. Define the GPR survey in Python using `elfe3d_gpr_io.inputs.survey.GPRSurvey.build()`.
2. Generate the mesh input and field input files with `survey.generate()`.
3. Run TetGen on the generated `.poly` file using `elfe3d_gpr_io.runner.run_tetgen()`.
4. Execute the Fortran solver with `elfe3d_gpr_io.runner.run_solver()`.
5. Read and visualize the solver output from `out_<experiment_name>`.

Example workflow
----------------

.. code-block:: python

   from pathlib import Path
   from elfe3d_gpr_io.runner import ProjectPaths, run_tetgen, run_solver
   from elfe3d_gpr_io.inputs.survey import GPRSurvey
   
   MASTER_PATH = Path('..') / 'elfe3D_GPR'
   MASTER_PATH = MASTER_PATH.resolve()
   
   paths = ProjectPaths(
       master_dir=MASTER_PATH,
       exec_rel='elfe3d_gpr',
       use_wsl=True,
   )
   
   survey = GPRSurvey.build(
       experiment_name='air_only',
       base_dir=MASTER_PATH,
       air_eps_r=1.0,
       air_sigma=1e-16,
       layer_thicknesses=[0.3],
       layer_eps_r=[1.0],
       layer_sigma=[1e-16],
       layer_mu_r=[1.0],
       layer_sigma_m=[0.0],
       # Note: use explicit frequency lists; legacy `ricker_central_f` behaviour is deprecated.
       f_list=[100e6],
       antenna_position=[0.0, 0.0, 0.025],
       num_receivers_inline=48,
   )
   
   survey.generate()
   run_tetgen(paths, survey.io.poly_file)
   run_solver(paths, survey)


Key components
--------------

- `elfe3d_gpr_io.runner` contains the execution layer and host path handling.
- `elfe3d_gpr_io.inputs.survey` is the central model builder.
- `elfe3d_gpr_io.inputs.writeinputfiles` writes the solver input files.
- `elfe3d_gpr_io.inputs.writetetgenpoly` writes the TetGen geometry input.