Workflow
========

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

Documentation syntax examples
-----------------------------

URLs
^^^^

Use standard Markdown links:

.. code-block:: markdown

   [Read the Docs](https://readthedocs.org)


Images
^^^^^^

Use Markdown image syntax or the Sphinx directive.

.. code-block:: markdown

   ![Model diagram](../_static/model-diagram.png)


.. code-block:: rst

   .. image:: ../_static/model-diagram.png
      :alt: Model diagram
      :align: center


Tables
^^^^^^

A simple Markdown table works well:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Field
     - Description
   * - ---
     - ---
   * - `air_eps_r`
     - Relative permittivity for the air layer
   * - `layer_sigma`
     - Electrical conductivity for each earth layer


LaTeX and math
^^^^^^^^^^^^^^

Inline math:

.. code-block:: markdown

   The wave number is $k = \omega \sqrt{\mu \epsilon}$.


Display math:

.. code-block:: markdown

   $$
   \nabla \times \mathbf{E} = -j \omega \mathbf{B}
   $$


If the parser supports MyST and Sphinx math, this will render as math in RTD.

Cross-references
^^^^^^^^^^^^^^^^

Use explicit labels and `:ref:` references when you want robust links between pages.

For example:

- See :ref:`python_interface` for the Python package API.
- See :doc:`quickstart` for a minimal run-through.
