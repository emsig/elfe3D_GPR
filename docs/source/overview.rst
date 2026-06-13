Overview
========

``elfe3D_GPR`` is a 3D edge-based **Finite Element** (FE) software for **Ground Penetrating Radar** (GPR) geophysical heterogeneous models in the frequency-domain. 
Unlike time-domain simulation codes, ``elfe3D_GPR`` produces field distributions directly in the frequency-domain based on the input frequencies of interest.

Origin
------

``elfe3D_GPR`` builds on the ``elfe3D`` software developed by Paula Rulff [ELFE3D]_, [RUL2023]_, [RUL2021]_ that simulates diffusive-field problems in Controlled-Source Electromagnetism. 
Both ``elfe3D`` and ``elfe3D_GPR`` are programmed in Fortran 90, that use ``tetgen`` [TETGEN]_, [SI2009]_ for producing unstructured tetrahedral mesh, and ``MUMPS`` [MUMPS]_ for solving the FE system of linear equations.

The major changes between ``elfe3D`` and ``elfe3D_GPR`` are as follows:

- The boundary value problem ``elfe3D`` solves is approximated for low-frequency diffusive fields. ``elfe3D_GPR`` now extends it to the full wave equation of electromagnetism. 
- A **Perfectly Matched Layer (PML)** is added to absorb outgoing waves from the truncated computational model domain. Specifically, it is a Uniaxial-PML that attenuates outgoing waves with an exact decay function [DIN2025]_, [FEN2019]_. 

These two major changes have been rigorously developed and described in the master's thesis written by Chaitanya Singh, and supervised by Paula Rulff and Evert Slob [SIN2025]_.

Since the thesis, ``elfe3D_GPR`` now also includes a Python I/O module for rapid and intuitive input model generation and output processing. You can read more here :doc:`inputs_and_models`.

Authors
-------

**Chaitanya Dinesh Singh**

- Implemented changes between ``elfe3D`` and ``elfe3D_GPR``.
- Email: chaitanya.singh@northumbria.ac.uk

Academic Supervisors
^^^^^^^^^^^^^^^^^^^^

- **Paula Rulff**: Author of ``elfe3D``.
- **Evert Slob** : Supervisor, author of MATLAB scripts for semi-analytical validation.

Contact
^^^^^^^

If you would like to discuss the theoretical background of ``elfe3D_GPR``, suggest changes to the software, or would like to contribute, 
you are welcome to contact Chaitanya at his email: chaitanya.singh@northumbria.ac.uk.

Citation
--------

If you publish results generated with ``elfe3D_GPR``, please give credit to the ``elfe3D_GPR`` developers 
by citing the thesis [SIN2025]_. 

Do not forget to acknowledge ``MUMPS`` [MUMPS]_ and ``tetgen`` [TETGEN]_ developers as well.

License
-------

``elfe3D_GPR`` is licensed under the Apache License, Version 2.0 (the 
\"License\"); you may not use any ``elfe3D_GPR`` files except in compliance 
with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software 
distributed under the License is distributed on an "AS IS" BASIS, 
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
See the License for the specific language governing permissions and 
limitations under the License.
