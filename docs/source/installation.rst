Installation
============

This page describes the core installation and build steps for ``elfe3D_GPR``.

System prerequisites
--------------------

- a modern Fortran compiler (Fortran 2008 or later compiler for ``.f90`` sources; tested on ``gfortran``).
- ``make``
- ``OpenBLAS`` (with ``libopenblas-dev``, for example).
- ``TetGen`` for mesh generation.
- ``MUMPS`` for the direct linear algebra solver.
- ``Python 3.10+`` for the Python I/O wrapper and docs build.

Build the Fortran solver
------------------------

The core solver source is located in ``elfe3D_GPR/``.
``elfe3D_GPR`` is written in modern Fortran and uses shared-memory parallelisation with OpenMP.
The system of equations is solved with a direct solver.

Install TetGen
^^^^^^^^^^^^^^

``TetGen`` can be downloaded from <https://wias-berlin.de/software/index.jsp?id=TetGen>. 
Otherwise, on Debian/Ubuntu, you can:

.. code-block:: bash

   sudo apt install tetgen


### Install MUMPS

``MUMPS`` is available at <https://mumps-solver.org>.
Build it from source and copy the required headers into the ``elfe3D_GPR/`` source tree.

.. code-block:: bash

   wget https://mumps-solver.org/MUMPS_5.7.3.tar.gz
   tar zxvf MUMPS_5.7.3.tar.gz
   cd MUMPS_5.7.3
   cp Make.inc/Makefile.debian.SEQ Makefile.inc
   sudo apt install libmetis-dev libparmetis-dev libscotch-dev libptscotch-dev libatlas-base-dev openmpi-bin libopenmpi-dev liblapack-dev libscalapack-openmpi-dev
   make all


Copy the required header files into the Fortran source directory:

.. code-block:: bash

   cp MUMPS/libseq/mpif.h elfe3D_GPR/.
   cp MUMPS/include/zmumps_root.h elfe3D_GPR/.
   cp MUMPS/include/zmumps_struc.h elfe3D_GPR/.


Update the ``MUMPS_LIB_DIR`` variable in ``elfe3D_GPR/Makefile`` to point to your installed MUMPS library directory:

.. code-block:: make

   MUMPS_LIB_DIR = /path/to/your/MUMPS_5.7.3/lib


Compile the solver:

.. code-block:: bash

   cd elfe3D_GPR
   make all


Set the OpenMP thread count to avoid oversubscription:

.. code-block:: bash

   export OMP_NUM_THREADS=<number_of_threads>


You can run the Fortran solver independent of the Python I/O module. Once the ``make`` process 
succeeds, you can simply:

.. code-block:: bash

   ./elfe3d


Currently, it will solve for the homogeneous air model (the simplest reference example) as the input 
files have been made available on release for verification of the fortran installation.


Verify Python setup
-------------------

Check the active Python environment and package tools:

.. code-block:: bash

   python --version
   python -m pip --version
   python -m pip install --upgrade pip setuptools wheel


Install the Python I/O module
-----------------------------

The Python I/O module is implemented in the `io/` folder and packaged under the namespace `elfe3d_gpr_io`.
Install it from the repository root:

.. code-block:: bash

   cd .. #to go up one directory from elfe3D_GPR into the repository root
   pip install -e .


The supported import namespace is used as:

.. code-block:: python

   from elfe3d_gpr_io.runner import ProjectPaths, run_tetgen, run_solver
   from elfe3d_gpr_io.inputs.survey import GPRSurvey

