.. index:: Waiwera; installation

******************
Installing Waiwera
******************

Installing Waiwera on Linux
===========================

.. also needs various other software packages to be installed (compilers, build tools and other libraries needed by PETSc):
.. gcc g++ gfortran
.. git pkg-config valgrind ninja make meson cmake
.. openmpi | mpich blas lapack bison flex
.. can install these using preconfig script appropriate for your Linux distibution (or adapt nearest one if using another distribution). Needs root privileges. Or use Ansible?
.. substitute other C, C++ and Fortran compilers if desired
.. if installing on compute cluster, generally load appropriate modules if they are available, or build yourself if they are not, or pip install if possible (meson, ninja)

.. then run Waiwera configure script (distro-independent) in Waiwera root directory:

.. code-block:: bash

   python config.py

.. this has some optional parameters: --debug or --release, --no_rpath, --prefix, --libdir, --petsc_revision

.. this will detect PETSc on your system (via pkg-config) if it exists, otherwise will download it (check out from git repo) configure and build in external/PETSc directory.
.. will also detect FSON and Zofu (similarly via pkg-config) if they exists, otherwise will check out and configure as subprojects
.. creates build/ directory to do the build. When configure complete, build using:

.. code-block:: bash

   ninja -C build

.. can run unit tests (or put this in separate testing section?) using:

.. code-block:: bash

   python unit_tests.py

.. optional parameters: module names to test. Will run on 1 .. 4 parallel processes (or as many as there are, if that is less than 4).

.. install Waiwera using:

.. code-block:: bash

   ninja -C build install

.. will install Waiwera executable to prefix/bin. This should be on your $PATH if you want to be able to execute Waiwera from any directory.

.. uninstall using:

.. code-block:: bash

   ninja -C build uninstall

Installing Waiwera on Mac OSX
=============================

.. use Docker image

Installing Waiwera on MS Windows
================================

.. use Docker image
