.. index:: Waiwera; running, running; Waiwera

***************
Running Waiwera
***************

How Waiwera is executed depends on whether it is being run as a Docker container, or as a natively built executable (Linux only).

.. index:: running; Docker, Docker; running, PyWaiwera; waiwera-dkr
.. _run_docker:

Running Waiwera using Docker
============================

The easiest way to run Waiwera via Docker is by using the script ``waiwera-dkr``, which is part of the `PyWaiwera <https://pypi.org/project/pywaiwera>`_ Python library.  Besides PyWaiwera, You will need `Python <https://www.python.org/>`_ and `Docker <https://www.docker.com/>`_ installed on your machine to be able to use this script. For more details, see :ref:`using_docker`.

What the ``waiwera-dkr`` script does
---------------------------------------

This script does three main things:

- checks if the Waiwera Docker container image has been installed, installs it if necessary, or updates it if a newer version is available
- runs Waiwera inside the Docker container
- manages sharing of files between the Docker container and the directory in which you run Waiwera

How to run the ``waiwera-dkr`` script
----------------------------------------

The script is run from the command line in the same way as any other program. The name of the JSON input file (see :ref:`waiwera_input`) for your simulation is specified as an argument, e.g. if the simulation has the filename ``model.json``, you can run it as follows:

.. code-block:: bash

   waiwera-dkr model.json

This would run the simulation in serial. Running in serial is only suitable for small problems. Waiwera is designed primarily for large problems that need to be run in parallel.

.. index:: running; number of processes

To run Waiwera using Docker in parallel, the number of parallel processes must be specified using the ``-np`` parameter, e.g.:

.. code-block:: bash

   waiwera-dkr -np 16 model.json

runs Waiwera in parallel with 16 processes.

.. index:: Docker; options

Optional parameters for the ``waiwera-dkr`` script
-----------------------------------------------------

Besides the ``-np`` option for specifying the number of processes, the ``waiwera-dkr`` script has some other optional parameters for controlling its behaviour. Details of all available options can be displayed using the ``--help`` (or ``-h``) option, e.g.:

.. code-block:: bash

   waiwera-dkr --help

or:

.. code-block:: bash

   waiwera-dkr -h

These options include:

- ``--noupdate`` (or ``-nu``): do not check for or download an updated Waiwera Docker image before running (the default behaviour is to check before each run, and download an updated image if there is one available)
- ``--update`` (or ``-u``): check for an updated Waiwera Docker image and download if available, and then exit (without running anything)
- ``--test_volume`` (or ``-tv``): test that the sharing of files between the current directory and the Docker container is working correctly, and then exit
- ``--examples`` (or ``-e``): create example models (files will be written into an ``examples`` sub-directory) and then exit (note any existing files will be overwritten).s
- ``--verbose`` (or ``-v``): output additional diagnostic message while running (for debugging Docker-related problems)
- ``--interactive`` (or ``-it``): start an interactive Linux terminal inside the Docker container. If a command is also specified then this will be run.

.. index:: Docker; file paths

File paths when running with Docker
-----------------------------------

The Waiwera JSON input file (see :ref:`waiwera_input`) contains some paths to other files, e.g. the mesh file (see :ref:`specifying_mesh`). There are two things to note about file paths when running using Docker:

- file paths must always be specified using POSIX (i.e. Linux-style) file path syntax, i.e. forward slashes for directory delimiters (not backslashes as on Windows), and any spaces in the file path (usually better avoided if possible) should be "escaped" by preceding them with backslashes. This is because Waiwera is run using Linux inside the Docker container.
- any files specified in the JSON input file name need to be in the directory that Waiwera is being run in, or a subdirectory of it. This is because those are the only directories that are shared with the Docker container.

The same considerations apply when running Waiwera using the ``waiwera-dkr`` script and specifying a path to the simulation input file on the command line. In general, when running with Docker it is recommmended to run from the directory containing the simulation input file. Then avoids the need to specify a path to your file, and simplifies the directories that need to be shared with Docker.

.. index:: Docker; Python

Running Waiwera via Docker from a Python script
-----------------------------------------------

It is also possible to use PyWaiwera to run Waiwera via Docker from
within a Python script. This is done by importing the ``pywaiwera``
package, creating a Docker environment, and using that to run the
Waiwera simulation, as in the following example:

.. code-block:: python

   import pywaiwera

   env = pywaiwera.docker.DockerEnv()
   env.run_waiwera('model.json')

.. index:: running; native Linux executable
.. _run_native:

Running the native Linux Waiwera executable
===========================================

The native Linux Waiwera executable (see :ref:`native_linux_build`) can be run either in serial or in parallel. To run in serial from the command line simply type ``waiwera`` followed by the name of the Waiwera JSON input file (see :ref:`waiwera_input`).

For example, if the JSON input file for the simulation has the filename ``model.json``, it can be run as follows:

.. code-block:: bash

   waiwera model.json

Running in serial is only suitable for small problems. Waiwera is designed primarily for large problems that need to be run in parallel.

Waiwera is run in parallel using MPI (Message Passing Interface), via the ``mpiexec`` command. As in serial, the ``waiwera`` command is followed by the JSON input filename, but all of this is preceded by the ``mpiexec`` command and the ``-np`` parameter to specify the number of processes the simulation is to be run on. For example:

.. code-block:: bash

   mpiexec -np 16 waiwera model.json

runs Waiwera on the input file ``model.json``, in parallel on 16 processes.

Running Waiwera with the ``--help`` option (or no arguments at all) instead of a filename prints a help message on basic usage to the console. The ``-h`` option does the same thing, but also prints out PETSc-related help.

Similarly, running Waiwera with the ``--version`` option instead of a filename prints the Waiwera version to the console, and the ``-v`` option does the same thing but also prints PETSc version information.

.. index:: running; number of processes

Choosing the number of parallel processes
=========================================

The appropriate number of processes depends on how many are available, and on the size of the problem. Waiwera will generally run faster when more processes are used, but there is usually a point at which there are diminishing returns from adding more processes. In more extreme cases, using more processes can even start to slow the simulation down.

This is because the processes need to communicate with each other, e.g. to communicate values at cells on the edges of the mesh partitions (see :ref:`mesh_partitioning`). There is a cost involved with this communication, which rises as the number of processes is increased. Eventually, if too many processes are used, the communication costs start to outweigh the benefits of increased parallelisation.

These considerations apply to most MPI programs. The `PETSc documentation <https://www.mcs.anl.gov/petsc/petsc-dev/docs/faq.html>`_ recommends that there should be an absolute minimum of 10,000 unknowns per process for good parallel performance, with at least 20,000 unknowns per process being preferable.

For Waiwera the number of unknowns per process is equal to the number of cells on each process multiplied by the number of unknowns per cell. The cells are usually divided approximately evenly between the processes, so the number of cells per process is approximately the total number of cells divided by the number of processes. The number of unknowns per cell depends on the :ref:`eos` module being used.

For example, the :ref:`water_air_energy_eos` EOS has three unknowns per cell. Supposing the mesh has 100,000 cells, this means the total number of unknowns in the problem is 300,000. Hence the maximum number of processes that should be used is approximately 30, with around 15 being preferable.

.. index:: running; parameters
.. _petsc_command_line_parameters:

PETSc command line parameters
=============================

When Waiwera is run, the main parameter it takes is the filename, which should follow the ``waiwera`` command (or ``waiwera-dkr`` if :ref:`run_docker`). However, it is also possible to control many PETSc-related aspects of the simulation by adding other command line parameters, which can be specified after the filename.

These PETSc command line parameters can be used, for example, to control the behaviour of the PETSc linear and non-linear solvers used by Waiwera, as well as many other options such as diagnostic or debugging output. Some of these options (e.g. the linear and non-linear solver parameters) can also be specified in the Waiwera JSON input file.

For example, if running a native Linux executable:

.. code-block:: bash

   mpiexec -np 16 waiwera model.json -log_view

again runs Waiwera in parallel on 16 processes, but also displays PETSc profiling information at the end of the run (data on how much time is spent in various parts of the code, etc.).

..
   TODO: how to do this when running with Docker

.. code-block:: bash

   waiwera-dkr -np 16 model.json -log_view

More information about specific PETSc command line parameters can be found in the `PETSc <https://www.mcs.anl.gov/petsc/>`_ documentation.

Run-time console output
=======================

As Waiwera runs, by default all the log messages (see :ref:`setup_logfile`) being written to the YAML log file are also echoed to the console output, so the progress of the simulation can be monitored. If this is not needed it can be disabled in the input JSON file (see :ref:`control_log_output`).
