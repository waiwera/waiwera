.. index:: Waiwera; running, running; Waiwera

***************
Running Waiwera
***************

Like all programs based on the `PETSc <https://www.mcs.anl.gov/petsc/>`_ library, Waiwera can be run either in serial or in parallel.

.. index:: running; serial

Running in serial
=================

Waiwera can be run in serial from the command line simply by typing ``waiwera`` followed by the name of the Waiwera JSON input file (see :ref:`waiwera_input`).

For example, if the JSON input file for the simulation has the filename ``model.json``, it can be run as follows:

.. code-block:: bash

   waiwera model.json

Running in serial is only suitable for small problems. Waiwera is primarily designed for large problems that need to be run in parallel (see :ref:`run_parallel`).

.. index:: running; parallel
.. _run_parallel:

Running in parallel
===================

Waiwera can be run in parallel using MPI (Message Passing Interface). Clearly, this requires that MPI is installed on the machine Waiwera is to be run on. There are various implementations of MPI available (e.g. `OpenMPI <https://www.open-mpi.org/>`_ or `MPICH <https://www.mpich.org/>`_), but which one is installed should not matter.

Waiwera is run in parallel using the MPI ``mpiexec`` command. As in serial, the ``waiwera`` command is followed by the JSON input filename, but all of this is preceded by the ``mpiexec`` command and the ``-np`` flag to specify the number of processes the simulation is to be run on. For example:

.. code-block:: bash

   mpiexec -np 16 waiwera model.json

runs Waiwera on the input file ``model.json``, in parallel on 16 processes.

.. index:: running; number of processes

Choosing the number of processes
--------------------------------

The appropriate number of processes depends on how many are available, and on the size of the problem. Waiwera will generally run faster when more processes are used, but there is usually a point at which there are diminishing returns from adding more processes. In more extreme cases, using more processes can even start to slow the simulation down.

This is because the processes need to communicate with each other, e.g. to communicate values at cells on the edges of the mesh partitions (see :ref:`mesh_partitioning`). There is a cost involved with this communication, which rises as the number of processes is increased. Eventually, if too many processes are used, the communication costs start to outweigh the benefits of increased parallelisation.

These considerations apply to most MPI programs. The PETSc documentation recommends that there should be an absolute minimum of 10,000 unknowns per process for good parallel performance, with at least 20,000 unknowns per process being preferable.

For Waiwera the number of unknowns per process is equal to the number of cells on each process multiplied by the number of unknowns per cell. The cells are usually divided approximately evenly between the processes, so the number of cells per process is approximately the total number of cells divided by the number of processes. The number of unknowns per cell depends on the :ref:`eos` module being used.

For example, the :ref:`water_air_energy_eos` EOS has three unknowns per cell. Supposing the mesh has 100,000 cells, this means the total number of unknowns in the problem is 300,000. Hence the maximum number of processes that should be used is approximately 30, with around 15 being preferable.

.. index:: running; parameters
.. _petsc_command_line_parameters:

PETSc command line parameters
=============================

When Waiwera is run, the main parameter it takes is the filename, which should follow the ``waiwera`` command. However, it is also possible to control many PETSc-related aspects of the simulation by adding other command line parameters, which can be specified after the filename.

These PETSc command line parameters can be used, for example, to control the behaviour of the PETSc linear and non-linear solvers used by Waiwera, as well as many other options such as diagnostic or debugging output. Some of these options (e.g. the linear and non-linear solver parameters) can also be specified in the Waiwera JSON input file.

More information about specific PETSc command line parameters can be found in the `PETSc <https://www.mcs.anl.gov/petsc/>`_ documentation.
