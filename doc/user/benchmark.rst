.. index:: testing; benchmark

*******************************
Running Waiwera benchmark tests
*******************************

Benchmark tests
===============

The ``/test/benchmark`` directory of the Waiwera source code contains
a suite of benchmark test problems which can be used to verify that
the code is working correctly. They can also serve as examples to help
users become familiar with how Waiwera simulations are set up.

The Waiwera results for each benchmark test are compared with results
generated using TOUGH2. For some tests, analytical or semi-analytical
solutions are also available, or results using other simulators have
previously been published, and these can also be used for comparison
purposes.

The benchmark tests are organised into several sub-directories under
the ``/test/benchmark`` directory, one for each category of
test. These include:

- **MINC**: benchmark tests with fractured rocks represented using the
  Multiple INTeracting Continua (MINC) approach (see :ref:`minc`)
- **Model Intercomparison Study**: a selection of test problems from
  the 1980 "Geothermal Model Intercomparison Study"
- **NCG**: tests for non-condensible gas equations of state (see
  :ref:`water_ncg_eos`)
- **sources**: tests for verifying correct behaviour of various source
  controls (see :ref:`source_controls`)

More details about some of these test problems can be found in various
`publications <https://waiwera.github.io/pubs/>`_ about Waiwera code
development.

.. index:: CREDO

CREDO
=====

The `CREDO <https://github.com/waiwera/credo>`_ library is used to
manage the Waiwera benchmark tests. CREDO is a Python toolkit for
running, analysing and benchmarking computational models, originally
developed as part of the `Underworld
<https://www.underworldcode.org/>`_ geodynamics code project.

Waiwera uses a fork of the original CREDO library, modified to include
modules specific to running and analysing Waiwera and TOUGH2 models.

CREDO can be installed from the Python Package Index (PyPI) using
``pip``, e.g.:

.. code-block:: bash

   pip install waiwera-credo

PyTOUGH
=======

Running the Waiwera benchmark tests also requires `PyTOUGH
<https://github.com/acroucher/PyTOUGH>`_, a Python library for
handling TOUGH2 simulations. This is needed for extracting the TOUGH2
results for each benchmark test, so they can be compared with the
Waiwera results.

See the `PyTOUGH website <https://github.com/acroucher/PyTOUGH>`_ for
instructions on how to install PyTOUGH.

Running the test suite
======================

The whole suite of benchmark tests can be run using the Python script
``benchmark_tests.py``, in the root directory of the Waiwera source
code:

.. code-block:: bash

   python benchmark_tests.py

By default, all tests will be run in serial. To run them in parallel,
the script takes a ``-np`` parameter which can be used to specify the
number of parallel processes, e.g.:

.. code-block:: bash

   python benchmark_tests.py -np 4

would run the benchmark test suite using four processes.

Running using Docker
====================

The above commands apply if running a native Linux executable (see
:ref:`native_linux_build`). If you are running Waiwera via Docker (see
:ref:`using_docker`), you can use the ``benchmark_tests.py`` script
using the ``--docker`` (or ``-d``) parameter, e.g.:

.. code-block:: bash

   python benchmark_tests.py --docker -np 2

This would run all the benchmark tests via Docker, using the
``waiwera-dkr`` script (see :ref:`run_docker`), and using two parallel
processes.

Note that if you are running Waiwera via Docker, using the
``waiwera-dkr`` script provided with PyWaiwera, you will probably not
have a copy of the Waiwera source code on your machine. In that case,
to run the benchmark tests you will first have to clone or download
the Waiwera source code from its `Github repository
<https://github.com/waiwera/waiwera>`_.

Test results
============

The result (pass or fail) for each benchmark test are written to the
console, together with the overall result for the whole suite.

In addition, an HTML page with similar information is written to
``/test/benchmark/test_summary.html``. This page can be viewed in a web
browser, and contains links to detailed results for each benchmark
test, including plots showing selected model results compared with
results from other simulators (and/or analytical solutions).

Running individual tests
========================

It is also possible to run individual benchmark tests, by navigating
to the desired test directory and running the Python script for that
test.

The test scripts all have filenames of the form ``test_*.py``. Like
the ``benchmark_tests.py`` script used for running the whole test
suite, the individual test scripts also take the same ``-np`` and
``--docker`` (or ``-d``) parameters for specifying the number of
parallel processes and whether Docker is used to run the Waiwera
models.

For example, in the ``/test/benchmark/minc/column`` directory, running:

.. code-block:: bash

   python test_minc_column.py -np 2

runs the MINC column benchmark test on two processes.

The test results are summarised on the console. The HTML page with
more detailed information about the test can be found in the test
directory at ``/output/TESTNAME/TESTNAME-report.html``, where
``TESTNAME`` is the name of the individual test.
