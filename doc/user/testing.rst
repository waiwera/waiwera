.. index:: Waiwera; testing, testing

***************
Testing Waiwera
***************

.. index:: testing; unit tests

Unit tests
----------

Unit tests are lower-level tests of individual subroutines in the code. These are used mostly during code development, but may also be run to help verify that Waiwera has been installed correctly.

Unit tests are carried out using the `FRUIT <https://sourceforge.net/projects/fortranxunit/>`_ framework for Fortran unit testing, via a separate Python interface called `FRUITPy <https://github.com/acroucher/FRUITPy>`_. You will need to install both FRUIT and FRUITPy before running the unit tests.

.. index:: FRUIT, FRUIT; installation, FRUITPy, FRUITPy; installation

Installing FRUIT and FRUITPy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. todo:: add instructions on installing FRUIT/FRUITPy

Running unit tests
^^^^^^^^^^^^^^^^^^

To run all unit tests, navigate to the main Waiwera directory and type:

.. code-block:: bash

  python unit_tests.py

It is also possible to run unit tests only for specific code modules, by listing the names of those code modules after the command above. For example, to test only the `IAPWS` module, type:

.. code-block:: bash

  python unit_tests.py IAPWS

or to test only the `face` and `cell` modules, type:

.. code-block:: bash

  python unit_tests.py face cell

If the tests have successfully passed, the output from FRUIT will appear something like this:

.. code-block:: bash

  All tests passed.
  Hit rate:
  asserts:  670 / 670 (100%)
  cases  :  86 / 86 (100%)

The precise numbers of asserts and cases will vary, depending on how many modules are being tested (and how many tests are included for the version of Waiwera you are running). If any tests fail, there will be output regarding which tests are not passing.

The tests will be repeated several times, using different numbers of processors each time (by default, using 1, 2, 3 and 4 processors). This can help show up any problems that occur only when running in parallel, or for specific numbers of processors.

.. index:: testing; benchmark tests

Benchmark tests
---------------

These are higher-level tests of the code as a whole, running complete flow simulation problems and comparing results against reference results (analytical solutions or output from other simulators). They may also be run to help verify correct installation of Waiwera.

Benchmark tests are carried out using a modified verison of the CREDO testing framework, originally developed for the `Underworld <http://www.underworldcode.org/>`_ geodynamics simulator.

.. index:: CREDO, CREDO; installation

Installing CREDO
^^^^^^^^^^^^^^^^
.. todo:: add instructions on installing CREDO

Running benchmark tests
^^^^^^^^^^^^^^^^^^^^^^^

To run all benchmark tests, navigate to the main Waiwera directory and type:

.. code-block:: bash

  python benchmark_tests.py

This will run the benchmark tests in serial. To run them in parallel with a specified number of processors, use the ``-np`` argument after the script name. For example, the benchmark tests can be run on two processors as follows:

.. code-block:: bash

  python benchmark_tests.py -np 2

As the tests are run, the path to each individual test script will be displayed, with a ``Pass`` or ``Fail`` after each one. When they are all finished, an overall ``Pass`` or ``Fail`` will be displayed.

In addition, an HTML page with a summary of the test results is written to the benchmark test directory: ``test/benchmark/test_summary.html``. This contains a list of all tests and a link to the results page for each one. These individual test results pages contain the results of all comparisons with reference solutions for each test, together with plots of relevant quantities.

It is also possible to run tests individually. The tests are all contained in sub-directories under the ``test/benchmark`` directory. Within each test directory there is a Python script named ``test_*.py`` (the exact name depending on the test). Executing this test script will run the test. These test scripts can also be run either in serial (the default) or in parallel by adding the ``-np`` argument after the script name, followed by the number of processors.
