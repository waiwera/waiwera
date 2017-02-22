***************
Testing Waiwera
***************

Unit tests
----------

Unit tests are lower-level tests of individual subroutines in the code. These are used mostly during code development, but may also be run to help verify that Waiwera has been installed correctly.

Unit tests are carried out using the `FRUIT <https://sourceforge.net/projects/fortranxunit/>`_ framework for Fortran unit testing, via a separate Python interface called `FRUITPy <https://github.com/acroucher/FRUITPy>`_. You will need to install both FRUIT and FRUITPy before running the unit tests.

.. todo:: add instructions on installing FRUIT/FRUITPy
.. todo:: add instructions on running unit tests

Benchmark tests
---------------

These are higher-level tests of the code as a whole, running complete flow simulation problems and comparing results against reference results (analytical solutions or output from other simulators). They may also be run to help verify correct installation of Waiwera.

Benchmark tests are carried out using a modified verison of the CREDO testing framework, originally developed for the `Underworld <http://www.underworldcode.org/>`_ geodynamics simulator.

.. todo:: add instructions on installing CREDO
.. todo:: add instructions on running benchmark tests
