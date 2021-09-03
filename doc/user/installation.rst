.. index:: Waiwera; installation, installation

******************
Installing Waiwera
******************

Waiwera can be run on most operating systems (e.g. Linux, Windows, macOS) using `Docker <https://www.docker.com/>`_, and for many users this is likely to be the easiest option. Linux users also have the option of building a native Waiwera executable (see :ref:`native_linux_build`).

.. index:: Docker
.. _using_docker:

Using Docker
============

Docker is a "`containerisation <https://www.docker.com/resources/what-container>`_" technology, in which a software application is deployed in a standard unit containing not only the code itself but all its dependencies as well. Because a container isolates the software from its environment, it should always run the same way, regardless of which kind of operating system it is run on.

.. index:: Docker; installation

Installing Docker
-----------------

To run Waiwera this way, first Docker must be installed on your system. You will need administrator privileges to install Docker.

Linux
.....

On Linux systems, Docker may be installed directly from your package manager. Alternatively (e.g. if you need a newer version than your package manager is able to provide) you can install it from the Docker repositories. The Docker website has instructions on how to do that for various Linux distributions (e.g. `Ubuntu <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_, `Debian <https://docs.docker.com/install/linux/docker-ce/debian/>`_). It is also important to follow the `post-install instructions <https://docs.docker.com/install/linux/linux-postinstall/>`_, particularly to add your user account to the ``docker`` group (so you can run Docker containers without root privileges).

Windows
.......

For running Waiwera on Windows using Docker, the Windows 10 Pro, Enterprise or Education versions are recommended. Users of these versions can install `Docker Desktop for Windows <https://docs.docker.com/docker-for-windows/install/>`_, which provides the best performance and also an easy-to-use graphical interface for configuring Docker.

Users of the Windows 7, 8 or 10 Home versions cannot use Docker Desktop, but can use the older `Docker Toolbox for Windows <https://docs.docker.com/toolbox/toolbox_install_windows/>`_ instead. This has a higher performance overhead than Docker Desktop, so Waiwera will run more slowly. Installing and configuring it are also not as convenient.

macOS
.....

For running Waiwera on macOS using Docker, macOS version 10.13 or newer is recommended, in which case `Docker Desktop for macOS <https://docs.docker.com/docker-for-mac/install/>`_ can be installed. This provides the best performance and also an easy-to-use graphical interface for configuring Docker.

Users of older versions of macOS cannot install Docker Desktop, but can use the older `Docker Toolbox for macOS <https://docs.docker.com/toolbox/toolbox_install_mac/>`_ instead. This has a higher performance overhead than Docker Desktop, so Waiwera will run more slowly. Installing and configuring it are also not as convenient.

.. index:: Docker; configuring

Configuring Docker
..................

Once Docker is installed, by default it should have access to all the processors on your system, a sensible amount of RAM, and the directories you are likely to want to run it in. However, it is possible to configure these settings if needed, for example if you have other drives you want to make available to Docker.

The Docker Desktop application on Windows and Mac OS provides a graphical interface for doing this.

.. index:: Docker; running, PyWaiwera; installation

Installing the PyWaiwera Python library
---------------------------------------

The simplest way to run Waiwera via Docker is to use the ``waiwera-dkr`` executable script which is part of `PyWaiwera <https://pypi.org/project/pywaiwera>`_, a Python library for simplifying tasks related to working with Waiwera. To use this, you will need to:

- have the `Python <https://www.python.org/>`_ scripting language installed on your machine
- install the PyWaiwera library

The easiest way to install PyWaiwera is using ``pip``, Python's package manager:

.. code-block:: bash

   pip install pywaiwera

or if you don't have permissions for installing system-wide Python packages, you can just install it locally inside your own user account:

.. code-block:: bash

   pip install --user pywaiwera

This will download and install PyWaiwera from the Python Package Index (`PyPI <https://pypi.org>`_). During the installation it may warn you that executable scripts are being installed to a directory that is not listed in your system's ``PATH`` environment variable. This means that you need to add this directory to your ``PATH`` if you want to be able to run it from anywhere on your machine.

You can also upgrade PyWaiwera to the latest version at any time by running:

.. code-block:: bash

   pip install --upgrade pywaiwera

.. index:: Docker; running

Running the Waiwera Docker container
------------------------------------

The `waiwera-dkr` script will check if the Waiwera Docker container image is already installed on your system, and if not, it will automatically install it before running your Waiwera model. It also handles the sharing of files between the Docker container and your system. For more details, see :ref:`run_docker`.

.. index:: Waiwera; building
.. _native_linux_build:

Native Linux build
==================

For building a native Waiwera executable on Linux, Waiwera uses the `Ansible <https://www.ansible.com/>`_ deployment system, which automates the build process. This includes checking if the necessary tools (e.g. compilers, build tools) are present on your system, installing them if they are not, building Waiwera's dependency libraries (e.g. `PETSc <https://www.mcs.anl.gov/petsc/>`_), and building Waiwera itself (which is carried out using the `Meson <https://mesonbuild.com/>`_ build system).

.. index:: Ansible

Install Ansible
---------------

First, Ansible itself must be installed. Ansible is Python-based, so it can be installed either via your system package manager (e.g. ``sudo apt install ansible`` on Debian-based systems), or via `PyPI <https://pypi.org/>`_ and `pip`. For more details, consult the Ansible `documentation <https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html>`_.

Note that the Waiwera build requires Ansible version 2.4 or later.

Download the Waiwera source code
--------------------------------

The Waiwera source code is version-controlled using `Git <https://git-scm.com/>`_ and hosted on GitHub. The easiest way to download the code is by cloning the Waiwera `Git repository <https://github.com/waiwera/waiwera>`_. You will need Git installed on your system first. To clone the Waiwera repository, make a new directory for the Waiwera code, and in it execute the following command:

.. code-block:: bash

   git clone https://github.com/waiwera/waiwera.git .

Alternatively, you can download a ZIP archive of the code `here <https://github.com/waiwera/waiwera/archive/master.zip>`_.

Build Waiwera
-------------
In the downloaded Waiwera source code repository, navigate to the install directory:

.. code-block:: bash

   cd install

To build and install Waiwera (and its various dependencies, if needed), execute the command:

.. code-block:: bash

   ansible-playbook ansible/install_local.yml

Root privileges are often required if, for example, additional software packages (e.g. compilers) need to be installed as part of the build process. This can be done by using the ``--ask-become-pass`` option:

.. code-block:: bash

   ansible-playbook ansible/install_local.yml --ask-become-pass

in which case you will be prompted for the ``sudo`` password at the start of the build.

Waiwera will be installed into the user's home directory (``~/bin/waiwera``) by default. Additional command line parameters can be used to customise the installation. For example:

.. code-block:: bash

   ansible-playbook ansible/install_local.yml -e "base_dir=/home/fred/waiwera/"

where ``base_dir`` is the desired Waiwera installation directory (in this example, ``/home/fred/waiwera/``).

Other example variables include :

* ``petsc_update=true`` will build a new version of PETSc even if an installed version is detected
    * defaults to ``false`` meaning PETSc will only be built if an installed version isn't detected
* ``waiwera_update=true`` will build Waiwera every time even a new version isn't pulled by git
    * defaults to ``false``
* ``zofu_build=true``
    * defaults to ``false`` and uses meson to build zofu
* ``fson_build=true``
    * defaults to ``false`` and uses meson to build zofu
* ``ninja_build=true``
    * defaults to ``false`` and only builds locally if no ninja install is detected

The following command builds Waiwera (and dependencies), but does not install it.

.. code-block:: bash

  ansible-playbook ansible/build.yml

.. index:: PETSc; linking

Linking to your own PETSc library
---------------------------------

The build system will download and build its own version of the PETSc library, if it does not detect one already present on your computer. Note that Waiwera currently requires PETSc version 3.15.2 or newer.

If you do have a suitably up-to-date version of PETSc on your machine, and want to link to it rather than having Waiwera build its own, you need to make sure your PETSc is discoverable by the ``pkg-config`` tool, before you start the build.

This means that you need a ``pkg-config`` file called ``PETSc.pc`` somewhere on your system. You can generally find such a file in the ``lib/pkgconfig`` subdirectory of your PETSc build directory. You also need to have the directory containing this file included in your ``PKG_CONFIG_PATH`` environment variable.

You can check the location of the PETSc library that ``pkg-config`` is going to detect using the command:

.. code-block:: bash

   pkg-config --libs PETSc

It is also possible to check the PETSc version using:

.. code-block:: bash

   pkg-config --modversion PETSc

.. index:: FSON; linking, Zofu; linking

Linking to other libraries
--------------------------

Similar considerations apply to the other two libraries that Waiwera links to (`FSON <https://github.com/josephalevin/fson>`_ and `Zofu <https://github.com/acroucher/zofu>`_). If you have pre-existing versions already built on your system that you want to link to, just make sure there are ``pkg-config`` files for them in locations on your ``PKG_CONFIG_PATH``. Otherwise, the Waiwera build system will build its own versions.

.. index:: Waiwera; setting path

Setting the path
----------------
To ensure Waiwera can be run from any directory, the Waiwera installation directory should be on the user's ``PATH``. If this is not the case then it can be added by editing the ``PATH`` environment variable. To make this change persistent, add the path modification to your login script (e.g. ``~/.profile``).

.. index:: testing; unit tests, Zofu

Running the unit tests
----------------------

You can check the Waiwera build by running the unit tests. The unit tests (which test individual routines in the Waiwera code) are created using the `Zofu <https://github.com/acroucher/zofu>`_ framework for Fortran unit testing, and run using Meson. In the Waiwera base directory, execute:

.. code-block:: bash

   python unit_tests.py

This will run the Waiwera unit tests on 1, 2, 3 and 4 processes (or up to the number of processes available, if that is less than 4).

It is also possible to run subsets of the unit tests by specifying the module names, e.g.:

.. code-block:: bash

  python unit_tests.py IAPWS

which tests only the `IAPWS` module, or:

.. code-block:: bash

  python unit_tests.py face cell

which tests only the `face` and `cell` modules.

If the tests have successfully passed, the unit test output will appear something like this:

.. code-block:: bash

  Ok:                   32
  Expected Fail:         0
  Fail:                  0
  Unexpected Pass:       0
  Skipped:               0
  Timeout:               0

The precise numbers of asserts and cases will vary, depending on how many modules are being tested (and how many tests are included for the version of Waiwera you are running). If any tests fail, there will be output regarding which tests are not passing.

Installing Waiwera on your system
---------------------------------

From the Waiwera root directory, the Waiwera executable can be
installed on your system as follows:

.. code-block:: bash

   ninja -C build install

It can subsequently be uninstalled using:

.. code-block:: bash

   ninja -C build uninstall

..
   Section on cluster install?

..
   --mpi_wrapper_compiler option in config?

..
   By default, parallel unit test runs will be carried out using the `mpiexec` command, with the number of processes specified using the `-np` option. These can be changed by passing the `exe` and `procs` parameters to the `unit_tests.py` script. For example, if you are running the tests on a compute cluster and need to submit them via the `Slurm <https://slurm.schedmd.com/>`_ workload manager, the unit tests might be run using a command like this:

   .. code-block:: bash

     python unit_tests.py mesh --exe "srun --qos=debug -A acc00100 --time=2:00 --mem-per-cpu=100" --procs "n"
