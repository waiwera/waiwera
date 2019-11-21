.. index:: Waiwera; installation

******************
Installing Waiwera
******************

Waiwera can be run on most operating systems (e.g. Linux, Windows, Mac OS) using `Docker <https://www.docker.com/>`_, and for many users this is likely to be the easiest option. Linux users also have the option of building a native Waiwera executable.

Using Docker
============

Docker is a "`containerisation <https://www.docker.com/resources/what-container>`_" technology, in which a software application is deployed in a standard unit containing not only the code itself but all its dependencies as well. Because a container isolates the software from its environment, it should always run the same way, regardless of which kind of operating system it is run on.

Installing Docker
-----------------

To run Waiwera this way, first Docker must be installed on your system. You will need administrator privileges to install Docker. The Docker website has instructions for installing it on `Windows <https://docs.docker.com/docker-for-windows/>`_ and `Mac OS <https://docs.docker.com/docker-for-mac/>`_. On Linux systems, Docker may be installed directly from your pacakge manager. Alternatively (e.g. if you need a newer version than your package manager is able to provide) you can install it from the Docker repositories. The Docker website has instructions on how to do that for various Linux distributions (e.g. `Ubuntu <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_, `Debian <https://docs.docker.com/install/linux/docker-ce/debian/>`_).

Once Docker is installed, by default it should have access to all the processors on your system, a sensible amount of RAM, and the directories you are likely to want to run it in. However, it is possible to configure these settings if needed. The Docker Desktop application on Windows and Mac OS provides a graphical interface for doing this. You may also need to add your user account to the 'docker' user group in order to be able to run Docker containers without administrator privileges. Check the `Docker documentation <https://docs.docker.com/>`_ for more details on configuring Docker.

For Windows users, Windows 10 is recommended for running Waiwera via Docker. On Windows 7, it should still work, but Docker uses different underlying technology (based on virtual machines) which is slower and less reliable.

Running Waiwera using Docker
----------------------------

Waiwera provides a Python script to simplify running the Waiwera Docker container. This script will check if the Waiwera Docker container image is already installed on your system, and if not, it will automatically install it before running your Waiwera model. It also handles the sharing of files between the Docker container and your system. For more details, see *TODO link to running Docker section*.

Native Linux build
==================

For building a native Waiwera executable on Linux, Waiwera uses the `Ansible <https://www.ansible.com/>`_ deployment system, which automates the build process. This includes checking if the necessary tools (e.g. compilers, build tools) are present on your system, installing them if they are not, building Waiwera's dependency libraries (e.g. `PETSc <https://www.mcs.anl.gov/petsc/>`_), and building Waiwera itself.

Install Ansible
---------------

First, Ansible itself must be installed. Ansible is Python-based, so it can be installed either via your system package manager (e.g. ``sudo apt install ansible`` on Debian-based systems), or via `PYPI <https://pypi.org/>`_ and pip. For more details, consult the Ansible `documentation <https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html>`_.

Download the Waiwera source code
--------------------------------

The Waiwera source code is version-controlled using `Git <https://git-scm.com/>`_ and hosted on GitHub. The easiest way to download the code is by cloning the Waiwera `Git repository <https://github.com/waiwera/waiwera>`_. You will need Git installed on your system first. To clone the Waiwera repository, make a new directory for the Waiwera code, and in it execute the following command:

.. code-block:: bash

   git clone git@github.com:waiwera/waiwera.git .

Alternatively, you can download a ZIP archive of the code `here <https://github.com/waiwera/waiwera/archive/master.zip>`_.

Build Waiwera
-------------

Finally, build Waiwera by executing:

*TODO. Need to go to install/ directory?*

.. code-block:: bash

   ansible-playbook etc etc

*TODO: does this also install as well as build?*

Running the unit tests
----------------------

You can check the Waiwera build by running the unit tests. In the Waiwera base directory, execute:

.. code-block:: bash

   python unit_tests.py

This will run the Waiwera unit tests on 1, 2, 3 and 4 processes (or up to the number of processes available, if that is less than 4).

*TODO: how to tell if the unit tests have passed*
