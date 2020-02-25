Waiwera Python Package
======================

PyWaiwera is the official Python package for use with the parallel open-source geothermal flow simulator, [Waiwera](https://waiwera.github.io/).

# Installation

The easiest way to install PyWaiwera is using `pip`:

    pip install pywaiwera

or if you do not have permissions to install system-wide Python packages, you can install it locally in your own user account:

    pip install --user pywaiwera

During the installation it may warn you that executable scripts are being installed to a directory that is not listed in your system's ``PATH`` environment variable. This means that you need to add this directory to your ``PATH`` if you want to be able to run the supplied scripts from anywhere on your machine.

PyWaiwera can be upgraded to the latest version at any time by running:

    pip install --upgrade pywaiwera

# Running Waiwera using Docker

PyWaiwera provides a console executable script `waiwera-dkr` for running Waiwera via [Docker](https://www.docker.com/).  It should be accessible directly as a console command.

To run a Waiwera model using Docker, please navigate to the input file's location and run the command:

    waiwera-dkr filename.json

where `filename.json` is the name of the model JSON input file.

Waiwera example problems to try can be obtained via the command:

    waiwera-dkr --examples

which will download a directory of example Waiwera benchmark test problems.

For more detailed help on the options available with `waiwera-dkr`, run:

    waiwera-dkr --help

# Using PyWaiwera in a script

It is also possible to use the PyWaiwera package from within Python scripts. For example, the following script imports PyWaiwera, creates a Docker environment and uses it to run Waiwera on the model input file `input.json`:

    import pywaiwera
    print(pywaiwera.__version__)
    env = pywaiwera.docker.DockerEnv()
    env.run_waiwera('input.json')
