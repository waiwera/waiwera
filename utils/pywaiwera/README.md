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

PyWaiwera provides a console executable script `waiwera-dkr` for running Waiwera via [Docker](https://www.docker.com/).  It should be accesible directly as a console command.

To run a Waiwera model using Docker, please navigate to the input file's location and run the command:

    waiwera-dkr filename.json

where `filename.json` is the name of the model JSON input file.

Waiwera example problems to try can be obtained via the command:

    waiwera-dkr --examples

which will download a directory of example Waiwera benchmark test problems.

For more detailed help on the options available with `waiwera-dkr`, run:

    waiwera-dkr --help

# Waiwera simulator

[Official website](https://waiwera.github.io/)

[Source code repository](https://github.com/waiwera/waiwera)

[User guide](https://waiwera.readthedocs.io/en/latest/)

# Using the PyWaiwera package

    import pywaiwera
    print(pywaiwera.__version__)
    env = pywaiwera.docker.DockerEnv()
    env.run_waiwera(['input.json'])

Or you can import it differently

    from pywaiwera import docker


# (Developer) Compile and upload package

    python setup.py sdist bdist_wheel
    twine upload --repository-url https://test.pypi.org/legacy/ dist/*

In GitHub I am trying to use workflow/action, to trigger building and uploading of the package.  This is done by creating a workflow file `.github/workflows/pythonpublish.yml`:

    name: Build and Upload Python Package

    on:
      push:
        branches:
        - py-package

    jobs:
      deploy:
        runs-on: ubuntu-latest
        steps:
        - uses: actions/checkout@v1
        - name: Set up Python
          uses: actions/setup-python@v1
          with:
            python-version: '3.x'
        - name: Install dependencies
          run: |
            python -m pip install --upgrade pip
            pip install setuptools wheel twine
        - name: Build and publish
          env:
            TWINE_USERNAME: __token__
            TWINE_PASSWORD: ${{ secrets.PYPI_PASSWORD }}
          run: |
            cd python
            python setup.py sdist bdist_wheel
            twine upload --repository-url https://test.pypi.org/legacy/ dist/*

See [PyPA page](https://packaging.python.org/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/) and [GitHub example code](https://github.com/marketplace/actions/pypi-publish).

Also needs to setup [secrets](https://help.github.com/en/actions/automating-your-workflow-with-github-actions/creating-and-using-encrypted-secrets#creating-encrypted-secrets).  Just add the value of `PYPI_PASSWORD`.

It is possible to check package version by `pip freeze`.  However it is slightly more difficult to have a builtin version number that can be displayed when the console script or module is used.  The following code would work, but requires user to have setuptools, which is often available, but not quranteed.

    import pkg_resources
    dist = pkg_resources.get_distribution("waiwera")
    VERSION = dist.version

# (Developer) Install PyWaiwera as developer

It is possible to install the PyWaiwera package while developing.

  pip install -e utils/pywaiwera

This will install the package in editable mode (develop mode) from a local project path or a version control url.
