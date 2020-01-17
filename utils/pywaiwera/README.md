Waiwera Python Package
======================

This is the official python package of the Waiwera simulator.

# Install

The easiest way to install is using pip:

    pip install --index-url https://test.pypi.org/simple/ waiwera

Console executable scripts are supplied with this package, please ensure the install location of python scripts is added to PATH.

# Running Waiwera using Docker

A console executable script `waiwera-dkr` is supplied with this package.  It should be accesible directly as a console command.

Waiwera examples can be obtained via the command:

    waiwera-dkr -e

To run a model, please navigate to the input file's location, and run the command:

    waiwera-dkr input.json

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


# (Developer) compile and upload package

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

