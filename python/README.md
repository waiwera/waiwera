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

# (Developer) compile and upload package

    python setup.py sdist bdist_wheel
    twine upload --repository-url https://test.pypi.org/legacy/ dist/*

    