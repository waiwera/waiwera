import setuptools
import re

with open("README.md", "r") as fh:
    long_description = fh.read()

# use official waiwera version number from fortran source
with open('../../src/version.F90', 'r') as fv:
    matches = re.findall('waiwera_version += +"(.+?)"', fv.read())
    if matches:
        wai_version = matches[0]
    else:
        raise Exception('Unable to find waiwera version string.')

setuptools.setup(
    name="pywaiwera",
    version=wai_version,
    author="Waiwera Project",
    author_email="waiwera.project@gmail.com",
    description="Python package for the Waiwera simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://waiwera.github.io",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
    entry_points={'console_scripts': ['waiwera-dkr = pywaiwera.docker:main']},
)
