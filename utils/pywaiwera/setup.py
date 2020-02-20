import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pywaiwera",
    version="1.1.0",
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
