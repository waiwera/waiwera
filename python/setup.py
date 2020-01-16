import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="waiwera", # Replace with your own username
    version="0.0.3",
    author="Example Author",
    author_email="cyeh015@aucklanduni.ac.nz",
    description="Python package for the Waiwera simulator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://waiwera.github.io",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
    entry_points={'console_scripts': ['waiwera-dkr = waiwera.docker:main']},
)
