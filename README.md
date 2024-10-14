Waiwera
=======

[![Unit and benchmark tests](https://github.com/waiwera/waiwera/actions/workflows/test.yml/badge.svg?branch=testing)](https://github.com/waiwera/waiwera/actions/workflows/test.yml) ![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/waiwera/waiwera) ![Docker Pulls](https://img.shields.io/docker/pulls/waiwera/waiwera?color=green) [![Documentation Status](https://readthedocs.org/projects/waiwera/badge/?version=latest)](https://waiwera.readthedocs.io/en/latest/?badge=latest)

Waiwera is a parallel, open-source geothermal flow simulator.

Waiwera features:

- numerical simulation of high-temperature subsurface flows, including robust phase changes
- parallel execution on shared- or distributed-memory computers and clusters
- use of [PETSc](https://petsc.org/) (Portable Extensible Toolkit for Scientific Computation) for parallel data structures, linear and non-linear solvers, etc.
- standard file formats for input ([JSON](http://www.json.org)) and output ([HDF5](https://portal.hdfgroup.org/display/HDF5/HDF5), [YAML](http://www.yaml.org/about.html))
- structured, object-oriented code written in Fortran 2003
- free, open-source license (GNU LGPL)

Waiwera was developed at the University of Auckland's [Geothermal Institute](http://www.geothermal.auckland.ac.nz/). Initial development was part of the "Geothermal Supermodels" research project, funded by the NZ Ministry of Business, Innovation and Employment ([MBIE](https://www.mbie.govt.nz/)), with additional support from [Contact Energy Ltd](https://contact.co.nz/).

The word *Waiwera* comes from the Māori language and means "hot water".

Further information can be found on the Waiwera [website](https://waiwera.github.io/) and in the [user guide](https://waiwera.readthedocs.io/).
