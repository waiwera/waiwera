# this requires setuptools
import pkg_resources
import os
import re

def get_pkg_version():
    """ use python's package info, if loaded as a packaged module """
    try:
        dist = pkg_resources.get_distribution("pywaiwera")
        return dist.version
    except pkg_resources.DistributionNotFound:
        return None

def get_f90_version():
    """ use official waiwera version number from fortran source, if loaded from
    PYTHONPATH """
    here = os.path.dirname(os.path.abspath(__file__))
    vf90 = os.path.abspath(os.path.join(here, '..', '..', '..', 'src', 'version.F90'))
    if not os.path.isfile(vf90):
        raise Exception('Unable to find waiwera/src/version.F90.')
    with open('../../src/version.F90', 'r') as fv:
        matches = re.findall('waiwera_version += +"(.+?)"', fv.read())
        if matches:
            return matches[0]
        else:
            raise Exception('Unable to parse version string from version.F90.')

__version__ = get_pkg_version()
if __version__ is None:
    __version__ = get_f90_version()



