# this requires setuptools
import pkg_resources

def get_pkg_version():
    dist = pkg_resources.get_distribution("waiwera")
    return dist.version

__version__ = get_pkg_version()

