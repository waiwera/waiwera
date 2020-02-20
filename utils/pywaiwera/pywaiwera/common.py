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
        return ''

__version__ = get_pkg_version()



