from importlib import metadata
import os
import re

def get_pkg_version():
    """ use python's package info, if loaded as a packaged module """
    try:
        return metadata.version("pywaiwera")
    except metadata.PackageNotFoundError:
        return ''

__version__ = get_pkg_version()



