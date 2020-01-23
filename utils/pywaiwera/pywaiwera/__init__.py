from .common import __version__

# allows use of `from pywaiwera import *`
__all__ = ['common', 'docker']
# allows use tof `import pywaiwera; pywaiwera.docker`
from . import *
