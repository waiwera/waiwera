# Build and run Waiwera tests. Optional command line arguments are the test module
# names required in the test, otherwise all test modules are used.

from __future__ import print_function
from sys import argv
import os
import multiprocessing
import subprocess

test_modules = [m + '_test' for m in argv[1:]]

os.chdir('build')

try:
    num_available_procs = multiprocessing.cpu_count()
except (NotImplementedError):
    num_available_procs = 1
num_procs = range(1, min(num_available_procs, 4) + 1)

for np in num_procs:
    print('-'*72)
    print('Processors:', np, '/', num_procs[-1])
    cmd = ['meson', 'test'] + test_modules + \
          ['--wrap="mpiexec -np ' + str(np) + '"']
    subprocess.call(" ".join(cmd), shell = True)
