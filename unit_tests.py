# Build and run Waiwera tests. Optional command line arguments are the test module
# names required in the test, otherwise all test modules are used.

from __future__ import print_function
from sys import argv
import os
import multiprocessing
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("tests", default = [], help = "test modules", nargs="*")
parser.add_argument("--exe", default = "mpiexec", help = "command for executing unit test")
parser.add_argument("--procs", default = "np", help = "option for specifying number of processes ")
args = parser.parse_args()

test_modules = [m + '_test' for m in args.tests]

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
          ['--wrap="' + ' '.join([args.exe, "-" + args.procs, str(np)]) + '"']
    subprocess.call(" ".join(cmd), shell = True)
