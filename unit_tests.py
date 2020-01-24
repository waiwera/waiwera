# Build and run Waiwera tests. Optional command line arguments are the
# test module names required in the test, otherwise all test modules
# are used. Named arguments specify the command for executing tests
# (e.g. mpiexec), number of processes, and name of build directory to
# run tests in.

from __future__ import print_function
import sys
import os
import multiprocessing
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("tests", default = [], help = "test modules", nargs="*")
parser.add_argument("--exe", default = "mpiexec", help = "command for executing unit test")
parser.add_argument("--procs", default = "np", help = "option for specifying number of processes")
parser.add_argument("--builddir", default = "build", help = "build subdirectory")
args = parser.parse_args()

test_modules = [m + '_test' for m in args.tests]

os.chdir(args.builddir)

try:
    num_available_procs = multiprocessing.cpu_count()
except (NotImplementedError):
    num_available_procs = 1
num_procs = range(1, min(num_available_procs, 4) + 1)

ret = 0
for np in num_procs:
    print('-'*72)
    print('Processors:', np, '/', num_procs[-1])
    cmd = ['meson', 'test'] + test_modules + \
          ['--wrap="' + ' '.join([args.exe, "-" + args.procs, str(np)]) + '"']
    ret += subprocess.call(" ".join(cmd), shell = True)

sys.exit(ret)
