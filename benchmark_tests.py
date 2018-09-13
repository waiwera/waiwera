# run all benchmark tests and produce a summary report
from __future__ import print_function

tests_path = "test/benchmark"

import os
import glob

orig_path = os.getcwd()
for path, dirs, files in os.walk(tests_path):
    os.chdir(path)
    scripts = glob.glob("test_*.py")
    for script in scripts:
        print('running test:', script)
        execfile(script)
    os.chdir(orig_path)
