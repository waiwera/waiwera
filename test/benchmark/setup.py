# set up reference solutions for all benchmark tests
from __future__ import print_function

import os
import glob

setupfilename = "setup.py"

orig_path = os.getcwd()
for path, dirs, files in os.walk(os.path.curdir):
    os.chdir(path)
    scripts = glob.glob("test_*.py")
    if scripts and setupfilename in files:
        print("Setting up", path)
        execfile(setupfilename)
    os.chdir(orig_path)
print("done.")
