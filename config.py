#!/usr/bin/env python

# configure Waiwera CMake build

from __future__ import print_function
import os
import argparse
import subprocess

env = os.environ.copy()
orig_path = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--debug", action = "store_true", help = "Debug mode")
parser.add_argument("--release", action = "store_true", help = "Release mode")
parser.add_argument("--fson_dir", help = "FSON library directory")
args = parser.parse_args()

if args.release: config_mode = "Release"
else: config_mode = "Debug" if args.debug else "Release"
print("Configure mode:", config_mode)

if args.fson_dir:
    if os.path.isdir(args.fson_dir):
        env["FSON_DIR"] = args.fson_dir
    else: raise Exception("Specified FSON library directory does not exist: " +
                          args.fson_dir)
unit_test_driver_source_filename = "test/unit/src/test_all.F90"
if not os.path.isfile(unit_test_driver_source_filename):
    open(unit_test_driver_source_filename, 'a').close()

os.chdir("build")

env["CC"] = "mpicc"; env["FC"] = "mpif90"
subprocess.Popen(["cmake", "..",
                  "-DCMAKE_BUILD_TYPE=" + config_mode],
                 env = env).wait()

os.chdir(orig_path)
