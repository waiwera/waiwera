#!/usr/bin/env python

# configure Waiwera Meson build

import os
import argparse
import subprocess
from fruit_config import write_fruit_pkgconfig_file
from fson_config import write_fson_pkgconfig_file

env = os.environ.copy()
orig_path = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--debug", action = "store_true", help = "debug mode")
parser.add_argument("--release", action = "store_true", help = "release mode")
parser.add_argument("--fson_dir", help = "FSON library directory")
parser.add_argument("--fruit_dir", help = "FRUIT library directory")
parser.add_argument("--reconfigure", action = "store_true", help = "reconfigure")
parser.add_argument("--no_rpath", action = "store_true", help = "do not set RPATH in executable")
args = parser.parse_args()

if args.release: build_type = "release"
else: build_type = "debugoptimized" if args.debug else "release"

def env_update(key, value, separator = ' '):
    if key in env: env[key] += separator + value
    else: env[key] = value

fflags = " ".join(["-fPIC",
                   "-ffree-line-length-none",
                   "-Wno-unused-dummy-argument",
                   "-Wno-unused-function",
                   "-Wno-return-type",
                   "-Wno-maybe-uninitialized"])
env_update('FFLAGS', fflags)

# set pkg-config path for finding dependency libraries:
if "PETSC_DIR" in env and "PETSC_ARCH" in env:
    petsc_pkgconfig_path = os.path.join(env["PETSC_DIR"], env["PETSC_ARCH"],
                                        "lib", "pkgconfig")
    env_update('PKG_CONFIG_PATH', petsc_pkgconfig_path, ':')

if args.fson_dir:
    if os.path.isdir(args.fson_dir):
        fson_pkgconfig_path = os.path.join(args.fson_dir, "dist", "pkgconfig")
        if not os.path.isdir(fson_pkgconfig_path): os.mkdir(fson_pkgconfig_path)
        fson_pkgconfig_filename = os.path.join(fson_pkgconfig_path, "FSON.pc")
        if not os.path.isfile(fson_pkgconfig_filename):
            write_fson_pkgconfig_file(args.fson_dir)
        env_update('PKG_CONFIG_PATH', fson_pkgconfig_path, ':')
    else: raise Exception("Specified FSON library directory does not exist: " +
                          args.fson_dir)

if args.fruit_dir:
    if os.path.isdir(args.fruit_dir):
        fruit_pkgconfig_path = os.path.join(args.fruit_dir, "pkgconfig")
        if not os.path.isdir(fruit_pkgconfig_path): os.mkdir(fruit_pkgconfig_path)
        fruit_pkgconfig_filename = os.path.join(fruit_pkgconfig_path, "FRUIT.pc")
        if not os.path.isfile(fruit_pkgconfig_filename):
            write_fruit_pkgconfig_file(args.fruit_dir)
        env_update('PKG_CONFIG_PATH', fruit_pkgconfig_path, ':')
    else: raise Exception("Specified FRUIT library directory does not exist: " +
                          args.fruit_dir)

# unit_test_driver_source_filename = "test/unit/src/test_all.F90"
# if not os.path.isfile(unit_test_driver_source_filename):
#     open(unit_test_driver_source_filename, 'a').close()

install_prefix = os.path.expanduser("~")
install_dir = os.path.join(install_prefix, "bin")
if not os.path.isdir(install_dir): os.mkdir(install_dir)

set_rpath = 'false' if args.no_rpath else 'true'

os.chdir("build")
env["CC"] = "mpicc"; env["FC"] = "mpif90"

if args.reconfigure:
    subprocess.Popen(["ninja", "reconfigure"],
                     env = env).wait()
else:
    subprocess.Popen(["meson",
                      "--buildtype", build_type, "..",
                      "--prefix", install_prefix,
                      "-Dset_rpath=" + set_rpath],
                     env = env).wait()

os.chdir(orig_path)