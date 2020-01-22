#!/usr/bin/env python

# configure Waiwera Meson build, and configure/build PETSc library if needed

from __future__ import print_function
import os
import argparse
import subprocess
from shutil import copyfile

env = os.environ.copy()
base_dir = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("--builddir", default = "build", help = "build subdirectory")
parser.add_argument("-d", "--debug", action = "store_true", help = "debug mode")
parser.add_argument("--release", action = "store_true", help = "release mode")
parser.add_argument("--no_rpath", action = "store_true", help = "do not set RPATH in executable")
parser.add_argument("--prefix", default = os.path.expanduser("~"), help = "prefix for installation path")
parser.add_argument("--libdir", default = "lib", help = "library installation directory")
parser.add_argument("--petsc_revision", default = "339143f02b", help = "PETSc git revision")
parser.add_argument("--mpi_wrapper_compiler", default = False, help = "Use MPI wrapper compiler (specified in FC environment variable)")
args = parser.parse_args()

if args.release: build_type = "release"
else: build_type = "debug" if args.debug else "release"

def env_update(key, value, separator = ' '):
    """Update environment variable if necessary"""
    if key in env:
        if value not in env[key]: env[key] += separator + value
    else: env[key] = value

def petsc_found():
    """Returns true if PETSc library is discoverable by pkg-config and
    PETSc environment variables are set"""
    return subprocess.Popen(["pkg-config", "--exists", "PETSc"]).wait() == 0

def petsc_mkdir():
    """Creates directory for building PETSc and returns the path to it"""
    petsc_dir = os.path.join(base_dir, "external", "PETSc")
    if not os.path.isdir(petsc_dir): os.makedirs(petsc_dir)
    return petsc_dir

def petsc_download(petsc_dir):
    """Downloads PETSc to external projects directory"""
    petsc_url = "https://gitlab.com/petsc/petsc.git"
    print("PETSc library not found: installing PETSc in", petsc_dir)
    os.chdir(petsc_dir)
    subprocess.Popen(["git", "clone",
                      petsc_url, "."]).wait()
    subprocess.Popen(["git", "checkout", args.petsc_revision]).wait()
    os.chdir(base_dir)

def petsc_config(petsc_dir):
    """Configure PETSc"""
    env["PETSC_DIR"] = petsc_dir
    env["PETSC_ARCH"] = build_type
    os.chdir(petsc_dir)
    download_pkgs = ["hdf5", "pnetcdf", "zlib", "netcdf",
                         "exodusii", "triangle", "ptscotch", "chaco"]
    if build_type == "release":
        flags = ["--with-debugging=0",
                 "COPTFLAGS='-O3 -march=native -mtune=native'",
                 "CXXOPTFLAGS='-O3 -march=native -mtune=native'",
                 "FOPTFLAGS='-O3 -march=native -mtune=native'"]
    else: flags = []
    subprocess.Popen(["./configure"] +
                     ["--download-" + pkg for pkg in download_pkgs] +
                     flags, env = env).wait()
    os.chdir(base_dir)

def petsc_build(petsc_dir):
    """Build PETSc"""
    print("Building PETSc:")
    os.chdir(petsc_dir)
    subprocess.Popen(["make", "all"], env = env).wait()
    os.chdir(base_dir)

def petsc_check(petsc_dir):
    """Run basic PETSc check"""
    os.chdir(petsc_dir)
    subprocess.Popen(["make", "check"], env = env).wait()
    os.chdir(base_dir)

def copy_petsc_pkgconfig():
    """Copy PETSc pkgconfig file into pkgconfig directory"""
    petsc_pkgconfig_file = os.path.join(env["PETSC_DIR"], env["PETSC_ARCH"],
                                        "lib", "pkgconfig", "PETSc.pc")
    dest_petsc_pkgconfig_file = os.path.join(pkg_config_dir, "PETSc.pc")
    if os.path.isfile(petsc_pkgconfig_file):
        if not os.path.exists(pkg_config_dir):
            os.makedirs(pkg_config_dir)
        print("Copying PETSc.pc to", pkg_config_dir)
        copyfile(petsc_pkgconfig_file, dest_petsc_pkgconfig_file)

# set pkg-config path:
pkg_config_dir = os.path.join(args.prefix, args.libdir, 'pkgconfig')
env_update('PKG_CONFIG_PATH', pkg_config_dir, ':')

if not petsc_found():
    petsc_dir = petsc_mkdir()
    petsc_download(petsc_dir)
    petsc_config(petsc_dir)
    petsc_build(petsc_dir)
    petsc_check(petsc_dir)
    copy_petsc_pkgconfig()

set_rpath = 'false' if args.no_rpath else 'true'

build_dir = args.builddir
if not os.path.isdir(build_dir): os.mkdir(build_dir)
os.chdir(build_dir)

meson_args = [
    "--buildtype", build_type, "..",
    "--prefix", args.prefix,
    "-Dlibdir=" + args.libdir,
    "-Dset_rpath=" + set_rpath]

if args.mpi_wrapper_compiler:
    meson_args.append("-Dmpi_wrapper_compiler=true")

subprocess.Popen(["meson"] + meson_args, env = env).wait()

os.chdir(base_dir)
