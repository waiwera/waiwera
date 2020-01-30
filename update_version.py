# update version in Waiwera code source and config files
# usage: python update_version.py version
# where version is a string containing the version, e.g. '1.2.3'.

from __future__ import print_function
import sys
import os

def output(found, filename, lines):
    """Writes lines to specified file if found is True, otherwise prints
    warning"""
    if found:
        with open(filename, 'wt') as f:
            f.writelines(lines)
    else:
        print('Could not find version in', filename)

def meson(version):
    """Update version in meson build file"""
    filename = 'meson.build'
    lines = []
    found = False
    with open(filename, 'rt') as f:
        for line in f:
            if line.strip().startswith('version:'):
                parts = line.split("'")
                newline = "'".join([parts[0], version, parts[-1]])
                lines.append(newline)
                found = True
            else: lines.append(line)
    output(found, filename, lines)

def ford(version):
    """Update version in FORD config file"""
    filename = 'devdoc.md'
    lines = []
    found = False
    with open(filename, 'rt') as f:
        for line in f:
            if line.strip().startswith('version:'):
                newline = 'version: ' + version + '\n'
                lines.append(newline)
                found = True
            else: lines.append(line)
    output(found, filename, lines)

def src(version):
    """Update version in F90 source code"""
    filename = os.path.join('src', 'version.F90')
    lines = []
    found = False
    with open(filename, 'rt') as f:
        for line in f:
            if 'waiwera_version' in line:
                parts = line.split('"')
                newline = '"'.join([parts[0], version, parts[-1]])
                lines.append(newline)
                found = True
            else: lines.append(line)
    output(found, filename, lines)

def userdoc(version):
    """Update version in user docs"""
    filename = os.path.join('doc', 'user', 'conf.py')
    lines = []
    found_version, found_release = False, False
    with open(filename, 'rt') as f:
        for line in f:
            if line.strip().startswith('version'):
                short_version = '.'.join(version.split('.')[:-1])
                newline = "version = u'" + short_version + "'\n"
                lines.append(newline)
                found_version = True
            elif line.strip().startswith('release'):
                newline = "release = u'" + version + "'\n"
                lines.append(newline)
                found_release = True
            else: lines.append(line)
    output(found_version and found_release, filename, lines)

version = sys.argv[1]
for update in [meson, ford, src, userdoc]:
    update(version)
