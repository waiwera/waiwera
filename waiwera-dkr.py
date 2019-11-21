#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import subprocess
import posixpath
import shlex
import sys

REPO = 'waiwera/waiwera'
TAG = 'latest'
WAIWERA_PATH = '/opt/waiwera'
VERSION = '0.1'

def convert_path_nt_to_posix(path):
    return '/{0}'.format(path.replace("\\","/").replace(":",""))

##########
# A simple python wrapper around the docker image for waiwera
#########

# docker --check


# list check for write
# bind succesfulu
# list

# function instead of main
# --version
# User guide

# SEGFAULTS user guide

# PETSC check for version
# Update
# file sharing directory

if __name__ == "__main__":
    """
    Args:
    """
    parser = argparse.ArgumentParser(description='Runs Waiwera, \
                        the open-source geothermal flow simulator')
    parser.add_argument('command', nargs='+', help='the command passed to \
                        waiwera')
    parser.add_argument('-np', '--num_processors', help='The number of \
                        processors to utilize, otherwise uses the docker \
                        default for your system')
    parser.add_argument('-r', '--repo',
                        default=REPO)
    parser.add_argument('-t', '--tag',
                        default=TAG)
    parser.add_argument('-i', '--image', help='The docker image to use \
                        e.g. waiwera/waiwera:latest')
    parser.add_argument('-it','--interactive',
                        help='starts an interactive terminal and does NOT run \
                        mpiexec by default',
                        action='store_true')
    parser.add_argument('-u','--noupdate',
                    help='stops the script pulling an image update',
                    action='store_true')


    args = parser.parse_args()

    current_path = os.getcwd()
    if sys.platform == 'win32':
        current_path = convert_path_nt_to_posix(current_path)

    data_path = '/data'

    if args.image == None:
        image = '{0}:{1}'.format(args.repo, args.tag)
    else:
        image = args.image

    if args.num_processors:
        np ='-np {0}'.format(args.num_processors)
    else:
        np = ''

    if args.interactive:
        it ='--interactive --tty'
        mpiexec = ''
        work_dir = ''
    else:
        it  = ''
        # Change the working directory to
        build_path = posixpath.join(WAIWERA_PATH, 'build', 'waiwera')
        mpiexec = 'mpiexec {0} {1}'.format(np, build_path)
        work_dir = '--workdir {0}'.format(data_path)

    command = ''
    for string in args.command:
        command = '{0} {1}'.format(command, string)

    if not args.noupdate:
        print('Checking for Waiwera update')
        pull_cmd = "docker pull {0}".format(image)
        os.system(pull_cmd)
        # print(subprocess.check_output(shlex.split(pull_cmd)))

    fo = open("uidcheck.txt", "wb")
    fo.close()

    #  docker run -v ${p}:/data -w /data waiwera-phusion-debian mpiexec -np $args[1] /home/mpirun/waiwera/dist/waiwera $args[0]
    print('Running Waiwera')
    run_cmd = "docker run --rm {0} --volume {1}:{2} {3} {4} \
                {5} {6}" \
                .format(it,
                    current_path,
                    data_path,
                    work_dir,
                    image,
                    mpiexec,
                    command)
    print(run_cmd)
    # print(run_cmd)
    os.system(run_cmd)
    #print(subprocess.check_output(shlex.split(run_cmd)))
