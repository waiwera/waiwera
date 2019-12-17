#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import subprocess
import posixpath
import shlex
import sys
import signal

REPO = 'waiwera/waiwera'
TAG = 'latest'
WAIWERA_PATH = '/opt/waiwera'
VERSION = '0.2'
CID_LEN = 12

##########
# A simple python wrapper around the docker image for waiwera
#########

def signal_handler(sig, frame):
        if os.path.isfile('.cid'):
            with open('.cid', 'r') as f:
                cid = f.readline().strip()[:CID_LEN]
            print('You pressed Ctrl+C! Killing Waiwera container %s' % cid)
            # os.system('docker ps')
            print('docker kill %s' % cid)
            os.system('docker kill %s' % cid)
            # os.system('docker ps')
            os.remove('.cid')
        sys.exit(1)

def convert_path_nt_to_posix(path):
    return '/{0}'.format(path.replace("\\","/").replace(":",""))

def waiwera_docker(args):
    current_path = os.getcwd()
    if sys.platform == 'win32':
        current_path = convert_path_nt_to_posix(current_path)

    data_path = '/data'

    if args.image == None:
        image = ['{0}:{1}'.format(args.repo, args.tag)]
    else:
        image = [args.image]

    if args.num_processors:
        np = ['-np', '{}'.format(args.num_processors)]
    else:
        np = ['']

    if args.interactive:
        it = ['--interactive', '--tty']
        work_dir = ['']
        mpiexec = ['']
    else:
        it  = ['']
        work_dir = ['--workdir', data_path]
        mpiexec = ['mpiexec'] + np + [posixpath.join(WAIWERA_PATH, 'build', 'waiwera')]

    if not args.noupdate:
        print('Checking for Waiwera update')
        pull_cmd = "docker pull {0}".format(image[0])
        os.system(pull_cmd)
        # print(subprocess.check_output(shlex.split(pull_cmd)))

    fo = open(".idcheck", "wb")
    fo.close()

    #  docker run -v ${p}:/data -w /data waiwera-phusion-debian mpiexec -np $args[1] /home/mpirun/waiwera/dist/waiwera $args[0]
    print('Running Waiwera')
    run_cmd = ['docker',
               'run',
               '--cidfile', '.cid',
               '--rm',
               '--volume', '{}:{}'.format(current_path, data_path),
               ] + it + work_dir + image + mpiexec + args.waiwera_args
    run_cmd = [c for c in run_cmd if c] # remove empty strings
    print(run_cmd)
    p = subprocess.Popen(run_cmd)
    ret = p.wait()
    with open('.cid', 'r') as f:
        cid = f.readline().strip()[:CID_LEN]
    if ret == 0:
        print('\nWaiwera finished running using Docker container {}.\n'.format(cid))
    else:
        print('\nError running Waiwera in Docker container {}.\n'.format(cid))
    os.remove(".idcheck")
    os.remove('.cid')

signal.signal(signal.SIGINT, signal_handler)

if __name__ == "__main__":
    """
    Args:
    """
    parser = argparse.ArgumentParser(description='Runs Waiwera, \
                        the open-source geothermal flow simulator')
    parser.add_argument('waiwera_args', nargs=argparse.REMAINDER,
                        help='the command passed to waiwera')
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
    waiwera_docker(args)
