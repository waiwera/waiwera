import argparse
import os
import subprocess
import posixpath
import shlex
import sys

REPO = 'waiwera/waiwera'
TAG = 'latest'
WAIWERA_PATH = '/opt/waiwera'

def convert_path_nt_to_posix(path):
    return '/{0}'.format(path.replace("\\","/").replace(":",""))

##########
# Author: Tim Harton
# Date: 4/11/19
# A simple python wrapper around the docker image for waiwera
#########

if __name__ == "__main__":
    """
    Args:
    """
    parser = argparse.ArgumentParser(description='Runs Waiwera, the open-source geothermal flow simulator.')
    parser.add_argument('command', nargs='+', help='the command passed to waiwera')
    parser.add_argument('-np', '--num_processors', help='The number of'\
        ' processors to utilize, otherwise uses the docker default for'\
        ' your system')
    parser.add_argument('-r', '--repo', default=REPO)
    parser.add_argument('-t', '--tag', default=TAG)
    parser.add_argument('-i', '--image', help='The docker image to use'
                        'e.g. waiwera/waiwera:latest')

    args = parser.parse_args()

    current_path = os.getcwd()
    if sys.platform == 'win32':
        current_path = convert_path_nt_to_posix(current_path)

    data_path = posixpath.join(WAIWERA_PATH, 'data')

    if args.image == None:
        image = '{0}:{1}'.format(args.repo, args.tag)
    else:
        image = args.image

    if args.num_processors:
        np ='-np {0}'.format(args.num_processors)
    else:
        np = ''

    command = ''
    for string in args.command:
        command = '{0} {1}'.format(command, string)

    pull_cmd = "docker pull {0}".format(image)
    output = subprocess.check_output(shlex.split(pull_cmd))

    run_cmd = "docker run --rm --volume {0}:{1} --workdir {1} {2} mpiexec {3} {4} {5}" \
            .format(current_path,
                    data_path,
                    image,
                    np,
                    WAIWERA_PATH,
                    command)

    output = subprocess.check_output(shlex.split(run_cmd))
