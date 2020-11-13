#!/usr/bin/env python
from __future__ import print_function

# Python 2 and 3 compatible, see: https://stackoverflow.com/a/5868543/2368167
try:
    input = raw_input
except NameError:
    pass

import argparse
import os
import subprocess
import posixpath
import shlex
import sys
import signal
import json
import glob
import re
import time
import shutil

from . import common

REPO = 'waiwera/waiwera'
TAG = 'latest'
WAIWERA_PATH = '/opt/waiwera'
WAIWERA_EXE = '/opt/waiwera/build/waiwera'
VERSION = common.__version__
CID_LEN = 12

### TODO, separate module in new package

def bytes2str(s):
    """ Handle suprocess routines often return bytes instead of str in python3,
    which is different from python2 """
    if hasattr(s, 'decode'):
        return s.decode()
    else:
        return s

def call(cmds, input=None, verbose=False, error_fmt='{}'):
    try:
        p = subprocess.Popen(cmds, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:
        if verbose: print(error_fmt.format(e))
        return False, '', str(e)
    out, err = p.communicate(input=input)
    out, err = bytes2str(out), bytes2str(err)
    retcode = p.wait()
    if retcode == 0:
        return True, out, err
    else:
        if verbose: print(error_fmt.format(err.strip()))
        return False, out, err

def check_tool(name, pth_env, exe_name, cmds, fmt='Found {}', default_pth='',
               verbose=False):
    if verbose: print('Searching for {}...'.format(name))
    if pth_env in os.environ:
        env_dir = os.environ[pth_env]
        fexe = os.path.join(env_dir, exe_name)
        if os.path.isfile(fexe):
            os.environ['PATH'] += os.pathsep + env_dir
            if verbose: print('    Please add {} into your PATH'.format(env_dir))
        else:
            raise Exception('Failed to find {0}'.format(name))
    else:
        fexe = os.path.join(default_pth, exe_name)
        if os.path.isfile(fexe):
            os.environ['PATH'] += os.pathsep + default_pth
            if verbose: print('    Please add {} into your PATH'.format(default_pth))
        else:
            raise Exception('Failed to find {0}'.format(name))
    v = subprocess.check_output(cmds).strip()
    v = bytes2str(v)
    if verbose: print(fmt.format(v))
    return True

class DockerEnv(object):
    """ Provides facilities to do basic health check on docker.  Also have
    additional helper methods to deal with the more fragile Docker Toolbox. """
    def __init__(self, check=True, verbose=True, toolbox=None, toolbox_folder_map=[]):
        """ check (default) is usually recommended, as it will also detect the
        type of docker (toolbox vs non-toolbok).

        For Docker Toolbox users, with check being skipped, toolbox_folder_map
        can be set set to further speed up.  The correct value can be obtained
        by .get_vbox_share().
        """
        super(DockerEnv, self).__init__()
        # None means hasn't been checked
        self.exists = None # True/False after self.check_running()
        self.running = None # True/False after self.check_running()
        self.info = None # {...} after self.check_running()
        self.is_toolbox = toolbox # True/False after self.check_running()
        self.folder_map = toolbox_folder_map # only needed with Docker Toolbox
        if toolbox_folder_map:
            self.is_toolbox = True
        if check:
            # basic check
            self.check_exists()
            self.check_running()
            if not self.running:
                if self.is_toolbox:
                    # basic repair
                    print('Trying to fix environment variables for Docker Toolbox...')
                    self.update_env(True)
                    self.check_running(verbose=True)
                    if not self.running:
                        # more aggressive repair for Docker Toolbox on Windows
                        if sys.platform == 'win32':
                            self.repair_win_toolbox(verbose=verbose)
                            self.check_running(verbose=True)
                        if not self.running:
                            raise Exception('Error, unable to repair Docker Toolbox!')
                else:
                    svr_err = 'n'.join(self.info['ServerErrors'])
                    raise Exception('Error, Docker is not running: {}'.format(svr_err) )
        if self.is_toolbox and len(self.folder_map) == 0:
            self.folder_map = self.get_vbox_share(verbose=True)

    def check_exists(self, verbose=True):
        """ docker --version """
        if verbose: print('Checking if Docker is installed...')
        # only check if docker is reachable, not necessary working
        self.exists, out, err = call(['docker', '--version'],
                            input=None, verbose=verbose,
                            error_fmt='    Cannot find docker command: \n    {0}')
        out, err = bytes2str(out), bytes2str(err)
        if self.exists and verbose:
            print('    {}'.format(out.strip()))
        return self.exists

    def check_running(self, verbose=True):
        """ docker info (only works if daemon/vm/machine also working). If
        works, will be able to tell if Toolbox/Desktop/Linux.
        """
        if self.exists is None:
            self.check_exists()
            if not self.exists:
                # refuse check running if does not exists at all
                return None
        if verbose: print('Checking if Docker is up and running...')
        # docker info only works if the daemon is running
        self.running, out, err = call(['docker', 'info', '--format', '{{json .}}'],
                            input=None, verbose=verbose,
                            error_fmt='    Docker is not working properly: \n    {0}')
        # Docker Desktop/Linux always able to return info,  but Docker Toolbox
        # will fail and return error message
        try:
            self.info = json.loads(out)
        except ValueError:
            self.info = {'STDOUT': out, 'STDERR': err}
            self.is_toolbox = True
        # docker info runs okay, but json may indicate server not running
        if self.running:
            if 'ServerErrors' in self.info:
                # docker info return ok, but docker server not working
                svr_err = 'n'.join(self.info['ServerErrors'])
                if verbose: print('    Docker daemon error: \n     {0}'.format(svr_err))
                self.running = False
            else:
                if verbose: print('    Running on {0}'.format(self.info['OperatingSystem']))
                # determine if running Docker Toolbox
                if 'boot2docker' in self.info['OperatingSystem'].lower():
                    self.is_toolbox = True
                else:
                    # 'docker desktop' or 'linux' can usually be found
                    self.is_toolbox = False
        return self.running

    def _find_vboxmanage(self, verbose=True):
        """ looking for virtualbox manage in VBOX_MSI_INSTALL_PATH or
        VBOX_INSTALL_PATH """
        if sys.platform == 'win32':
            vb = check_tool('VirtualBox: VBoxManage',
                       'VBOX_MSI_INSTALL_PATH', 'VBoxManage.exe',
                       ['VBoxManage', '--version'], '    Found VirtualBox {}',
                       os.path.join(os.environ['ProgramFiles'], 'Oracle', 'VirtualBox'),
                       verbose=verbose)
            if not vb:
                check_tool('VirtualBox: VBoxManage',
                           'VBOX_INSTALL_PATH', 'VBoxManage.exe',
                           ['VBoxManage', '--version'], '    Found VirtualBox {}',
                           os.path.join(os.environ['ProgramFiles'], 'Oracle', 'VirtualBox'),
                           verbose=verbose)
        elif sys.platform == 'darwin':
            msg = 'Please ensure vboxmanage '


    def repair_win_toolbox(self, verbose=True):
        """ Windows + Toolbox only.  Based on the Docker Quick Start bash script
        provided with Docker Toolbox on Windows.

        TODO: Toolbox on Mac should be similar, but most users on Mac should
        have been updated to new enough version with Docker Desktop.
        """
        ## looking for virtualbox manage in VBOX_MSI_INSTALL_PATH or VBOX_INSTALL_PATH
        self._find_vboxmanage()

        ## docker machine name
        vm = self.toolbox_vm()

        ## checking has vm for docker-machine, create one if not
        if verbose: print('Checking if VM {0} exists...'.format(vm))
        vms = subprocess.check_output(['VBoxManage', 'list', 'vms']).strip()
        vms = bytes2str(vms)
        if '"{0}"'.format(vm) not in vms:
            if verbose: print('    VM {0} does not exist, cleanup and re-create one...'.format(vm))
            subprocess.call(['docker-machine', 'rm', '-f', vm])
            # subprocess.call(['rmdir', '/s', '/q', '%HOMEDRIVE%%HOMEPATH%\\.docker\\machine\\machines\\{0}'.format(vm)], shell=True)
            shutil.rmtree('{0}{1}/.docker/machine/machines/{2}'.format(
                              os.environ['HOMEDRIVE'], os.environ['HOMEPATH'], vm),
                              ignore_errors=True)
            os.environ['PROXY_ENV'] = ''
            if 'HTTP_PROXY' in os.environ:
                os.environ['PROXY_ENV'] += '--engine-env HTTP_PROXY=%HTTP_PROXY%'
            if 'HTTPS_PROXY' in os.environ:
                os.environ['PROXY_ENV'] += '--engine-env HTTPS_PROXY=%HTTPS_PROXY%'
            if 'NO_PROXY' in os.environ:
                os.environ['PROXY_ENV'] += '--engine-env NO_PROXY=%NO_PROXY%'
            if 'PROXY_ENV' in os.environ:
                subprocess.call(['docker-machine', 'create', '-d', 'virtualbox'] +
                                os.environ['PROXY_ENV'].split(' ') + [vm])
            else:
                subprocess.call(['docker-machine', 'create', '-d', 'virtualbox', vm])

        ## checking status of the machine, if not running, start docker-machine
        if verbose: print('Checking if docker-machine is running...')
        vmstatus = subprocess.check_output(['docker-machine', 'status', vm]).strip()
        vmstatus = bytes2str(vmstatus)
        if vmstatus != 'Running':
            if verbose: print('    docker-machine is not running, starting...')
            # start
            ok, out, err = call(['docker-machine', 'start'],
                                input=None, verbose=False)
            if not ok:
                if verbose: print(out)
                if verbose: print(err)
                if verbose: print('ERROR! Unable to start docker-machine.')
                exit(1)
            # regenerate-certs
            ok, out, err = call(['docker-machine', 'regenerate-certs', vm],
                                input='y\n', verbose=False)
            if not ok:
                if verbose: print(out)
                if verbose: print(err)
                if verbose: print('ERROR! Unable to start docker-machine.')
                exit(1)
            if verbose: print('    Wait 5 seconds for docker-machine to finish starting up...')
            time.sleep(5)
        ## final checking docker-machine status
        vmstatus = subprocess.check_output(['docker-machine', 'status', vm]).strip()
        vmstatus = bytes2str(vmstatus)
        vmip = subprocess.check_output(['docker-machine', 'ip', vm]).strip()
        vmip = bytes2str(vmip)
        if verbose: print('    Docker is {0}, configured to use the {1} machine with IP {2}'.format(vmstatus, vm, vmip))

        ## this sets current users environment (within this script)
        self.update_env(verbose=verbose)

    def conf_vbox_cpu(self, vm, verbose=True):
        """
        VBoxManage modifyvm %vm% --cpus %NUMBER_OF_PROCESSORS%
        VBoxManage modifyvm default --memory !memory!
        """
        vm = self.toolbox_vm()
        self._find_vboxmanage(verbose=verbose)
        ncpu = multiprocessing.cpu_count()
        subprocess.call(['VBoxManage', 'modifyvm', vm, '--cpus', str(ncpu)])

    def conf_vbox_share(self, vm, to_mount=[], start_after_config=True, verbose=True):
        """ configure VM to use more resources, and shared folders """
        vm = self.toolbox_vm()
        self._find_vboxmanage(verbose=False)
        if verbose: print('Configure shared folders in Docker Toolbox...')
        subprocess.call(['docker-machine', 'start'])

        # TODO help user mount all drives
        # mapping = dict(self.folder_map)
        # to_mount = []
        # for drvltr in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
        #     drv = drvltr+':\\'
        #     if drv not in mapping and os.path.exists(drv):
        #         ymnt = input('    Drive {} not shared, do you want to mount it? (Y/[N])?'.format(drv))
        #         if 'y' in ymnt.lower():
        #             to_mount.append((drvltr, drvltr))

        if to_mount:
            subprocess.call(['docker-machine', 'stop'])
            for drvltr,sname in to_mount:
                hpath = '{0}\\'.format(drvltr)
                # hpath = '{0}\\\\'.format(drvltr)
                # VBoxManage sharedfolder add %1 --name "!lc_drive!" --hostpath "%%I:\\" --automount
                ok, out, err = call(['VBoxManage', 'sharedfolder', 'add', vm,
                    '--name', sname, '--hostpath', hpath,
                    '--automount'],
                    input=None, verbose=True, error_fmt='    {}')
                if not ok:
                    if verbose: print('ERROR! Failed to add share folder in ({0}: --name {1} --hostpath {2}:'.format(
                    vm, sname, hpath))
                    if verbose: print(out + err)
                    raise Exception(out + err)
            if start_after_config:
                if verbose: print('    Restart docker-machine after adding share.')
                subprocess.call(['docker-machine', 'start'])
                self.update_env(verbose=False)

    def get_vbox_share(self, verbose=True):
        """ obtain vm shared drive mapping. VM name can be found in dkr_info['Name']
        """
        def parse_vminfo(vminfo):
            """ extract shared folders from vboxmanage showvminfo default """
            # pat = re.compile("Name: '(?P<name>.+?)', Host path '(?P<host>.+?)' (?P<comment>.+)")
            pat = re.compile("Name: '(?P<name>.+?)', Host path: '(?P<host>.+?)' (?P<comment>.+)")
            maps = []
            found, empty = False, 0
            for line in vminfo.split('\n'):
                if line.startswith('Shared folders:'):
                    if '<none>' in line:
                        break
                    found = True
                    continue
                if found:
                    line = line.strip()
                    if line:
                        ms = pat.search(line)
                        m = ms.groupdict()
                        maps.append((m['host'], m['name']))
                        if verbose: print('    {0} -> {1} ({2})'.format(m['host'], m['name'], m['comment']))
                    else:
                        empty += 1
                    if empty == 2:
                        found = False
            if (not maps) and verbose: print('    <None>')
            return maps
        if verbose: print('Inspecting shared folders (Docker Toolbox)...')
        vm = self.toolbox_vm()
        self._find_vboxmanage(verbose=False)
        vminfo = subprocess.check_output(['vboxmanage', 'showvminfo', vm])
        vminfo = bytes2str(vminfo)
        self.folder_map = parse_vminfo(vminfo)
        return self.folder_map

    def convert_path_vbox_share(self, path):
        """ might return None if not found in mapping, this should work for both tool """
        for m,share in reversed(self.folder_map):
            pp = os.path.join(path, '')  # needed for dirs same prefix
            mm = os.path.join(m, '')
            if sys.platform == 'win32':
                # had a case where drive letter case went funny
                pp = pp.lower()
                mm = mm.lower()
            if pp.startswith(mm):
                # print('{0} matches {1}'.format(pp, mm))
                return '{0}'.format(share) + '/' + path[len(m):]

    def convert_path_nt_to_posix(self, path):
        return '/{0}'.format(path.replace("\\","/").replace(":","").replace('//','/'))

    def toolbox_vm(self):
        """ do I really need this? """
        if self.info is not None:
            if 'Name' in self.info:
                return self.info['Name']
        elif 'DOCKER_MACHINE_NAME' in os.environ:
            if os.environ['DOCKER_MACHINE_NAME']:
                return os.environ['DOCKER_MACHINE_NAME']
        else:
            ok, out, err = call(['docker-machine', 'active', '--timeout', '1'])
            if ok:
                return out
        return 'default'

    def update_env(self, verbose=True):
        """ docker-machine env

        This sets current environment for docker to communicate with docker-
        machine.  NOTE, env only in this script.  Follow on-screen instruction
        if user wish to have the valid env within the current shell.
        """
        if verbose: print('Updating environment variables...')
        ok, out, err = call(['docker-machine', 'env', '--shell=cmd'])
        if not ok:
            print('    {}'.format(err.strip()))
            return
        for line in out.split('\n'):
            if line.startswith('SET'):
                vn, vv = re.split(' |=', line)[1:]
                os.environ[vn] = vv
                if verbose > 1: print("    os.environ['{0}'] = {1}".format(vn, vv))
        ok, out, err = call(['docker-machine', 'env'])
        for line in out.split('\n'):
            if line.startswith('REM ') or line.startswith('#'):
                if verbose: print('    ' + ' '.join(line.split()[1:]))

    def volume_path(self, path=None):
        if path is None:
            path = os.getcwd()
        # toolbox convert vbox mounts
        _cwd = path
        if self.is_toolbox:
            _cwd = self.convert_path_vbox_share(path)
            if _cwd is None and sys.platform == 'win32':
                print('You are running Docker Toolbox which requires directory to be shared.')
                drv, pth = os.path.splitdrive(path)
                ymnt = input('    Do you want drive {} to be shared (Y/[N])?'.format(drv))
                if 'y' in ymnt.lower():
                    drv = drv + '\\'
                    sname = drv[0].upper()
                    self.conf_vbox_share(self.info['Name'], to_mount=[(drv, sname)],
                                           start_after_config=True)
                    self.folder_map = self.get_vbox_share(verbose=True)
                else:
                    return None
            _cwd = self.convert_path_vbox_share(path)
        # windows convert drive and slashes
        cwd = _cwd
        if sys.platform == 'win32':
            cwd = self.convert_path_nt_to_posix(_cwd)
        return cwd

    def run_command(self, current_path, data_path, interactive_cmd, work_dir,
                    image, native_cmd):
        """Returns run command for given run parameters."""
        run_cmd = ['docker',
                   'run',
                   '--cidfile', '.cid',
                   '--rm',
                   '--volume', '{}:{}'.format(current_path, data_path),
                   ] + interactive_cmd + work_dir + image + native_cmd
        run_cmd = [c for c in run_cmd if c] # remove empty strings
        return run_cmd

    def run_ls_test(self, verbose=True):
        """ run docker test if --volume works by list file/folders in current
        directory. """
        if verbose: print('Testing bind mount (--volume) current working directory...')
        cwd = self.volume_path()
        lsout = subprocess.check_output([
            'docker', 'run',
            '--rm',
            '-v', '{0}:/data'.format(cwd),
            'alpine',
            'ls', '-1', '/data',
            ])
        # python2 returns strings, python2 returns bytes
        # in python2 b'' == '', in python3 it matters
        lsout = bytes2str(lsout)
        fs1 = [f.rstrip('\\') for f in lsout.splitlines() if f]
        fs2 = [f for f in glob.glob('*')]
        if set(fs1) != set(fs2):
            msg = str(fs1) + ' != ' + str(fs2)
            msg += "\n    Error bind mount current directory '{}'!".format(cwd)
            if __name__ == '__main__':
                print(msg)
                exit(1)
            raise Exception(msg)
        # print(fs1)
        if verbose: print('    Bind mount OK! Found {0} files in current directory:\n'\
                          '        {1}'.format(len(fs1), '\n        '.join(sorted(list(fs1)))))
        return True

    def run_docker_pull(self, image=None, repo=REPO, tag=TAG):
        """ TODO: should we do an docker images to show what available? """
        if image == None:
            image = ['{0}:{1}'.format(repo, tag)]
        else:
            image = [image]

        print('Checking for Waiwera update')
        pull_cmd = ['docker', 'pull'] + image
        p = subprocess.Popen(pull_cmd)
        ret = p.wait()

    def run_copy_examples(self, image=None, repo=REPO, tag=TAG, verbose=False):
        current_path = self.volume_path()
        data_path = '/data'

        if image == None:
            image = ['{0}:{1}'.format(repo, tag)]
        else:
            image = [image]

        script_name = '.copy_examples.py'
        work_dir = ['--workdir', data_path]
        cmd = ['/usr/bin/python', script_name]
        # cmd = ['pwd']

        fo = open(".idcheck", "wb")
        fo.close()

        # TODO: replace the .sh sript with this .py one, copies only some files
        py_lines = '\n'.join([
            "import os",
            "import glob",
            "import shutil",
            "DEST_PATH = '/data/examples'",
            "SOURCE_PATH = '/opt/waiwera/test/benchmark'",
            "COPY_TYPES = ['g*.dat', '*.json', 'g*.exo', 'g*.msh']",
            "def md(pth):",
            "    if not os.path.exists(pth):",
            "        os.mkdir(pth)",
            "for root, dirs, files in os.walk(SOURCE_PATH):",
            "    if os.path.basename(root) not in ['run', 'data']:",
            "        to_dir = root.replace(SOURCE_PATH, DEST_PATH)",
            "        md(to_dir)",
            "        if 'run' in dirs:",
            "            for fn in COPY_TYPES:",
            "                for f in glob.glob(os.path.join(root, 'run', fn)):",
            "                    print('Copying file: {}'.format(f))",
            "                    shutil.copy2(f, to_dir)",
            ])

        with open(script_name, 'w') as fo:
            fo.write(py_lines)

        it = ['']
        run_cmd = self.run_command(current_path, data_path, it, work_dir,
                    image, cmd)
        if verbose: print('Docker command:', run_cmd)
        print('Running Docker...')
        # TODO: window+git bash+toolbox need shell=True to handle path with space
        p = subprocess.Popen(run_cmd)
        ret = p.wait()
        with open('.cid', 'r') as f:
            cid = f.readline().strip()[:CID_LEN]
        if ret == 0:
            print('\nFinished running using Docker container {}.\n'.format(cid))
        else:
            print('\nError copying examples from Docker container {}.\n'.format(cid))
        os.remove(".idcheck")
        os.remove('.cid')
        os.remove(script_name)

    def run_waiwera(self, filename='', waiwera_args=[], image=None,
                    repo=REPO, tag=TAG, num_processes=None, interactive=False,
                    noupdate=False, verbose=False, dryrun=False, get_output=False):
        """ run waiwera

        If interactive is True, docker withh run with --iteractive --tty using
        waiwera_args as command to execute. If waiwera_args is left empty, a
        simple '/bin/bash' will be used.

        If dryrun is True, final command will not be run but returned instead.
        This is useful for testing.

        If get_output is True, the STDOUT from docker run will be returned as
        string.  Note this option has no effect when dryrun is True (the actual
        command will be returned instead).
        """
        current_path = self.volume_path()
        data_path = '/data'
        waiwera_args = [a for a in waiwera_args if a] # remove empty strings

        if image == None:
            image = ['{0}:{1}'.format(repo, tag)]
        else:
            image = [image]

        if num_processes:
            np = ['-np', '{}'.format(num_processes)]
        else:
            np = ['']

        if not noupdate:
            print('Checking for Waiwera update')
            pull_cmd = ['docker', 'pull'] + image
            p = subprocess.Popen(pull_cmd)
            ret = p.wait()

        if interactive:
            it = ['--interactive', '--tty']
            work_dir = ['']
            if len(waiwera_args) == 0:
                native_cmd = ['/bin/bash']
            else:
                native_cmd = waiwera_args
        else:
            it  = ['']
            work_dir = ['--workdir', data_path]
            native_cmd = ['mpiexec'] + np + [WAIWERA_EXE] + [filename] + waiwera_args

        #  docker run -v ${p}:/data -w /data waiwera-phusion-debian mpiexec -np $args[1] /home/mpirun/waiwera/dist/waiwera $args[0]
        run_cmd = self.run_command(current_path, data_path, it, work_dir,
                                   image, native_cmd)
        if verbose:
            print('Docker command:', run_cmd)
        if dryrun:
            return run_cmd

        fo = open(".idcheck", "wb")
        fo.close()
        print('Running Docker...')
        # TODO: window+git bash+toolbox need shell=True to handle path with space
        if get_output:
            p = subprocess.Popen(run_cmd, stdout=subprocess.PIPE)
            ret = p.wait()
            output = p.stdout.read()
        else:
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
        if get_output:
            return output

def in_directory(file, directory):
    """ checks if a file is within a directory (or its sub-directories).

    https://stackoverflow.com/q/3812849/2368167
    make both absolute, true if the common prefix of both is equal to
    directory e.g. /a/b/c/d.rst and directory is /a/b, the common prefix is
    /a/b
    """
    directory = os.path.join(os.path.realpath(directory), '')
    file = os.path.realpath(file)
    return os.path.commonprefix([file, directory]) == directory

def input_file_ok(f):
    """ ensures input.json ok:
        1. may or may not be ok, but decided to disencourage users to use abs path
        2. has to be in sub-folder because Docker --volume bind mount
        3. has to be linux path sep
    """
    ok = True
    if os.path.isfile(f):
        if os.path.isabs(f):
            print('Error, waiwera-dkr (using Docker) does not support absolute file path!')
            ok = False
        if not in_directory(f, os.getcwd()):
            print('Error, waiwera-dkr (using Docker) requires input file to be contained inside current (or sub-directory) directory!')
            ok = False
        if '\\' in f:
            print('Error, waiwera-dkr (using Docker) requires linux-style relative path!')
            ok = False
    else:
        print('Error, simulation input file not found: {0}'.format(f))
        ok = False
    return ok

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

signal.signal(signal.SIGINT, signal_handler)

def main():
    """
    Args:
    """
    examples = "\n".join([
        "examples:",
        "  waiwera-dkr input.json",
        "  waiwera-dkr -np 2 input.json",
        "  waiwera-dkr --tag testing input.json",
        "  waiwera-dkr --interactive",
        ])
    parser = argparse.ArgumentParser(description='Runs Waiwera, the parallel open-source geothermal flow simulator, via Docker',
                        epilog=examples, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('filename', metavar='FILENAME', nargs='?',
                        help='Waiwera input file name (must come after optional arguments)',
                        default='')
    parser.add_argument('waiwera_args', metavar='...', nargs=argparse.REMAINDER,
                        help='additional arguments passed to Waiwera (must come last)')

    # without default=1, mpiexec will utilise max num of processes within
    # container, this can cause slow down small models, and crash very small
    # models from examples
    parser.add_argument('-np', '--num_processes', help='the number of \
                        processors to utilize, default is 1 (serial)',
                        default=1)

    parser.add_argument('-r', '--repo',
                        default=REPO)
    parser.add_argument('-t', '--tag',
                        default=TAG)
    parser.add_argument('-i', '--image', help='the docker image to use \
                        e.g. waiwera/waiwera:latest')
    parser.add_argument('-it','--interactive',
                        help='starts an interactive terminal (and does not run Waiwera)',
                        action='store_true')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-nu','--noupdate',
                    help='do not check for an updated Waiwera image before running',
                    action='store_true')
    group.add_argument('-u','--update',
                    help='pull an updated Waiwera image and exit',
                    action='store_true')

    parser.add_argument('-tv','--test_volume',
                    help='test docker --volume (bind mount) with current directory and exit',
                    action='store_true')
    parser.add_argument('-v','--verbose',
                    help='print additional diagnostic messages while running',
                    action='store_true')
    parser.add_argument('-e','--examples',
                    help='create example models (./examples) at current directory and exit',
                    action='store_true')
    parser.add_argument('--version',
                    help='show PyWaiwera version, together with Waiwera version',
                    action='store_true')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # early check (faster) on inputfile
    if not args.interactive:
        if args.filename:
            inok = input_file_ok(args.filename) # first one always input.json
            if not inok:
                exit(1)

    if args.interactive:
        # interactive assumes remainder arguments are commands to run container
        #   waiwera-dkr -it
        #     => docker run --interactive --tty waiwera/waiwera:latest /bin/bash
        #   waiwera-dkr -it touch x
        #     => docker run --interactive --tty waiwera/waiwera:latest touch x
        if args.filename:
            args.waiwera_args = [args.filename] + args.waiwera_args
            args.filename = ''


    dkr = DockerEnv(check=True)
    # print('docker.exist', dkr.exists, 'dkr.running', dkr.running, 'dkr.is_toolbox', dkr.is_toolbox)
    if args.test_volume:
        dkr.run_ls_test()
        exit(0)

    if args.version:
        accept_kws = ['waiwera_args', 'image', 'repo', 'tag',
                      'num_processes', 'noupdate', 'verbose']
        kw_run_waiwera = {k:v for k,v in vars(args).items() if k in accept_kws}
        kw_run_waiwera['waiwera_args'] = ['--version']
        kw_run_waiwera['get_output'] = True
        v = dkr.run_waiwera(**kw_run_waiwera)
        if args.image == None:
            image = '{0}:{1}'.format(args.repo, args.tag)
        else:
            image = args.image
        print('PyWaiwera {}'.format(common.__version__))
        print('Waiwera {} (Docker image {})'.format(v.strip(), image))
        exit(0)

    if args.update:
        accept_kws = ['image', 'repo', 'tag']
        kws = {k:v for k,v in vars(args).items() if k in accept_kws}
        dkr.run_docker_pull(**kws)
        exit(0)

    if args.examples:
        accept_kws = ['image', 'repo', 'tag', 'verbose']
        kws = {k:v for k,v in vars(args).items() if k in accept_kws}
        dkr.run_copy_examples(**kws)
        exit(0)

    if args.filename or args.interactive:
        accept_kws = ['filename', 'waiwera_args', 'image', 'repo', 'tag',
                      'num_processes', 'interactive', 'noupdate', 'verbose']
        kw_run_waiwera = {k:v for k,v in vars(args).items() if k in accept_kws}
        dkr.run_waiwera(**kw_run_waiwera)

