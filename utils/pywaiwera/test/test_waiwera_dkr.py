import unittest

class TestRunWaiwera(unittest.TestCase):
    def setUp(self):
        """ expected useage most of the time """
        from pywaiwera import docker
        self.env = docker.DockerEnv(check=False, verbose=False,)

    def test_run_waiwera(self):
        # waiwera-dkr input.json => docker run mpiexec -np 1 waiwera input.json
        args = {
            'filename': 'input.json',
            'num_processes': 1,
            'interactive': False,
            'noupdate': True,
            'dryrun': True,
        }
        cmd = self.env.run_waiwera(**args)
        self.assertEqual(cmd[:6],
                         ['docker', 'run', '--cidfile', '.cid',
                          '--rm', '--volume'])
        # skip the 7th (volume binding depends in os and docker environment)
        self.assertEqual(cmd[7:],
                         ['--workdir', '/data', 'waiwera/waiwera:latest',
                          'mpiexec', '-np', '1', '/opt/waiwera/build/waiwera',
                          'input.json'])

        # waiwera-dkr input.json => docker run mpiexec -np 2 waiwera input.json a b
        args = {
            'filename': 'input.json',
            'waiwera_args': ['a', 'b'],
            'num_processes': 2,
            'interactive': False,
            'noupdate': True,
            'dryrun': True,
        }
        cmd = self.env.run_waiwera(**args)
        self.assertEqual(cmd[:6],
                         ['docker', 'run', '--cidfile', '.cid',
                          '--rm', '--volume'])
        # skip the 7th (volume binding depends in os and docker environment)
        self.assertEqual(cmd[7:],
                         ['--workdir', '/data', 'waiwera/waiwera:latest',
                          'mpiexec', '-np', '2', '/opt/waiwera/build/waiwera',
                          'input.json', 'a', 'b'])

    def test_img_repo(self):
        # waiwera-dkr --tag testing input.json => docker run waiwera/waiwera:testing mpiexec -np 1 waiwera input.json
        args = {
            'filename': 'input.json',
            'num_processes': 1,
            'interactive': False,
            'noupdate': True,
            'dryrun': True,
            'tag': 'testing',
        }
        cmd = self.env.run_waiwera(**args)
        self.assertEqual(cmd[:6],
                         ['docker', 'run', '--cidfile', '.cid',
                          '--rm', '--volume'])
        # skip the 7th (volume binding depends in os and docker environment)
        self.assertEqual(cmd[7:],
                         ['--workdir', '/data', 'waiwera/waiwera:testing',
                          'mpiexec', '-np', '1', '/opt/waiwera/build/waiwera',
                          'input.json'])

        # waiwera-dkr input.json => docker run wai/wai:latestmpiexec -np 1 waiwera input.json
        args = {
            'filename': 'input.json',
            'num_processes': 1,
            'interactive': False,
            'noupdate': True,
            'dryrun': True,
            'tag': 'testing', # overwriten by image
            'image': 'wai/wai:latest',
        }
        cmd = self.env.run_waiwera(**args)
        self.assertEqual(cmd[:6],
                         ['docker', 'run', '--cidfile', '.cid',
                          '--rm', '--volume'])
        # skip the 7th (volume binding depends in os and docker environment)
        self.assertEqual(cmd[7:],
                         ['--workdir', '/data', 'wai/wai:latest',
                          'mpiexec', '-np', '1', '/opt/waiwera/build/waiwera',
                          'input.json'])

    def test_interactive(self):
        #   waiwera-dkr -it
        #     => docker run --interactive --tty waiwera/waiwera:latest /bin/bash
        args = {
            'filename': None,
            'waiwera_args': [],
            'num_processes': 2,
            'interactive': True,
            'noupdate': True,
            'dryrun': True,
        }
        cmd = self.env.run_waiwera(**args)
        self.assertEqual(cmd[:6],
                         ['docker', 'run', '--cidfile', '.cid',
                          '--rm', '--volume'])
        # skip the 7th (volume binding depends in os and docker environment)
        self.assertEqual(cmd[7:],
                         ['--interactive', '--tty', 'waiwera/waiwera:latest',
                          '/bin/bash'])

        #   waiwera-dkr -it touch x
        #     => docker run --interactive --tty waiwera/waiwera:latest touch x
        args = {
            'filename': None,
            'waiwera_args': ['python', '--version'],
            'num_processes': 2,
            'interactive': True,
            'noupdate': True,
            'dryrun': True,
        }
        cmd = self.env.run_waiwera(**args)
        self.assertEqual(cmd[:6],
                         ['docker', 'run', '--cidfile', '.cid',
                          '--rm', '--volume'])
        # skip the 7th (volume binding depends in os and docker environment)
        self.assertEqual(cmd[7:],
                         ['--interactive', '--tty', 'waiwera/waiwera:latest',
                          'python', '--version'])


if __name__ == '__main__':
    unittest.main()
