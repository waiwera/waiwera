import unittest

class TestImport(unittest.TestCase):
    def test_from(self):
        """ expected useage most of the time """
        from pywaiwera import docker
        docker.DockerEnv

    def test_dot(self):
        """ allow user to access submodules with . """
        import pywaiwera
        pywaiwera.__version__
        pywaiwera.docker.DockerEnv
        pywaiwera.common.__version__

        def should_fail():
            pywaiwera.DockerEnv
        self.assertRaises(AttributeError, should_fail)

    def test_star_allowed(self):
        """ from pywaiwera import * should see submodule as current locals """
        c = compile('\n'.join([
            'from pywaiwera import *',
            'docker.DockerEnv',
            ]), 'test', 'exec')
        exec(c)

    def test_star_should_not_see(self):
        """ DockerEnv should be in the scope of pywaiwera.docker"""
        def should_fail():
            from pywaiwera import DockerEnv
        self.assertRaises(ImportError, should_fail)

    def test_star_not_allowed(self):
        """ from pywaiwera import * should see sub-modules only """
        c = compile('\n'.join([
            'from pywaiwera import *',
            'DockerEnv',
            ]), 'test', 'exec')
        try:
            exec(c)
            raise Exception('Should raise NameError, but runs ok.')
        except NameError:
            # expected to fail
            pass

if __name__ == '__main__':
    unittest.main()
