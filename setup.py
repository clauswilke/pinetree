#! /usr/bin/env python3

import os
import re
import sys
import sysconfig
import platform
import subprocess

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand
from setuptools.command.install import install
from setuptools.command.develop import develop
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        self.debug = True
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg, '--target', 'all']

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
        # Move test binary into lib.xxxx directory
        subprocess.check_call(['mv', './pinetree_test', extdir + '/'], cwd=self.build_temp)


class CustomInstall(install):
    """Customized setuptools install command - prints a friendly greeting."""
    def run(self):
        install.run(self)
        # A hacky way of getting the pinetree_test executable to somewhere
        # that the testrunner can find it
        lib_dir = self.install_lib
        scripts_dir = self.install_scripts
        subprocess.check_call(['mv', lib_dir+'/pinetree/pinetree_test', scripts_dir + '/'])

class CatchTestCommand(TestCommand):

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # Run catch tests
        subprocess.check_call(['pinetree_test'])
        print("\nC++ tests complete, now running Python tests...\n")

        # Run python unittest tests
        subprocess.check_call([sys.executable, '-m', 'unittest', 'discover', '--start-directory', 'tests'])

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pinetree',
    version='0.0.1',
    description='a stochastic simulation of gene expression',
    long_description=readme,
    author='Benjamin Jack',
    author_email='benjamin.r.jack@gmail.com',
    url='https://github.com/benjaminjack/pysinthe',
    license=license,
    packages=['pinetree'],
    install_requires=['pyaml',
                      'voluptuous'],
    ext_modules=[CMakeExtension('pinetree/pinetree')],
    cmdclass=dict(install=CustomInstall,
                  build_ext=CMakeBuild,
                  test=CatchTestCommand),
    zip_safe=False,
    scripts = ['bin/pinetree_run.py', 
               'bin/pinetree_batch.py', 
               'bin/parse_genbank.py']
)
