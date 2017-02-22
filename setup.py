# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pysinthe',
    version='0.1.0',
    description='a stochastic simulation of gene expression',
    long_description=readme,
    author='Benjamin Jack',
    author_email='benjamin.r.jack@gmail.com',
    url='https://github.com/benjaminjack/pysinthe',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
