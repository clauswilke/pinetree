#! /usr/bin/env python3

from skbuild import setup

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
    package_dir={'pinetree': 'src/pinetree'},
    zip_safe=False,
    test_suite='tests',
    entry_points={
        'console_scripts': [
            'pinetree_test=pinetree:pinetree_test'
        ]
    }
)
