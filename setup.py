#! /usr/bin/env python3

from setuptools import find_packages

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build')
    sys.exit(1)

with open('README.md', encoding='utf-8') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

_VERSION = "0.1.4"

setup(
    name='pinetree',
    version=_VERSION,
    description='stochastic simulation of gene expression with site-specific translation rates',
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Benjamin Jack',
    author_email='benjamin.r.jack@gmail.com',
    url='https://github.com/benjaminjack/pinetree',
    download_url='https://github.com/benjaminjack/pinetree/archive/v' + _VERSION + '.tar.gz',
    license='MIT',
    keywords=['gene', 'codon', 'transcription',
              'translation', 'biology', 'stochastic'],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    zip_safe=False,
    test_suite='tests',
    cmake_args=['-DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.9',
                '-DCMAKE_OSX_ARCHITECTURES:STRING=x86_64'],
    setup_requires=['cmake'],
    install_requires=['cmake']
)
