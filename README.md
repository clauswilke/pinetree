![pinetree](docs/pinetree-logo.png?raw=true)

# pinetree 
[![Build Status](https://travis-ci.org/benjaminjack/pinetree.svg?branch=master)](https://travis-ci.org/benjaminjack/pinetree)
[![Documentation Status](https://readthedocs.org/projects/pinetree/badge/?version=latest)](http://pinetree.readthedocs.io/en/latest/?badge=latest)

A flexible gene expression simulator with codon-specific translation rates.

**WARNING:** Pinetree is under active development and not ready for public consumption. Use at your own risk.

## Requirements

Pinetree requires Python, CMake, and a modern C++ compiler. Python 3 is recommended.

## Installation

To install the latest stable version of pinetree from PyPI, run the following:

```
pip install pinetree 
```

The latest development build may be installed from GitHub as follows:

```   
git clone https://github.com/benjaminjack/pinetree.git
pinetree/setup.py install
```

## Documentation

Full documentation is available [here](http://pinetree.readthedocs.io/).

You may also build the documentation from the source code. Building the documentation requires sphinx.

```
pinetree/setup.py build_sphinx
```
