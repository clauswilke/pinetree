![pinetree](docs/pinetree-logo.png?raw=true)

# pinetree 
[![Build Status](https://travis-ci.org/benjaminjack/pinetree.svg?branch=master)](https://travis-ci.org/benjaminjack/pinetree)
[![Documentation Status](https://readthedocs.org/projects/pinetree/badge/?version=latest)](http://pinetree.readthedocs.io/en/latest/?badge=latest)

A flexible gene expression simulator with codon-specific translation rates.

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

## Reproducing plots from manuscript

This repository contains scripts to reproduce the simulations and plots from the manuscript that describes Pinetree. R and the R packages `cowplot`, `readr`, `dplyr`, and `stringr` are required to generate plots. Run the following to reproduce the plots from the manuscript:

```
python3 ./examples/three_genes.py
python3 ./examples/three_genes_recoded.py
Rscript plots.R
```

To simulate a bacteriophage T7 infection, run the following script.

```
# WARNING: This simulation takes approximately 2-3 hours to complete
python3 ./examples/phage_model.py
```



