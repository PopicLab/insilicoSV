# insilicoSV: a framework for structural variant simulation 

## Overview

insilicoSV generates synthetic diploid genome sequences, given a reference genome and a configuration file.

It supports the following functionality:

* 19 types of structural variants (simple and complex), indels, and SNPs
* a grammar to define custom structural rearrangement signatures
* random, context-aware (e.g., in repeat regions), or fixed-mode genome placement

## Installation

Prerequisite: Python 3.8 and up - [Install](https://www.python.org/downloads/)

Installation using pip:

* `$ pip install insilicosv`

<!-- Installation using conda:

* Install and configure [bioconda](https://bioconda.github.io/)
* Install insilicosv with `conda install insilicosv` -->

## To Run
The recommended workflow for running insilicoSV is as follows:
1. Create a new directory
2. Populate a simulation config file and place it in the directory
3. Run insilicoSV providing the config file as input:
```
insilicosv <config_yml>
```
4. Results will be produced in this directory

The config file syntax is detailed in the [Input guidelines](docs/input_guidelines.md) section of the
documentation.

## Documentation
For documentation outlining the different features of insilicoSV along with usage examples and data resources, please refer to the following sections:
<!-- toc -->
- [Input guidelines](docs/input_guidelines.md)
- [Example Use Cases](docs/example_use_cases.md)
- [SV Grammar](docs/sv_grammar.md)
- [Example SV Visualizations](docs/example_sv_visualizations.md)
- [Tutorial Jupyter notebook](docs/demo_notebook.md)


## Authors
Chris Rohlicek - crohlice@broadinstitute.org

Nick Jiang - nickj@berkeley.edu

Ilya Shlyakhter - ilya@broadinstitute.org

Victoria Popic - vpopic@broadinstitute.org
