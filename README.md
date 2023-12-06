# insilicoSV (WIP): a framework for structural variant simulation 

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

Installation using conda:

* Install and configure [bioconda](https://bioconda.github.io/)
* Install insilicosv with `conda install insilicosv`

## To Run
```
insilicosv <config_yml>
```

The config file syntax is detailed below.  Outputs go into the directory where the config file is stored.

## Documentation
For documentation outlining the different features of insilicoSV along with usage examples and data resources, please refer to the wiki:
<!-- toc -->
- [Input guidelines](docs/input_guidelines.md)
- [Example Use Cases](docs/example_use_cases.md)
- [SV Grammar](docs/sv_grammar.md)
- [Benchmark Genomes](docs/benchmark_genomes.md)
- [Additional Utilities](docs/automated_pipelines_and_additional_utilities.md)
- [Example SV Visualizations](docs/example_sv_visualizations.md)
- [Tutorial Jupyter notebook](docs/demo_notebook.md)


## Authors
Chris Rohlicek - crohlice@broadinstitute.org

Nick Jiang - nickj@berkeley.edu

Ilya Shlyakhter - ilya@broadinstitute.org

Victoria Popic - vpopic@broadinstitute.org
