# insilicoSV (WIP): a framework for structural variant simulation 

## Overview

insilicoSV generates synthetic diploid genome sequences, given a reference genome and a configuration file.

It supports the following functionality:

* 19 types of structural variants (simple and complex), indels, and SNPs
* a grammar to define custom structural rearrangement signatures
* random, context-aware (e.g., in repeat regions), or fixed-mode genome placement

## Installation

Prerequisite: Python 3.6 and up - [Install](https://www.python.org/downloads/)

* Clone the repository: ```git clone git@github.com:PopicLab/insilicoSV.git ```

* `$ pip install -r install/requirements.txt`

* Set the ```PYTHONPATH``` as follows: ```export PYTHONPATH=${PYTHONPATH}:/path/to/insilicoSV```


## To Run
```
python insilicosv/simulate.py <config.yaml> 
```

## Documentation
For documentation outlining the different features of insilicoSV along with usage examples and data resources, please refer to the wiki:
<!-- toc -->
- [Input guidelines](https://github.com/PopicLab/insilicoSV/wiki#input-guidelines)
- [Example Use Cases](https://github.com/PopicLab/insilicoSV/wiki#example-use-cases)
- [SV Grammar](https://github.com/PopicLab/insilicoSV/wiki/SV-Grammar)
- [Benchmark Genomes](https://github.com/PopicLab/insilicoSV/wiki/Benchmark-Genomes)
- [Automated Pipelines and Additional Utilities](https://github.com/PopicLab/insilicoSV/wiki/Automated-pipelines-and-additional-utilities)
- [Example SV Visualizations](https://github.com/PopicLab/insilicoSV/wiki/Example-SV-visualizations)


## Authors
Chris Rohlicek - crohlice@broadinstitute.org

Nick Jiang - nickj@berkeley.edu
