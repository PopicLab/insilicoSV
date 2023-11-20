# insilicoSV
insilicoSV is a software to design and simulate simple and complex structural variants, both novel and known. 

## Requirements  (Prerequisites)
* Python 3.6 and up - [Install](https://www.python.org/downloads/)

## Installation

`$ pip install -r requirements.txt`

## To Run
```
samtools faidx <ref.fna>   # if index file already produced, skip this line
python simulate.py <ref.fna> <par.yaml> <prefix>
```
**NB:** Before running will also need to set `PYTHONPATH` to point to the insilicoSV directory with `export PYTHONPATH={/path/to/insilicoSV}`

## Documentation
For documentation outlining the different features of insilicoSV along with usage examples and data resources, please refer to the wiki:
<!-- toc -->
- [Input guidelines](https://github.com/PopicLab/insilicoSV/wiki#input-guidelines)
- [Example Use Cases](https://github.com/PopicLab/insilicoSV/wiki#example-use-cases)
- [SV Grammar](https://github.com/PopicLab/insilicoSV/wiki/SV-Grammar)
- [Benchmark Genomes](https://github.com/PopicLab/insilicoSV/wiki/Benchmark-Genomes)
- [Automated Pipelines and Additional Utilities](https://github.com/PopicLab/insilicoSV/wiki/Automated-pipelines-and-additional-utilities)
- [Example SV Visualizations](https://github.com/PopicLab/insilicoSV/wiki/Example-SV-visualizations)

## How to Contribute
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change. Please make sure to update tests as appropriate. If you'd like to contribute, please fork the repository and make changes as you'd like.

Steps to contribute:
1. Fork this repository (link to your repository)
2. Create your feature branch (git checkout -b feature/insilicoSV)
3. Commit your changes (git commit -am 'Add some feature')
4. Push to the branch (git push origin feature/insilicoSV)
5. Create a new Pull Request

## Running tests

To run all tests, in the repository root run 
```
./run_tests.sh
```

## Authors
Chris Rohlicek - crohlice@broadinstitute.org

Nick Jiang - nickj@berkeley.edu
