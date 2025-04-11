The WDL pipeline, ```insilicosv_workflow.wdl``` (in the ```workflows/``` directory), is provided to streamline 
SV truthset generation. It supports: (1) genome simulation using ```insilicoSV``, (2) read simulation for 
Illumina, PacBio, and ONT platforms (as well as custom long read simulation), (3) read alignment, (4) and SV visualization.

To install all the pipeline requirements:
* Install and configure [bioconda](https://bioconda.github.io/)
* Create a new conda environment: `$> conda create -n insilicosv_wdl --file workflows/install/insilicosv_workflow_requirements.txt `
* Activate the environment: `$> conda activate insilicosv_wdl`
* Install ```insilicoSV```: `$> pip install .`

To run the pipeline:
* Create an output folder: `$> mkdir test_truthset`
* Create and populate an `insilicoSV` YAML config file inside `test_truthset`
* Copy the JSON input file: `$> cp workflows/configs/wdl_config.json test_truthset/test.json`
* Populate all the input parameters in `test_truthset/test.json` (set `insilicosv_workflow.configYAML` to the `insilicoSV` YAML file)
* Run the pipeline: `$> cromwell run -i test_truthset/test.json workflows/insilicosv_workflow.wdl`

Pipeline parameters (JSON):
* `insilicosv_workflow.configYAML` [required]:  path to insilicoSV YAML config
* `insilicosv_workflow.hapCoverage` [required]: coverage per haplotype
* `insilicosv_workflow.short` [required]: generate Illumina paired-end reads (using DWGSIM)
* `insilicosv_workflow.hifi` [required]: generate PacBio HiFi reads (using PBSIM3)
* `insilicosv_workflow.ont` [required]: generate ONT reads (using PBSIM3)
* `insilicosv_workflow.ref` [required]: path to reference FASTA
* `insilicosv_workflow.threads` [required]: number of threads to use for read alignment and sorting
* `insilicosv_workflow.customLrPreset` [optional]: custom PBSIM3 parameters, will generate a custom long-read dataset if set
* `insilicosv_workflow.customLrAlnPreset` [optional]: custom minimap2 parameters to use with the custom long-read dataset
* `insilicosv_workflow.igvArgs` [optional]: custom igv-reports parameters (overwrites the default setting)
