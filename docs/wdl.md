### Pipeline 1: basic single-genome simulation

The WDL pipeline, `insilicosv_basic_workflow.wdl` (in the `workflows/` directory), is provided to streamline 
SV truthset generation. It supports: (1) single genome simulation using `insilicoSV`, (2) read simulation for 
Illumina, PacBio, and ONT platforms (as well as custom long read simulation), (3) read alignment, (4) and SV visualization.

To install all the pipeline requirements:
* Install and configure [bioconda](https://bioconda.github.io/)
* Create a new conda environment: `$> conda create -n insilicosv_wdl --file workflows/install/insilicosv_workflow_requirements.txt `
* Activate the environment: `$> conda activate insilicosv_wdl`
* Install ```insilicoSV```: `$> pip install .`

To run the pipeline:
* Create an output folder: `$> mkdir test_truthset`
* Create and populate an `insilicoSV` YAML config file inside `test_truthset` (for example, see ```workflows/configs/demo_insilicosv_config.wdl```)
* Copy the JSON input file: `$> cp workflows/configs/basic_wdl_config.json test_truthset/test.json`
* Populate all the input parameters in `test_truthset/test.json` (set `insilicosv_basic.configYAML` to the `insilicoSV` YAML file)
* Run the pipeline: `$> cromwell run -i test_truthset/test.json workflows/insilicosv_basic_workflow.wdl`

Pipeline parameters (JSON):
* `insilicosv_basic.configYAML` [required]:  path to insilicoSV YAML config
* `insilicosv_basic.hapCoverage` [required]: coverage per haplotype
* `insilicosv_basic.short` [required]: generate Illumina paired-end reads (using DWGSIM)
* `insilicosv_basic.hifi` [required]: generate PacBio HiFi reads (using PBSIM3)
* `insilicosv_basic.ont` [required]: generate ONT reads (using PBSIM3)
* `insilicosv_basic.ref` [required]: path to reference FASTA
* `insilicosv_basic.threads` [required]: number of threads to use for read alignment and sorting
* `insilicosv_basic.customLrPreset` [optional]: custom PBSIM3 parameters, will generate a custom long-read dataset if set
* `insilicosv_basic.customLrAlnPreset` [optional]: custom minimap2 parameters to use with the custom long-read dataset
* `insilicosv_basic.igvArgs` [optional]: custom igv-reports parameters (overwrites the default setting)

### Pipeline 2: lineage-based (multi-)genome simulation

To simulate genome evolution across multiple time points or a mixture of related genomes resulting from a branched 
lineage tree (e.g., cancer subclones), we provide a configurable WDL workflow which can be automatically generated using
the ```workflows/generate_lineage_wdl.py``` script based on the desired maximum lineage tree depth.

To generate the lineage WDL workflow (for an example depth of N=4):
* Install all the pipeline requirements as described in the section above and activate the created conda environment
* cd workflows/
* python generate_lineage_wdl.py 4

The resulting WDL pipeline ```insilicosv_lineage_workflow_4.wdl``` will be generated in the ```workflows/``` directory. 

To run the pipeline (from the top-level `insilicosv` directory):
* Create an output folder: `$> mkdir lineage_out` 
* Create and populate multiple `insilicoSV` YAML config files inside `lineage_out` - one for each time point (for example, see ```workflows/configs/demo_insilicosv_config.wdl```)
* Copy the JSON input file: `$> cp workflows/configs/lineage_wdl_config.json lineage_out/inputs.json`
* Populate all the input parameters in `lineage_out/inputs.json` 
* Run the pipeline: `$> cromwell run -i lineage_out/inputs.json workflows/insilicosv_lineage_workflow_4.wdl`

Pipeline parameters (JSON):
* `GenomeMix.outDir` [required]: path to the output directory
* `GenomeMix.reference` [required]: path to the reference genome FASTA (e.g. hg38)
* `GenomeMix.configs` [required]: list of `insilicoSV` YAML config files associated with each time point 
* `GenomeMix.genomes` [required]: list of genomes represented as an ordered list of time points (e.g. [[0], [0, 1, 3, 4], [0, 2]]); captures the desired lineage tree
* `GenomeMix.totalCoverage` [optional]: total coverage of the generated mixed-genome read dataset
* `GenomeMix.prevalence` [optional]: list of cellular prevalence values for each genome (e.g. [0.7, 0.1, 0.2]), must add to 1
* `GenomeMix.short` [optional]: generate Illumina paired-end reads (using DWGSIM)
* `GenomeMix.hifi` [optional]: generate PacBio HiFi reads (using PBSIM3)
* `GenomeMix.ont` [optional]: generate ONT reads (using PBSIM3)
* `GenomeMix.threads` [optional]: number of threads to use for read alignment and sorting

Note: setting `totalCoverage` > 0 will trigger read simulation and alignment; in this case the `prevalence` parameter is required, 
as well as at least one sequencing platform should be enabled by setting its corresponding flag to True.

Note: each intermediate time point will result in the simulation of a new genome using `insilicoSV` and the corresponding YAML file; 
however, only the leaf/terminal genomes listed in the `genomes` parameter will be mixed in the final sample using the 
ratios configured with the `prevalence` parameter.


