# Demo Jupyter notebook

We provide a tutorial Jupyter notebook to demonstrate the workflow of
using insilicoSV to generate a synthetic diploid genome and simulating
reads for downstream analysis. In the notebook, we provide an example
configuration file including a small set of simple and complex
variants, and we include calls to DWGSIM and PBSIM3 to generate short
and long reads respectively. Additionally, we include utilities to
parse insilicoSV output and plot the size distributions of the simulated SVs, as well as generate IGV pileup images of the aligned
reads at the sites of the simulated variants.

To run the notebook, first create a conda environment with all dependencies
used in the notebook.  The environment definition is provided in
`workflows/insilicosv-demo-env.reqs.txt`.  After setting up
[bioconda](https://bioconda.github.io/), create the environment with the command

```
conda create -n insilicosv-demo-env --file workflows/insilicosv-demo-env.reqs.txt
```

Activate the environment with
```
conda activate insilicosv-demo-env
```

To ensure that the current insilicosv version is installed, run
```
pip install --upgrade insilicosv
```

and launch the notebook with `jupyter notebook workflows/demo.ipynb` .
