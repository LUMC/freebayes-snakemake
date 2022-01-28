[![Continuous Integration](https://github.com/LUMC/freebayes-snakemake/actions/workflows/ci.yml/badge.svg)](https://github.com/LUMC/freebayes-snakemake/actions/workflows/ci.yml)
[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)
![GitHub release](https://img.shields.io/github/v/release/LUMC/freebayes-snakemake)
![Commits since latest release](https://img.shields.io/github/commits-since/LUMC/freebayes-snakemake/latest)

# freebayes-snakemake
Example of a snakemake project

## Installation
Download the repository from github
```bash
git clone https://github.com/LUMC/freebayes-snakemake.git
```

Install and activate the
[conda](https://docs.conda.io/en/latest/miniconda.html)
environment.
```bash
conda env create --file environment.yml
conda activate freebayes-snakemake
```

## Settings priority
There are three levels where configuration options are set, in decreasing order
of priority.
1. Flags passed to snakemake using `--config`, or in the specified
   `--configfile`.
2. Setting specified in the PEP project configuration, under the key
   `freebayes-snakemake`
3. The default settings for the pipeline, as specified in the `common.smk` file

## Pipeline settings
The pipeline only requires a single PEP configuration file, which specifies a
csv [sample
table](https://github.com/LUMC/freebayes-snakemake/blob/main/tests/pep/samples.csv) and a reference
([example](https://github.com/LUMC/freebayes-snakemake/blob/main/tests/pep/project_config.yaml)).

If you have multiple read pairs per sample, you can also specify a [subsample
table](https://github.com/LUMC/freebayes-snakemake/blob/main/tests/pep/subsamples.csv),
with one line for each read pair
([example](https://github.com/LUMC/freebayes-snakemake/blob/main/tests/pep/subsamples.csv)).
