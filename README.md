# SIEVE_benchmark_pipeline

This repository harbours the benchmarking framework, built upon [Snakemake](https://snakemake.readthedocs.io/en/stable/), for [SIEVE](https://github.com/szczurek-lab/SIEVE). 

## Required packages

This framework contains scripts written in Python 3 and R. Apart from an environment containing snakemake (we recommend using Conda), the following packages should be installed before running the pipeline:

- For Python3:
  - paramiko >= 2.8.0
- For R:
  - base >= 4.0
  - stringr
  - scales
  - dplyr
  - optparse
  - ape
  - phangorn
  - tools

## Configuration

A few files should be configured before running.

1. The pipeline is mainly configured in `config.yaml`. There are some entries to be filled by users, which are marked by a phrase `# TO BE SET` , followed either by `[MANDATORY]` (must be set) or by `[OPTIONAL]` (could be ignored). Users can search for the phrase to set everything up efficiently.

   - The template configuration files for SIEVE are under `templates/`.

   - For the key `[configFiles][SIEVE_simulator]`, a configuration file for the data simulator [SIEVE_simulator](https://github.com/szczurek-lab/SIEVE_simulator) should be specified. Those simulated scenarios used in the SIEVE paper are listed in `simulation_configs/`. For details, check the paper.

2. In `run/run.sh`, users can set a few things, e.g., the name of the conda environment containing snakemake (by default, `snakemake`), the number of cores to use and their ranges, etc.

### Advanced configuration

Since SiFit requires a large amount of memory even working on a small dataset, the framework supports running SiFit alone on another server (referred to as `the remote server`) with the help of git. For this to work, a few things should be noted and configured:

- The machine you plan to run the pipeline (referred to as `the local machine`) and `the remote server` should meet one of the following conditions:
  - They are in the same local network. 
  - If they are not in the same local network, `the local machine` must have a public IP address for `the remote server` to clone the git repo. However, `the remote server` can behind a gate server with a public IP address, specified by the key `[servers][jumpServerOfRemote]` in `config.yaml`.
- In `Snakefile`, comment out `include: "sifit_local.snake"`, and uncomment `# include: "sifit_remote.snake"`. 
- Set up `run/run_remote_sifit_true_monovar.sh` similarly to `run/run.sh`.

## Run the pipeline

With everything set up, users can run the pipeline under the root directory of this repo simply with

```bash
$ source run/run.sh Snakefile
```

or manually with

```bash
$ conda activate snakemake
$ snakemake --use-conda --cores {NUM} -kp
```
