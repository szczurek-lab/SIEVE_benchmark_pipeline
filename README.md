# SIEVE_benchmark_pipeline

This repository harbours the benchmarking framework, built upon [Snakemake](https://snakemake.readthedocs.io/en/stable/), for [SIEVE](https://github.com/szczurek-lab/SIEVE). 

## Installation

### Docker (Recommended)

We have configured a docker file containing everything needed to run the benchmarking pipeline. This is recommended for running on server or HPC.

To acqire the docker image, pull from Docker Hub with

```bash
$ docker pull senbaikang/sieve_benchmark:0.2
```

or build from Dockerfile in the root of this repository with

```bash
$ docker build -t sieve_benchmark:0.2 .
```

### Manually (Alternative)

#### Conda environment

This framework contains scripts written in Python 3 and R. We have provided a conda environment file (`environment.yml`), where a list of conda channels and python packages with specific versions are specified. By default, this environment is named `snake`. To create it, please follow the commands below:

```bash
$ git clone https://github.com/szczurek-lab/SIEVE_benchmark_pipeline.git
$ cd SIEVE_benchmark_pipeline
$ mamba env create -f environment.yml
$ mamba activate snake
```

#### Required packages

The following packages should additionally be installed before running the pipeline:

- For R:
  - base >= 4.0
  - stringr
  - scales
  - dplyr
  - optparse
  - ape
  - phangorn

## Configuration

A few files should be configured before running.

1. The pipeline is mainly configured in `config.yaml`. There are some entries to be filled by users, which are marked by a phrase `# TO BE SET` , followed either by `[MANDATORY]` (must be set) or by `[OPTIONAL]` (could be ignored). Users can search for the phrase to set everything up efficiently.

   - The template configuration files for SIEVE are under `templates/`.

   - For the key `[configFiles][SIEVE_simulator]`, a configuration file for the data simulator [SIEVE_simulator](https://github.com/szczurek-lab/SIEVE_simulator) should be specified. Those simulated scenarios used in the SIEVE paper are listed in `simulation_configs/`. For details, check the paper.

2. In `run/run.sh`, users can set a few things, e.g., the name of the conda environment containing snakemake (by default, `snake`), the number of cores to use and their ranges, etc. If you plan to use docker, please skip this step.

### Advanced configuration (Yet unsuporrted for docker)

Since SiFit requires a large amount of memory even working on a small dataset, the framework supports running SiFit alone on another server (referred to as `the remote server`) with the help of git. For this to work, a few things should be noted and configured:

- The machine you plan to run the pipeline (referred to as `the local machine`) and `the remote server` should meet one of the following conditions:
  - They are in the same local network. 
  - If they are not in the same local network, `the local machine` must have a public IP address for `the remote server` to clone the git repository. However, `the remote server` can behind a gate server with a public IP address, specified by the key `[servers][jumpServerOfRemote]` in `config.yaml`.
- In `Snakefile`, comment out `include: "sifit_local.snake"`, and uncomment `# include: "sifit_remote.snake"`. 
- Set up `run/run_remote_sifit_true_monovar.sh` similarly to `run/run.sh`.

## Run the pipeline

### Benchmarking SIEVE by default

#### With docker

The docker image only contains executables of all the benchmarked tools. To run the pipeline, you need to mount the local directory to this repository containing the snakemake rules and supporting scripts to the docker container under `/root/data`:

```bash
$ docker run --name sieve_benchmark -v /local/path/to/SIEVE_benchmark_pipeline:/root/data senbaikang/sieve_benchmark:0.2
```

The console output of snakemake will appear in the terminal. To run the pipeline in the background, add `-d` to the command above before the image name, and access the console outputs through:
  
```bash
$ docker logs sieve_benchmark
```

#### Manually

With everything set up, users can run the pipeline by default with rules defined in `Snakefile` under the root of this repository simply with

```bash
$ source run/run.sh
```

or with

```bash
$ conda activate snake
$ snakemake --use-conda --cores {NUM} -kp
```

### Benchmarking of the efficiency

#### With docker

If benchmarking of the efficiency is of the concern, the snakemake file containing the corresponding rules should be used. Hence, the default commands specified in the docker image must be overwritten. To do so, run the docker container with the following command:

```bash
$ docker run --name sieve_benchmark -v /local/path/to/SIEVE_benchmark_pipeline:/root/data senbaikang/sieve_benchmark:0.2 snakemake --use-conda --cores all -s efficiency_benchmark.snake --rerun-triggers mtime -kp
```

#### Manually

A list of rules defined in `efficiency_benchmark.snake` is readily available and can be run with

```shell
$ source run/run.sh efficiency_benchmark.snake
```

or manually with

```shell
$ conda activate snake
$ snakemake -s efficiency_benchmark.snake --use-conda --cores {NUM} -kp
```

