# references

Creation of references, PoN, and background

![CI](https://github.com/hydra-genetics/references/actions/workflows/ci.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

## :speech_balloon: Introduction

## :heavy_exclamation_mark: Dependencies

To run this workflow, the following tools need to be available:

![python](https://img.shields.io/badge/python-3.12-blue)
[![snakemake](https://img.shields.io/badge/snakemake-9.0.0-blue)](https://snakemake.readthedocs.io/en/stable/)
[![apptainer](https://img.shields.io/badge/apptainer-1.4.5-blue)](https://apptainer.org/)

## :school_satchel: Preparations

### Sample data

1. Add all sample ids to `samples.tsv` in the column `sample`.
2. Add all sample data information to `units.tsv`. Each row represents a `fastq` file pair with
corresponding forward and reverse reads. Also indicate the sample id, run id and lane number, adapter.

### Reference data

1. You need a ...

## :white_check_mark: Testing

The workflow repository contains a small test dataset `.tests/integration` which can be run like so:

```bash
cd .tests/integration
snakemake -s ../../workflow/Snakefile -j1 --configfile config.yaml --software-deployment-method apptainer
```

## :rocket: Usage

The workflow is designed for WGS data meaning huge datasets which require a lot of compute power. For
HPC clusters, it is recommended to use a cluster profile and run something like:

```bash
snakemake -s /path/to/workflow/Snakefile --profile my-awesome-profile
```

## :judge: Rule Graph

![rule_graph](https://raw.githubusercontent.com/path.../rulegraph.svg)
