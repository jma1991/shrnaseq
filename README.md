# shrnaseq

Snakemake workflow of shRNA-seq and CRISPR-Cas9 genetic screen analysis using edgeR

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.5-brightgreen.svg)](https://snakemake.github.io)
[![GitHub-Actions](https://github.com/zifornd/shrnaseq/workflows/Tests/badge.svg?branch=master)](https://github.com/zifornd/shrnaseq/actions)

## Table of Contents

* [Overview](#Overview)
* [Installation](#Installation)
* [Usage](#Usage)
* [Shiny application](#shiny-application)
* [Acknowledgements](#Acknowledgements)

## Overview

This workflow is used to analysis shRNA-seq and CRISPR/cas9 genetic screens. It uses edgeR to perform data processing, quality control, differential expression analysis, and gene set testing, with batch correction implemented. 

## Installation

Install snakemake using the mamba package manager:

```bash
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake

$ mamba activate snakemake
```
    
Pull the workflow to your project directory:

```bash
$ git pull https://github.com/zifornd/shrnaseq
```
    
## Usage

Configure the workflow by editing the `config.yaml` file:

```console
$ nano config/config.yaml
```

Run the workflow: 

```console
$ snakemake --use-conda --cores all
```
For further details on Snakemake, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

## Shiny application

The output of this workflow can be visualised in Shiny from RStudio, see the [shRNAseq shiny documentation](https://github.com/zifornd/shrnaseq-shiny).

## Documentation

See the [Documentation](workflow/DOCUMENTATION.md) file for configuration and output information.

## Acknowledgements

This workflow is based on the following research article:

```
Dai Z, Sheridan JM, Gearing LJ et al. edgeR: a versatile tool for the analysis of shRNA-seq and CRISPR-Cas9 genetic screens [version 2; peer review: 3 approved]. F1000Research 2014, 3:95 (https://doi.org/10.12688/f1000research.3928.2)
```