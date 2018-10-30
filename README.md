[![singularity](https://img.shields.io/badge/singularity-%3E%3D%202.4.2-blue.svg)](http://singularity.lbl.gov/)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.2-brightgreen.svg)](https://www.nextflow.io/)

**heinzlab/smrna-seq-pipeline** is a bioinformatics best-practice analysis pipeline used for small RNA sequencing data at the Heinz Lab at UCSD. The pipeline is adapted from the nf-core set of Nextflow pipelines.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

This pipeline can be run on any Unix based environment that has Nextflow installed.

## What is here?
* **README.md:** is what you are reading, which has a complete walkthrough of building and running the container.
* **Singularity:** includes the build recipe for the main Singularity container. Can more or less be copied over for implementation of other genomics pipelines.
* **nextflow.config:** the main nextflow config file with profile and default parameter information.

## Installation
### NextFlow installation
See https://github.com/SciLifeLab/NGI-NextflowDocs for instructions on how to install and configure Nextflow.

### Pipeline installation
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub when running `heinzlab/smrna-seq-pipeline` is specified as the pipeline name.
