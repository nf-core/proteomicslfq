# nf-core/proteomicslfq

**Proteomics label-free quantification (LFQ) analysis pipeline using OpenMS and MSstats, with feature quantification, feature summarization, quality control and group-based statistical analysis.**.

[![Build Status](https://travis-ci.com/nf-core/proteomicslfq.svg?branch=master)](https://travis-ci.com/nf-core/proteomicslfq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/proteomicslfq.svg)](https://hub.docker.com/r/nfcore/proteomicslfq)

## Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


## Documentation
The nf-core/proteomicslfq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
    * [Reference genomes](docs/configuration/reference_genomes.md)  
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits
nf-core/proteomicslfq was originally written by [Julianus Pfeuffer](https://github.com/jpfeuffer), [Lukas Heumos](github.com/zethson), [Timo Sachsenberg](https://github.com/timosachsenberg) and [Leon Bichmann](https://github.com/Leon-Bichmann).
