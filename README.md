# oist-assembler
**oist-assembler** is pipeline that assembles genomes from long-read Oxford Nanopore sequences. It uses the data flow language [Nextflow](https://www.nextflow.io) to coordinate the steps in the pipeline, making it portable and highly configurable. The pipeline involves:
- Basecalling with Guppy
- Read adapter trimming and quality control
- _De novo_ or reference-guided assembly with different assemblers and parameters
- Automated assessment of the assembly

## Installation

    nextflow pull luslab/oist-assembler

To run this pipeline, you will need Nextflow and Docker (or Singularity). If you are planning to use GPU-accelerated basecalling, you'll also need a CUDA-capable GPU and the NVIDIA Container Toolkit.

## Test

In this repository clone:

    nextflow run main.nf -profile <oist,crick,docker,singularity> --input test_workflows/assembly_1.metadata.csv

Note that you can run the pipeline on a local machine with either docker or singularity, but currently, the computational requirements for running the test suite are relatively high.

## Use

Directly from GitHub:

    nextflow run luslab/oist-assembler -profile <oist,crick,docker,singularity> --input [metadata.csv]

In this repository clone:

    nextflow run main.nf -config [pipeline.config] --input [metadata.csv]

## Something Something
:tada:
