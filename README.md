# oist-assembler
This pipeline assembles genomes starting from long-read Oxford Nanopore sequences. It uses the data flow language _nextflow_ to coordinate the steps in the pipeline, which involves:
- Read-level quality control
- Basecalling with Guppy
- Read adapter trimming

- _De novo_ or reference-guided assembly with different assemblers and parameters
- Automated assessment of the assembly

## Installation

    nextflow pull luslab/oist-assembler

To run this pipeline, you will need Nextflow and Docker (or Singularity). If you are planning to use GPU-accelerated basecalling, you'll also need a CUDA-capable GPU and the NVIDIA Container Toolkit.

## Test

In this repository clone:

    nextflow run main.nf -profile <oist,crick> --input test_workflows/assembly_1.metadata.csv

## Use

Directly from GitHub:

    nextflow run luslab/oist-assembler -profile <oist,crick> --input [metadata.csv]

In this repository clone:

    nextflow run main.nf -config [pipeline.config] --input [metadata.csv]

## Something Something
:tada:
