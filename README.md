# oist-assembler
This pipeline assembles genomes starting from long-read Oxford Nanopore sequences. It uses the data flow language _nextflow_ to coordinate the steps in the pipeline, which involves:
- Automated QC of long-read genomic sequences (and Illumina reads for polishing)
- Assessment of contamination
- _De novo_ or reference-guided assembly with different assemblers and parameters
- Automated assessment of the assembly

## Installation
Docker, singularity, nextflow, etc.

## Use
nextflow run [pipeline.nf] --input ./[metadata.csv]

## Something Something
:tada:
