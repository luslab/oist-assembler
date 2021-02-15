#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { last_db } from "../../../luslab-modules/tools/last/main.nf"
include { map_reads_uniquely_to_genome } from '../main.nf'

genome = [
    [[sample_id:"genome"], "../../test_workflows/assembly_1_data/reference_genome_library/GCF_000146045.2_R64_genomic.fa"],
]

reads = [
    [[sample_id:"genome"], "../../test_workflows/assembly_1_data/test_data/nanopore_ERR1883389.0.85.fastq.gz"],
]

channel
    .from(genome)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_genome}

channel
    .from(reads)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_reads}


workflow {
    last_db(params.modules.last_db, ch_genome)
    map_reads_uniquely_to_genome(
        last_db.out,
        ch_reads)
}
