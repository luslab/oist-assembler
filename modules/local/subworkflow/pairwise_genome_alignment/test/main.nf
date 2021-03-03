#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { last_db } from "../../../luslab-modules/tools/last/main.nf"
include { pairwise_genome_alignment } from '../main.nf'

genome1 = [
    [[sample_id:"genome1"], "../../test_workflows/assembly_1_data/reference_genome_library/GCF_000146045.2_R64_genomic.fa"],
]

genome2 = [
    [[sample_id:"genome2"], "../../test_workflows/assembly_1_data/reference_genome_library/GCF_000146045.2_R64_genomic.fa"],
]

channel
    .from(genome1)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_genome1}

channel
    .from(genome2)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_genome2}


workflow {
    last_db(params.modules.last_db, ch_genome1)
    pairwise_genome_alignment(
        last_db.out,
        ch_genome2)
}
