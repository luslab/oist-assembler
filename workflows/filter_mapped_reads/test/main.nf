#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { filter_mapped_reads }          from '../main.nf'
include { last_db }                      from '../../../luslab-modules/tools/last/main.nf'
include { map_reads_uniquely_to_genome } from '../../map_reads_uniquely_to_genome/main.nf'

dups = [
    [[sample_id:"genome"], "../../test_workflows/assembly_1_data/test_data/purge_dups.bed"],
]

genome = [
    [[sample_id:"genome"], "../../test_workflows/assembly_1_data/reference_genome_library/assembly.fasta"],
]

reads = [
    [[sample_id:"genome"], "../../test_workflows/assembly_1_data/test_data/nanopore_ERR1883389.0.85.fastq"],
]

channel
    .from(dups)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_dups}

channel
    .from(genome)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_genome}

channel
    .from(reads)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_reads}

def opts = params.modules.last_filter_maf.clone()
opts.publish_dir = "filtered_reads"

workflow {
    last_db(params.modules.last_db, ch_genome)
    map_reads_uniquely_to_genome(last_db.out, ch_reads)
    filter_mapped_reads(opts, ch_dups, map_reads_uniquely_to_genome.out.maf, ch_reads)
}
