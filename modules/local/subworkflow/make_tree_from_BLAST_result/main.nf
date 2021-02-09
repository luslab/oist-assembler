#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { blast_makeblastdb } from '../../../luslab-nf-modules/tools/blast/main.nf'
include { blast_blastn } from '../../../luslab-nf-modules/tools/blast/main.nf'
include { blast_asn_to_tab } from '../../../luslab-nf-modules/tools/blast/main.nf'
include { subset_fasta_by_BLAST_results } from '../../process/subset_fasta_by_BLAST_results.nf'
include { cdhit_nucl } from '../../../luslab-nf-modules/tools/cdhit/main.nf'
include { }

workflow sort_index_bam {
    take: tuple_meta_fasta
    take:

    main:
        // def modules {
        //     'samtools_sort' {
        //         args             = ""
        //         suffix           = "_sorted.bam"
        //         publish_dir      = "samtools_sort"
        //         publish_results  = "none"
        //     }
        //     'samtools_index' {
        //         args             = ""
        //         suffix           = ".bam.bai"
        //         publish_dir      = "samtools_index"
        //         publish_results  = "none"
        //     }
        // }

        // Define workflow parameters
        params.modules['samtools_sort'].publish_results = 'none'
        params.modules['samtools_index'].publish_results = 'none'
        params.modules['samtools_sort'].publish_dir = '_'
        params.modules['samtools_index'].publish_dir = '_'

        // Sort bam
        samtools_sort( params.modules['samtools_sort'], tuple_meta_bam )

        // Index bam
        samtools_index( params.modules['samtools_index'], samtools_sort.out.bam )

    emit:
        bam_bai = samtools_index.out.bam
}
