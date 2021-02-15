#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { filter_mapped_reads as filter_mapped_reads_process } from "./local_process.nf"

workflow filter_mapped_reads {
    take:
        opts
        dups
        align
        reads
    main:
        filter_mapped_reads_process(
            opts,
            dups,
            align,
            reads)
}
