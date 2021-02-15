#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { last_train ;
          last_align } from "../../luslab-modules/tools/last/main.nf"
include { last_filter_one_to_many as last_filter_maf } from "./local_process.nf"

def last_train_opts                        = params.modules.last_train.clone()
def last_align_opts                        = params.modules.last_align.clone()
def last_filter_maf_opts                   = params.modules.last_filter_maf.clone()

last_train_opts.args                       += "-Q0"
last_align_opts.args                       += "-Q0"

last_train_opts.publish_dir                = "read_align/align"
last_align_opts.publish_dir                = "read_align/align"
last_filter_maf_opts.publish_dir           = "read_align/align"

workflow map_reads_uniquely_to_genome {
    take:
        index
        reads
    main:
        last_train(
            last_train_opts,
            index,
            reads)
        last_align(
            last_align_opts,
            index,
            last_train.out.par,
            reads)
        last_filter_maf(
            last_filter_maf_opts,
            last_align.out.maf)
    emit:
        last_filter_maf.out.maf
}