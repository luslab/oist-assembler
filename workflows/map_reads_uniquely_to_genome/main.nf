#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { last_train ;
          last_align ;
          last_filter_one_to_many as last_filter_maf ;
          last_convert_maf } from "../../local/modules/last/main.nf"

def last_train_opts                        = params.modules.last_train.clone()
def last_align_opts                        = params.modules.last_align.clone()
def last_filter_maf_opts                   = params.modules.last_filter_one_to_many.clone()
def last_convert_maf_opts                  = params.modules.last_convert_maf.clone()

last_train_opts.args                       += "-Q0 --revsym"
last_align_opts.args                       += "-Q0 -E0.05 -C2"
last_convert_maf_opts.args                 += "--noheader"
last_filter_maf_opts.args                  += "-m1e-20"

last_convert_maf_opts.suffix               = "tab"

last_train_opts.publish_dir                = "read_align/align"
last_align_opts.publish_dir                = "read_align/align"
last_filter_maf_opts.publish_dir           = "read_align/align"
last_convert_maf_opts.publish_dir          = "read_align/align"

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
        last_convert_maf(
            last_convert_maf_opts,
            last_filter_maf.out.maf)
    emit:
        maf = last_filter_maf.out.maf
        tab = last_convert_maf.out
}
