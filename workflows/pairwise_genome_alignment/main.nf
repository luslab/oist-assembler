#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { last_train ;
          last_align ;
          last_filter_maf ;
          last_convert_maf ;
          last_dotplot as last_dotplot_many2many ;
          last_dotplot as last_dotplot_one2one   } from "../../luslab-modules/tools/last/main.nf"

def last_train_opts                        = params.modules.last_train.clone()
def last_align_opts                        = params.modules.last_align.clone()
def last_filter_maf_opts                   = params.modules.last_filter_maf.clone()
def last_convert_maf_opts                  = params.modules.last_convert_maf.clone()
def last_dotplot_opts_many2many            = params.modules.last_dotplot.clone()
def last_dotplot_opts_one2one              = params.modules.last_dotplot.clone()

last_train_opts.publish_dir                = "pair_align/align"
last_align_opts.publish_dir                = "pair_align/align"
last_filter_maf_opts.publish_dir           = "pair_align/align"
last_convert_maf_opts.suffix               = "gff"
last_convert_maf_opts.args                 = "-J 2e5"
last_dotplot_opts_many2many.publish_dir    = "pair_align/plotMany2Many"
last_dotplot_opts_one2one.publish_dir      = "pair_align/plotOne2One"

workflow pairwise_genome_alignment {
    take:
        index
        genome
    main:
        last_train(
            last_train_opts,
            index,
            genome)
        last_align(
            last_align_opts,
            index,
            last_train.out.par,
            genome)
        last_filter_maf(
            last_filter_maf_opts,
            last_align.out.maf)
        last_convert_maf(
            last_convert_maf_opts,
            last_filter_maf.out.maf)
        last_dotplot_many2many(
            last_dotplot_opts_many2many,
            last_align.out.tab)
        last_dotplot_one2one(
            last_dotplot_opts_one2one,
            last_filter_maf.out.tab)
    emit:
        last_filter_maf.out.maf
}
