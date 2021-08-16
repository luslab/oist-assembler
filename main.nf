#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline run params
-------------------------------------------------------------------------------------------------------------------------------*/


params.modules.blast_makeblastdb.dbtype = "nucl"

params.modules.busco_genome.args = "--lineage_dataset metazoa_odb10"
def busco_genome0_opts = params.modules.busco_genome.clone()
busco_genome0_opts.publish_dir = "busco0"
def busco_genome0p_opts = params.modules.busco_genome.clone()
busco_genome0p_opts.publish_dir = "busco0p"
def busco_genome1_opts = params.modules.busco_genome.clone()
busco_genome1_opts.publish_dir = "busco1"

def filter_mapped_reads_opts = params.modules.last_filter_one_to_one.clone()
filter_mapped_reads_opts.publish_dir = "filtered_reads"


def minimap2_paf_flye_opts = params.modules.minimap2_paf.clone()
minimap2_paf_flye_opts.publish_dir = "minimap2_paf_flye"
def minimap2_paf_purged_opts = params.modules.minimap2_paf.clone()
minimap2_paf_purged_opts.publish_dir = "minimap2_paf_purged"

//params.modules["guppy_basecaller"].flowcell = "FLO-MIN106"
//params.modules["guppy_basecaller"].kit = "SQK-RAD002"
//params.modules["flye"].genome_size = "0.05m"
params.modules.last_db.args = "-uRY32 -R01"

params.modules.tantan_to_GFF3.publish_dir = "tantan"

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { check_max; build_debug_param_summary; luslab_header; check_params } from "./luslab-modules/tools/luslab_util/main.nf" /** required **/
include { fastq_metadata } from "./luslab-modules/tools/metadata/main.nf"
include { filter_mapped_reads } from "./workflows/filter_mapped_reads/main.nf"
include { filtlong } from "./luslab-modules/tools/filtlong/main.nf"
include { map_reads_uniquely_to_genome } from "./workflows/map_reads_uniquely_to_genome/main.nf"
include { minionqc } from "./luslab-modules/tools/minionqc/main.nf"
include { pairwise_genome_alignment as align_to_self ;
          pairwise_genome_alignment as align_to_reference } from "./workflows/pairwise_genome_alignment/main.nf"
include { porechop } from "./luslab-modules/tools/porechop/main.nf"
include { tantan ;
          tantan_to_GFF3 } from "./luslab-modules/tools/tantan/main.nf"
include { flye } from "./local/modules/flye/main.nf"
include { racon } from "./luslab-modules/tools/racon/main.nf"
include { purge_haplotigs } from "./luslab-modules/tools/purge_haplotigs/main.nf"
include { purge_dups } from "./luslab-modules/tools/purge_dups/main.nf"
include { repeatmodeler_database } from "./luslab-modules/tools/repeatmodeler/main.nf"
include { repeatmodeler_model } from "./luslab-modules/tools/repeatmodeler/main.nf"
include { repeatmasker } from "./luslab-modules/tools/repeatmodeler/main.nf"
include { minimap2_index                      ;
          minimap2_paf as minimap2_paf_flye   ;
          minimap2_paf as minimap2_paf_purged ;
          minimap2_sam                        } from "./luslab-modules/tools/minimap2/main.nf"
include { busco_genome as busco_genome0  ;
          busco_genome as busco_genome0p ;
          busco_genome as busco_genome1  } from "./luslab-modules/tools/busco/main.nf"
include { augustus_run_custom } from "./luslab-modules/tools/augustus/main.nf"
include { infernal_cmscan } from "./luslab-modules/tools/infernal/main.nf"
include { blast_makeblastdb } from "./luslab-modules/tools/blast/main.nf"
include { blast_blastn } from "./luslab-modules/tools/blast/main.nf"
include { blast_asn_to_tab } from "./luslab-modules/tools/blast/main.nf"
include { mafft } from "./luslab-modules/tools/mafft/main.nf"
include { emboss_seqret } from "./luslab-modules/tools/emboss/main.nf"
include { phyml } from "./luslab-modules/tools/phyml/main.nf"

include { last_db } from "./luslab-modules/tools/last/main.nf"

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/

// Show banner and param summary
log.info luslab_header()
if(params.verbose)
    log.info build_debug_param_summary()

// Check inputs
check_params(["input"])

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    // Collect metadata.
    // The metadata reads in a .csv file, which is structured to contain a "data1" column corresponding
    // to the nanopore reads, and then any number of additional labeled features beyond the first column.
    fastq_metadata(params.input)

    // Separate out the various data from the metadata header
    ch_illumina = fastq_metadata.out
        .map { row -> [ [sample_id:row[0].sample_id], row[0].illumina1, row[0].illumina2] }
    genome_size = fastq_metadata.out
        .map { row -> row[0].genome_size}
    augustus_model = fastq_metadata.out
        .map { row -> row[0].ref_augustus_model }
    repeat_library = fastq_metadata.out
        .map { row -> row[0].ref_repeat_library }
    rna_library = fastq_metadata.out
        .map { row -> row[0].ref_rna_library }
    ref_marker_gene_library = fastq_metadata.out
        .map { row -> row[0].ref_marker_gene_library }

    // Do some statistics on the input reads
    //minionqc(params.modules["minionqc"], fastq_metadata.out.metadata)

    //minionqc(params.modules["minionqc"], guppy_basecaller.out.sequencing_summary)

    // Optionally filter the reads with filtlong
    if (params.with_filtlong == true) {
        params.modules.filtlong.args = "--target_bases ${params.filtlong_target}"
        filtlong(params.modules["filtlong"], fastq_metadata.out)
        flye_input_reads = filtlong.out.fastq
    } else {
        flye_input_reads = fastq_metadata.out
    }

    flye(params.modules["flye"], flye_input_reads)

    // Index the assembly
    last_db(params.modules["last_db"], flye.out.fasta)
    blast_makeblastdb(params.modules["blast_makeblastdb"], flye.out.fasta)

    // Remap the reads on the assembly
    minimap2_paf_flye(minimap2_paf_flye_opts, flye.out.fasta, fastq_metadata.out)
    if ( params.with_remap_reads == true | params.with_purge_reads == true ) {
        map_reads_uniquely_to_genome(last_db.out, fastq_metadata.out)
    }
    // Assess the assembly
    busco_genome0(busco_genome0_opts, flye.out.fasta)

    // Optionally purge duplicates and assess the result
    if (params.with_purge_dups == true | params.with_purge_reads == true) {
        purge_dups(params.modules["purge_dups"], flye.out.fasta, minimap2_paf_flye.out.paf)
        genome_assembly = purge_dups.out.purged_fasta
        busco_genome0p(busco_genome0p_opts, genome_assembly)
        minimap2_paf_purged(minimap2_paf_purged_opts, genome_assembly, fastq_metadata.out)
        minimapped_reads = minimap2_paf_purged
    } else {
        genome_assembly = flye.out.fasta
        minimapped_reads = minimap2_paf_flye
    }

    if (params.with_purge_reads == true) {
        filter_mapped_reads(filter_mapped_reads_opts,
                            purge_dups.out.bed,
                            map_reads_uniquely_to_genome.out.tab,
                            fastq_metadata.out)
    }

    // Analyse tandem repeats in the assembly
    tantan(params.modules["tantan"], genome_assembly)
    tantan_to_GFF3(params.modules["tantan_to_GFF3"], tantan.out)

    // Align the assembly to itself and optionally to a reference
    align_to_self(last_db.out, genome_assembly)
    if (params.with_reference) {
    channel
        .from(params.with_reference)
        .map { filename -> file(filename, checkIfExists: true) }
        .map { row -> [[sample_id:"ref_genome"], row] }
        .set {ch_reference}
    align_to_reference(last_db.out, ch_reference)
    }

    // Polish with Racon and assess with BUSCO
    racon(params.modules["racon"], fastq_metadata.out, minimapped_reads.out.paf, genome_assembly)
    busco_genome1(busco_genome1_opts, genome_assembly)
}

workflow.onError {
    println "Don't panic. The pipeline has failed, yielding the following error: ${workflow.errorMessage}"
}

workflow.onComplete {
    // Print a nice little message :-)
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*------------------------------------------------------------------------------------*/
