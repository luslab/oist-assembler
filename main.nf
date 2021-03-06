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
def busco_genome1_opts = params.modules.busco_genome.clone()
busco_genome1_opts.publish_dir = "busco1"

//params.modules["guppy_basecaller"].flowcell = "FLO-MIN106"
//params.modules["guppy_basecaller"].kit = "SQK-RAD002"
//params.modules["flye"].genome_size = "0.05m"
params.modules.last_db.args = "-uNEAR -R01"

params.modules.tantan_to_GFF3.publish_dir = "tantan"

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { check_max; build_debug_param_summary; luslab_header; check_params } from "./modules/luslab-nf-modules/tools/luslab_util/main.nf" /** required **/
include { fastq_metadata } from "./modules/luslab-nf-modules/tools/metadata/main.nf"
include { filtlong } from "./modules/luslab-nf-modules/tools/filtlong/main.nf"

include { minionqc } from "./modules/luslab-nf-modules/tools/minionqc/main.nf"
include { pairwise_genome_alignment as align_to_self ;
          pairwise_genome_alignment as align_to_reference } from "./modules/local/submodule/pairwise_genome_alignment/main.nf"
include { porechop } from "./modules/luslab-nf-modules/tools/porechop/main.nf"
include { tantan ;
          tantan_to_GFF3 } from "./modules/luslab-nf-modules/tools/tantan/main.nf"
include { flye } from "./modules/luslab-nf-modules/tools/flye/main.nf"
include { racon } from "./modules/luslab-nf-modules/tools/racon/main.nf"
include { purge_haplotigs } from "./modules/luslab-nf-modules/tools/purge_haplotigs/main.nf"
include { repeatmodeler_database } from "./modules/luslab-nf-modules/tools/repeatmodeler/main.nf"
include { repeatmodeler_model } from "./modules/luslab-nf-modules/tools/repeatmodeler/main.nf"
include { repeatmasker } from "./modules/luslab-nf-modules/tools/repeatmodeler/main.nf"
include { minimap2_index } from "./modules/luslab-nf-modules/tools/minimap2/main.nf"
include { minimap2_paf} from "./modules/luslab-nf-modules/tools/minimap2/main.nf"
include { minimap2_sam } from "./modules/luslab-nf-modules/tools/minimap2/main.nf"
include { busco_genome as busco_genome0 } from "./modules/luslab-nf-modules/tools/busco/main.nf"
include { busco_genome as busco_genome1 } from "./modules/luslab-nf-modules/tools/busco/main.nf"
include { augustus_run_custom } from "./modules/luslab-nf-modules/tools/augustus/main.nf"
include { infernal_cmscan } from "./modules/luslab-nf-modules/tools/infernal/main.nf"
include { blast_makeblastdb } from "./modules/luslab-nf-modules/tools/blast/main.nf"
include { blast_blastn } from "./modules/luslab-nf-modules/tools/blast/main.nf"
include { blast_asn_to_tab } from "./modules/luslab-nf-modules/tools/blast/main.nf"
include { mafft } from "./modules/luslab-nf-modules/tools/mafft/main.nf"
include { emboss_seqret } from "./modules/luslab-nf-modules/tools/emboss/main.nf"
include { phyml } from "./modules/luslab-nf-modules/tools/phyml/main.nf"

include { last_db } from "./modules/luslab-nf-modules/tools/last/main.nf"

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

    // Do the assembly
    flye(params.modules["flye"], flye_input_reads)
    // Remap the reads on the assembly
    minimap2_paf(params.modules["minimap2_paf"], flye.out.fasta, fastq_metadata.out)
    // Assess the assembly
    busco_genome0(busco_genome0_opts, flye.out.fasta)
    // Index the assembly
    last_db(params.modules["last_db"], flye.out.fasta)
    blast_makeblastdb(params.modules["blast_makeblastdb"], flye.out.fasta)
    // Analyse tandem repeats in the assembly
    tantan(params.modules["tantan"], flye.out.fasta)
    tantan_to_GFF3(params.modules["tantan_to_GFF3"], tantan.out)

    // Align the assembly to itself and optionally to a reference
    align_to_self(last_db.out, flye.out.fasta)
    if (params.with_reference) {
    channel
        .from(params.with_reference)
        .map { filename -> file(filename, checkIfExists: true) }
        .map { row -> [[sample_id:"ref_genome"], row] }
        .set {ch_reference}
    align_to_reference(last_db.out, ch_reference)
    }

    // Polish with Racon and assess with BUSCO
    justMinimapPaf = minimap2_paf.out.paf.map { row -> row[1] }
    justFlyeAssembly = flye.out.fasta.map { row -> row[1] }
    racon(params.modules["racon"], fastq_metadata.out, justMinimapPaf, justFlyeAssembly)
    busco_genome1(busco_genome1_opts, flye.out.fasta)
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
