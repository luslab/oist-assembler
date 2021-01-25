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

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { check_max; build_debug_param_summary; luslab_header; check_params } from "../luslab-modules/tools/luslab_util/main.nf" /** required **/
include { fastq_metadata } from "../luslab-modules/tools/metadata/main.nf"

include { minionqc } from "../luslab-modules/tools/minionqc/main.nf"
include { porechop } from "../luslab-modules/tools/porechop/main.nf"
include { flye } from "../luslab-modules/tools/flye/main.nf"
include { racon } from "../luslab-modules/tools/racon/main.nf"
include { purge_haplotigs } from "../luslab-modules/tools/purge_haplotigs/main.nf"
include { repeatmodeler_database } from "../luslab-modules/tools/repeatmodeler/main.nf"
include { repeatmodeler_model } from "../luslab-modules/tools/repeatmodeler/main.nf"
include { repeatmasker } from "../luslab-modules/tools/repeatmodeler/main.nf"
include { minimap2_index } from "../luslab-modules/tools/minimap2/main.nf"
include { minimap2_paf} from "../luslab-modules/tools/minimap2/main.nf"
include { minimap2_sam } from "../luslab-modules/tools/minimap2/main.nf"
include { busco_genome as busco_genome0 } from "../luslab-modules/tools/busco/main.nf"
include { busco_genome as busco_genome1 } from "../luslab-modules/tools/busco/main.nf"
include { augustus_run_custom } from "../luslab-modules/tools/augustus/main.nf"
include { infernal_cmscan } from "../luslab-modules/tools/infernal/main.nf"
include { blast_makeblastdb } from "../luslab-modules/tools/blast/main.nf"
include { blast_blastn } from "../luslab-modules/tools/blast/main.nf"
include { blast_asn_to_tab } from "../luslab-modules/tools/blast/main.nf"
include { mafft } from "../luslab-modules/tools/mafft/main.nf"
include { emboss_seqret } from "../luslab-modules/tools/emboss/main.nf"
include { phyml } from "../luslab-modules/tools/phyml/main.nf"

// The following lines are used to import LAST commands with reasonable presets for
// doing the tasks laid out in the process names. The arguments of these presets are in:
//     luslab-modules/tools/last/test/last.config
include { last_db } from "../luslab-modules/tools/last/main.nf"
include { last_train as last_train_genome_to_genome_distant } from "../luslab-modules/tools/last/main.nf"
include { last_train as last_train_reads_to_genome } from "../luslab-modules/tools/last/main.nf"
include { last_align as last_align_genome_to_genome_distant } from "../luslab-modules/tools/last/main.nf"
include { last_align as last_align_reads_to_genome } from "../luslab-modules/tools/last/main.nf"
include { last_convert_maf_to_sam } from "../luslab-modules/tools/last/main.nf"
include { last_filter_maf as last_filter_maf_distant } from "../luslab-modules/tools/last/main.nf"
include { last_dotplot as last_dotplot_distant } from "../luslab-modules/tools/last/main.nf"

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
    // Enumerate the genomes in the reference_genome_library directory
    ref_genome_library = fastq_metadata.out
        .map { row -> row[0].ref_genome_library }
        .map { dir -> file(dir).listFiles() }
    ref_marker_gene_library = fastq_metadata.out
        .map { row -> row[0].ref_marker_gene_library }

    // Do some statistics on the input reads
    //minionqc(params.modules["minionqc"], fastq_metadata.out.metadata)

    //minionqc(params.modules["minionqc"], guppy_basecaller.out.sequencing_summary)
    // Do the assembly
    flye(params.modules["flye"], fastq_metadata.out)
    // Remap the reads on the assembly
    minimap2_paf(params.modules["minimap2_paf"], flye.out.fasta, fastq_metadata.out)
    // Assess the assembly
    busco_genome0(busco_genome0_opts, flye.out.fasta)
    // Index the assembly
    last_db(params.modules["last_db"], flye.out.fasta)
    blast_makeblastdb(params.modules["blast_makeblastdb"], flye.out.fasta)
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
