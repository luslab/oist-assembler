#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline run params
-------------------------------------------------------------------------------------------------------------------------------*/
params.modules['guppy_basecaller'].flowcell = "FLO-MIN106"
params.modules['guppy_basecaller'].kit = "SQK-RAD002"
params.modules['flye'].genome_size = "0.05m"

/*-----------------------------------------------------------------------------------------------------------------------------
Module inclusions
-------------------------------------------------------------------------------------------------------------------------------*/

include { check_max; build_debug_param_summary; luslab_header; check_params } from './luslab-modules/tools/luslab_util/main.nf' /** required **/
include { guppy_basecaller } from './luslab-modules/tools/guppy/main.nf'
include { fastq_metadata } from './luslab-modules/tools/metadata/main.nf'
include { minionqc } from './luslab-modules/tools/minionqc/main.nf'
include { nanoplot} from './luslab-modules/tools/nanoplot/main.nf'
include { flye } from './luslab-modules/tools/flye/main.nf'

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
	// Collect metadata
	fastq_metadata(params.input)
	// fastq_metadata.out.metadata | view
	// Call bases with guppy
	guppy_basecaller(params.modules['guppy_basecaller'], fastq_metadata.out.metadata)
	// Do some statistics on the basecalled data
	//minionqc(fastq_metadata.out.metadata, guppy_basecaller.out.report)
	//nanoplot(params.modules['nanoplot'], guppy_basecaller.out.fastq)
	// Do the assembly
	flye(params.modules['flye'], guppy_basecaller.out.fastq)
}

workflow.onComplete {
}

/*------------------------------------------------------------------------------------*/
