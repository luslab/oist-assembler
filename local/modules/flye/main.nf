#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Process definition
process flye {
    label "high_cores"
    label "max_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "https://www.dropbox.com/s/mpzgq5mqjtf31gm/Flye-flye.2.8.3-b1763.sif?dl=1"

    input:
        val opts
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.sample_id}/assembly.fasta"), emit: fasta
	tuple val(meta), path("${meta.sample_id}/assembly_graph.gfa"), emit: gfa
        tuple path("${meta.sample_id}/assembly_info.txt"), path("${meta.sample_id}/flye.log"), emit: log

    script:

    args = ""

    if(opts.args && opts.args != "") {
        ext_args = opts.args
        args += " " + ext_args.trim()
    }

    //Build the command line options
    flye_command = "flye $args --genome-size ${opts.genome_size} \
			--threads ${task.cpus} \
			--out-dir ${meta.sample_id} \
			--nano-corr $reads"

	//SHELL
    """
    ${flye_command}
    """
}
