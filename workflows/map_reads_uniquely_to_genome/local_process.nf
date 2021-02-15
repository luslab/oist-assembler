#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process last_filter_one_to_many {
    label "min_cores"
    label "avg_mem"
    label "regular_queue"

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container "quay.io/biocontainers/last:1060--h8b12597_0"

    input:
        val opts
        tuple val(meta), path(unfiltered_maf)

    output:
        tuple val(meta), path("*.filtered.maf"), emit: maf

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        prefix = opts.suffix ? "${meta.sample_id}${opts.suffix}" : "${meta.sample_id}"

        last_rerarrange_command = "last-split ${unfiltered_maf} > ${unfiltered_maf.simpleName}.filtered.maf"

        if (params.verbose){
            println ("[MODULE] LAST one-to-one rearrangement command: " + last_rerarrange_command)
        }

        //SHELL
        """
        ${last_rerarrange_command}
        """
}
