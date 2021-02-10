#!/usr/bin/env nextflow

// Process def
process subset_fasta_by_BLAST_results {
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

    container 'ubuntu:16.04'

    input:
        val opts
        tuple val(meta), path(fasta), path(blast_tab)

    output:
        tuple val(meta), path("*.extracted_hits.fasta"), emit: fasta
        path("*.extracted_hits.tsv"), optional:true, emit: tsv

    script:

        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ' ' + ext_args.trim()
        }

        // If there is a BLAST filter applied and filtered_table is true
        if(opts.filter != "" && opts.filtered_table){
            subset_command = "python3 --fasta ${fasta} --blastoutput ${blast_tab} --filter ${opts.filter} --output ${fasta.simpleName}.extracted_hits.fasta --filtered_table ${fasta.simpleName}.extracted_hits.tsv"
        // If a BLAST filter is applied, but the filtered table is not desired
        } else if (opts.filter != "") {
            subset_command = "python3 --fasta ${fasta} --blastoutput ${blast_tab} --filter ${opts.filter} --output ${fasta.simpleName}.extracted_hits.fasta"
        } else {
            subset_command = "python3 --fasta ${fasta} --blastoutput ${blast_tab} --output ${fasta.simpleName}.extracted_hits.fasta"
        }

        if (params.verbose){
            println ("[MODULE] subset command: " + subset_command)
        }

        // SHELL
        """
        ${subset_command}
        """
}
