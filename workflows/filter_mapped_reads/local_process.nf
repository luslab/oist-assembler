#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process filter_mapped_reads {
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
        tuple val(meta), path(dups)
        tuple val(meta), path(filtered_alignment)
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("reads_kept.fq"), emit: fastq

    script:
        args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

        // List removed scaffolds
        make_grep_exprs  = "grep HAP\$ ${dups} | cut -f 1 | sed 's/\$/\$/' > scaf_rem.txt"
        // List reads matching removed scaffolds
        make_read_list   = "maf-convert sam ${filtered_alignment} | sed 1d | cut -f 1,3 | grep -f scaf_rem.txt | cut -f 1 | uniq > reads_rem.txt"
        // Functions to convert FASTQ to and from one-line intermediary format.
        to_one_line      = "toOneLine() { paste <(sed -n 1~4p \$1) <(sed -n 2~4p \$1) <(sed -n 3~4p \$1) <(sed -n 4~4p \$1) ; }"
        to_four_lines    = "toFourLines() { sed 's/\\t/\\n/g' ; }"
        // Remove reads matching haplotypes scaffolds from original reads.
        filter_HAP_reads = "toOneLine ${reads} | grep -v -f reads_rem.txt | toFourLines > reads_kept.fq"

        if (params.verbose){
            println ("[MODULE] filter_mapped_reads make grep expressions command: " + make_grep_exprs)
            println ("[MODULE] filter_mapped_reads make read list: " + make_read_list)
            println ("[MODULE] filter_mapped_reads shell function FASTQ to one line: " + to_one_line)
            println ("[MODULE] filter_mapped_reads shell function one line to FASTQ: " + to_four_lines)
            println ("[MODULE] filter_mapped_reads filter out reads matching haplotypes: " + filter_HAP_reads)
        }

        //SHELL
        """
        ${make_grep_exprs}
        ${make_read_list}
        ${to_one_line}
        ${to_four_lines}
        ${filter_HAP_reads}
        """
}
