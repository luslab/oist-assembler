#!/usr/bin/env nextflow
/*
=================================================================================
                                OIST-assembler
=================================================================================
 oist-assembler genome assembly pipeline
 #### Homepage / Documentation
 https://github.com/luslab/oist-assembler
---------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*-------------------------------------------------------------------------------
 * Print help message
 *-------------------------------------------------------------------------------
 */

include { check_max; build_debug_param_summary; luslab_header; check_params } from "./modules/luslab-nf-modules/tools/luslab_util/main.nf" /** required **/

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run luslab/oist-assembler --input samplesheet.csv -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

/*-------------------------------------------------------------------------------
 * Print parameter summary
 *-------------------------------------------------------------------------------
 */
log.info luslab_header()
if(params.verbose) {
    log.info build_debug_param_summary()
}

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

/*-------------------------------------------------------------------------------
 * Parameter checks
 *-------------------------------------------------------------------------------
 */
// Check inputs
check_params(["input"])

/*-------------------------------------------------------------------------------
 * Main workflow
 *-------------------------------------------------------------------------------
 */


workflow {
    /*
     *  SUBWORKFLOW: Run main oist-assembly analysis pipeline
     */
    include { ASSEMBLE } from './assemble.nf'
    ASSEMBLE ()
}



/*-------------------------------------------------------------------------------*/
