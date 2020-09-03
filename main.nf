#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * ampatchlab/nf-pdx: Patient-derived xenograft Nextflow pipeline
 *
 * Copyright (C) 2020 QIMR Berghofer Medical Research Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


nextflow.preview.dsl=2

import nextflow.config.ConfigParser

nextflow_config = file( "${baseDir}/nextflow.config" ).text
parsed_config = new ConfigParser().setIgnoreIncludes( true ).parse( nextflow_config )
defaults = parsed_config.params

check_params()


/*
 * Modules
 */

include { parse_input_csv } from './functions/input_csv_parsers.nf' params( params )
include { germline_pdx } from './workflows/germline_pdx.nf' params( params )


/*
 * Params
 */

// Cutadapt adapter files
params.r1_adapter_file = params.adapters in params.adapter_files
    ? params.adapter_files[ params.adapters ].r1_adapters
    : "${baseDir}/resource-adapters/null-1.fa"
params.r2_adapter_file = params.adapters in params.adapter_files
    ? params.adapter_files[ params.adapters ].r2_adapters
    : "${baseDir}/resource-adapters/null-2.fa"

// Reference genome files
params.human_ref_fasta = params.human_genome_assembly in params.genome_assemblies
    ? params.genome_assemblies[ params.human_genome_assembly ].ref_fasta
    : null
params.human_vep_cache = params.human_genome_assembly in params.genome_assemblies
    ? params.genome_assemblies[ params.human_genome_assembly ].vep_cache
    : null
params.mouse_ref_fasta = params.mouse_genome_assembly in params.genome_assemblies
    ? params.genome_assemblies[ params.mouse_genome_assembly ].ref_fasta
    : null
params.mouse_vep_cache = params.mouse_genome_assembly in params.genome_assemblies
    ? params.genome_assemblies[ params.mouse_genome_assembly ].vep_cache
    : null


/*
 * Workflows
 */

workflow {

    input_readgroups = parse_input_csv( params.csv )

    adapters = [ params.r1_adapter_file, params.r2_adapter_file ]

    human_params = [
        params.human_genome_assembly,
        params.human_ref_fasta,
        params.human_vep_cache,
    ]

    mouse_params = [
        params.mouse_genome_assembly,
        params.mouse_ref_fasta,
        params.mouse_vep_cache,
    ]

    workflow_inputs = [
        input_readgroups,
        adapters,
        *human_params,
        *mouse_params,
        params.multiqc_config,
    ]

    if( params.sample_pairs ) {

        false
    }
    else {

        germline_pdx( *workflow_inputs )
    }
}

workflow.onComplete {

    log.info "Workflow completed at: ${workflow.complete}"
    log.info "Time taken: ${workflow.duration}"
    log.info "Execution status: ${workflow.success ? 'success' : 'failed'}"
    log.info "Output directory: ${params.publish_dir}"
}

workflow.onError {

    log.info "Execution halted: ${workflow.errorMessage}"
}


/*
 * Functions
 */

def check_params() {

    if( params.help ) {
        usage()
        exit 0
    }

    if( params.version ) {
        log.info( workflow.manifest.version )
        exit 0
    }
}

def usage() {

    log.info"""
    Usage:
        nextflow run -profile <profile> -revision <revision> ampatchlab/nf-pdx [options]


    Nextflow execution options:

        -profile STR
            Nextflow configuration profile to use. Available profiles include:
            'awsbatch', 'conda', 'docker' and 'singularity'

        -revision STR
            Git branch/tag (version) of this workflow to use

        -work-dir DIR
            Directory where intermediate result files are stored

        -help
            Show additional execution options and exit


    Required params:

        --csv FILE
            Comma-separated list of sample and readgroup inputs


    Adapter trimming params:

        --adapters STR
            The adapters to trim [Either: ${defaults.adapter_files.keySet().join(", ")}; Default: ${defaults.adapters}]

        --r1_adapter_file FILE
            Override the R1 adapter file with FILE [Default: ${defaults.r1_adapter_file ?: null}]

        --r2_adapter_file FILE
            Override the R2 adapter file with FILE [Default: ${defaults.r2_adapter_file ?: null}]


    Output params:

        --publish_dir DIR
            Path where the results will be published [Default: ${defaults.publish_dir}]

        --publish_mode STR
            The mode used to publish files to the target directory [Default: ${defaults.publish_mode}]


    Standard params:

        --help
            Show this message and exit

        --version
            Show the pipeline version and exit
    """.stripIndent()
}
