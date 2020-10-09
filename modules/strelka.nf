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


/*
 * Params
 */

params.publish_dir = './results'
params.publish_everything = false
params.publish_mode = 'copy'

params.publish_strelka_germline_wf = false
params.publish_strelka_somatic_wf = false

params.exome = false


/*
 * Processes
 */

process strelka_germline_wf {

    tag { sample }

    label 'strelka'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_strelka_germline_wf,
        mode: params.publish_mode,
        saveAs: { fn ->
            [sample, fn.substring(fn.lastIndexOf('/')+1)].join('.')
        }
    )

    input:
    tuple val(sample), path(indexed_bam_files)
    path indexed_fasta
    path bed

    output:
    tuple val(sample), path("results/variants/variants.vcf.gz{,.tbi}"), emit: variants
    tuple val(sample), path("results/stats/runStats.{tsv,xml}"), emit: stats

    script:
    def call_regions = bed.name != 'null' ? /--callRegions "${bed}"/ : ''
    def exome = params.exome ? '--exome' : ''

    def avail_mem = task.memory ? "-g ${task.memory.toGiga()}" : 'unlimited'

    def bam_file_list = indexed_bam_files
        .findAll { infile -> infile.name.endsWith('.bam') }
        .collect { /--bam "${it}"/ }
        .join(' \\\n'+' '*8)

    """
    configureStrelkaGermlineWorkflow.py \\
        ${bam_file_list} \\
        --referenceFasta "${indexed_fasta.first()}" \\
        --runDir . \\
        ${call_regions} \\
        ${exome}
    python runWorkflow.py \\
        -m local \\
        -j ${task.cpus} \\
        ${avail_mem}
    """
}

process strelka_somatic_wf {

    tag { analysis }

    label 'strelka'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_strelka_somatic_wf,
        mode: params.publish_mode,
        saveAs: { fn ->
            [analysis, fn.substring(fn.lastIndexOf('/')+1)].join('.')
        }
    )

    input:
    tuple val(analysis), path(indexed_test_bam), path(indexed_control_bam), path(indexed_indel_candidates)
    path indexed_fasta
    path bed

    output:
    tuple val(analysis), path("results/variants/somatic.indels.vcf.gz{,.tbi}"), emit: somatic_indels
    tuple val(analysis), path("results/variants/somatic.snvs.vcf.gz{,.tbi}"), emit: somatic_snvs
    tuple val(analysis), path("results/stats/runStats.{tsv,xml}"), emit: stats

    script:
    def call_regions = bed.name != 'null' ? /--callRegions "${bed}"/ : ''
    def exome = params.exome ? '--exome' : ''

    def avail_mem = task.memory ? "-g ${task.memory.toGiga()}" : 'unlimited'

    """
    configureStrelkaSomaticWorkflow.py \\
        --tumorBam "${indexed_test_bam.first()}" \\
        --normalBam "${indexed_control_bam.first()}" \\
        --referenceFasta "${indexed_fasta.first()}" \\
        --indelCandidates "${indexed_indel_candidates.first()}" \\
        --runDir . \\
        ${call_regions} \\
        ${exome}
    python runWorkflow.py \\
        -m local \\
        -j ${task.cpus} \\
        ${avail_mem}
    """
}
