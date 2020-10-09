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

params.publish_manta_somatic_wf = false

params.exome = false


/*
 * Processes
 */

process manta_somatic_wf {

    tag { analysis }

    label 'manta'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_manta_somatic_wf,
        mode: params.publish_mode,
        saveAs: { fn ->
            [analysis, fn.substring(fn.lastIndexOf('/')+1)].join('.')
        }
    )

    input:
    tuple val(analysis), path(indexed_test_bam), path(indexed_control_bam)
    path indexed_fasta
    path bed

    output:
    tuple val(analysis), path("results/variants/candidateSmallIndels.vcf.gz{,.tbi}"), emit: candidate_small_indels
    tuple val(analysis), path("results/variants/candidateSV.vcf.gz{,.tbi}"), emit: candidate_sv
    tuple val(analysis), path("results/variants/diploidSV.vcf.gz{,.tbi}"), emit: diploid_sv
    tuple val(analysis), path("results/variants/somaticSV.vcf.gz{,.tbi}"), emit: somatic_sv
    tuple val(analysis), path("results/stats/alignmentStatsSummary.txt"), emit: alignment_stats_summary
    tuple val(analysis), path("results/stats/svCandidateGenerationStats.{tsv,xml}"), emit: sv_candidate_generation_stats
    tuple val(analysis), path("results/stats/svLocusGraphStats.tsv"), emit: sv_locus_graph_stats

    script:
    def call_regions = bed.name != 'null' ? /--callRegions "${bed}"/ : ''
    def exome = params.exome ? '--exome' : ''

    def avail_mem = task.memory ? "-g ${task.memory.toGiga()}" : 'unlimited'

    """
    configManta.py \\
        --runDir . \\
        --tumorBam "${indexed_test_bam.first()}" \\
        --normalBam "${indexed_control_bam.first()}" \\
        --referenceFasta "${indexed_fasta.first()}" \\
        ${call_regions} \\
        ${exome}
    python runWorkflow.py \\
        -m local \\
        -j ${task.cpus} \\
        ${avail_mem}
    """
}
