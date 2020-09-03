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

params.publish_samtools_faidx = false
params.publish_samtools_sort = false
params.publish_samtools_merge = false
params.publish_samtools_stats = false


/*
 * Processes
 */

process samtools_faidx {

    tag { fasta.name }

    label 'samtools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_samtools_faidx,
        mode: params.publish_mode,
    )

    input:
    path fasta

    output:
    path "${fasta}.fai"

    """
    samtools faidx "${fasta}"
    """
}

process samtools_sort {

    tag { bam.name }

    label 'samtools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_samtools_sort,
        mode: params.publish_mode,
    )

    input:
    path bam

    output:
    path "${bam.baseName}.sorted.bam{,.bai}"

    script:
    def tmp_prefix = "${bam.baseName}.sorted"
    def out_files = [ "${tmp_prefix}.bam", "${tmp_prefix}.bam.bai" ].join('##idx##')

    def avail_mem = task.memory ? task.memory.toGiga().intdiv(task.cpus) : 0
    def mem_per_thread = avail_mem ? "-m ${avail_mem}G" : ''

    """
    samtools sort \\
        --write-index \\
        ${mem_per_thread} \\
        -o "${out_files}" \\
        -T "${tmp_prefix}" \\
        -@ "${task.cpus - 1}" \\
        "${bam}"
    """
}

process samtools_merge {

    tag { sample }

    label 'samtools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_samtools_merge,
        mode: params.publish_mode,
    )

    input:
    tuple sample, path(bams)

    output:
    path "${sample}.bam{,.bai}"

    script:
    def input_bams = bams.collect { /"$it"/ }.join(' ')
    def out_files = [ "${sample}.bam", "${sample}.bam.bai" ].join('##idx##')

    """
    samtools merge \\
        --write-index \\
        "${out_files}" \\
        ${input_bams}
    """
}

process samtools_stats {

    tag { bam.name }

    label 'samtools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_samtools_stats,
        mode: params.publish_mode,
    )

    input:
    path bam

    output:
    path "${bam.baseName}.{stats,flagstat,idxstats}"

    """
    samtools stats "${bam}" > "${bam.baseName}.stats"
    samtools flagstat "${bam}" > "${bam.baseName}.flagstat"
    samtools idxstats "${bam}" > "${bam.baseName}.idxstats"
    """
}
