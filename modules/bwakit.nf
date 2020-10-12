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

params.publish_bwa_index = false
params.publish_bwa_mem = false


/*
 * Processes
 */

process bwa_index {

    tag { fasta }

    label 'bwakit'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_bwa_index,
        mode: params.publish_mode,
    )

    input:
    path fasta

    output:
    path "${fasta}.{amb,ann,bwt,pac,sa}"

    """
    bwa index "${fasta}"
    """
}

process bwa_mem {

    tag { sample == readgroup ? sample : "${sample}:${readgroup}" }

    label 'bwakit'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_bwa_mem,
        mode: params.publish_mode,
    )

    input:
    tuple sample, readgroup, path(reads)
    path bwa_index

    output:
    tuple sample, path("${readgroup}.aln.bam")

    script:
    def task_cpus = task.cpus > 1 ? task.cpus - 1 : task.cpus

    def idxbase = bwa_index.first().baseName
    def fastq_files = reads.collect { /"$it"/ }.join(' ')

    """
    bwa mem \\
        -t ${task_cpus} \\
        -R '@RG\\tID:${readgroup}\\tSM:${sample}' \\
        "${idxbase}" \\
        ${fastq_files} |
    samtools view \\
        -1 \\
        -o "${readgroup}.aln.bam" \\
        -
    """
}
