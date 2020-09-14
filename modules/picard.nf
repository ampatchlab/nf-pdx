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

params.publish_mark_duplicates = false

params.max_file_handles = 1000


/*
 * Processes
 */

process mark_duplicates {

    tag { bam.name }

    label 'picard'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_mark_duplicates,
        mode: params.publish_mode,
    )

    input:
    path bam

    output:
    path "${bam.baseName}.markdup.bam{,.bai}", emit: alignments
    path "${bam.baseName}.markdup.bam.md5", emit: md5sum
    path "${bam.baseName}.metrics.txt", emit: metrics

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 0

    def Xmx = avail_mem >= 8 ? "-Xmx${avail_mem - 1}G" : ''
    def Xms = avail_mem >= 8 ? "-Xms${avail_mem.intdiv(2)}G" : ''

    """
    picard \\
        -Djava.io.tmpdir="\${PWD}/tmp" \\
        -Dpicard.useLegacyParser=false \\
        -XX:+UseSerialGC \\
        ${Xms} \\
        ${Xmx} \\
    MarkDuplicates \\
        -INPUT "${bam}" \\
        -OUTPUT "${bam.baseName}.markdup.bam" \\
        -MAX_FILE_HANDLES "${params.max_file_handles}" \\
        -METRICS_FILE "${bam.baseName}.metrics.txt" \\
        -ASSUME_SORT_ORDER coordinate \\
        -CREATE_INDEX TRUE \\
        -CREATE_MD5_FILE TRUE \\
        -VALIDATION_STRINGENCY LENIENT
    mv \\
        "${bam.baseName}.markdup.bai" \\
        "${bam.baseName}.markdup.bam.bai"
    """
}
