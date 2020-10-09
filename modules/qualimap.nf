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

params.publish_qualimap = false

params.num_reads = 1000
params.num_windows = 400


/*
 * Processes
 */

process qualimap {

    tag { bam.name }

    label 'qualimap'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_qualimap,
        mode: params.publish_mode,
    )

    input:
    path bam
    path gff

    output:
    path "${bam.baseName}"

    script:
    def avail_mem = task.memory ? "${task.memory.toGiga()}g" : "1200m"
    def JAVA_OPTS = "-XX:+UseSerialGC -Xms32m -Xmx${avail_mem}"

    def feature_file = gff.name != 'null' ? /--feature-file "${gff}"/ : ''

    """
    export JAVA_OPTS="${JAVA_OPTS}"
    qualimap bamqc \\
        -bam "${bam}" \\
        -nr ${params.num_reads} \\
        -nt ${task.cpus} \\
        -nw ${params.num_windows} \\
        -outdir "${bam.baseName}" \\
        ${feature_file} \\
        --paint-chromosome-limits \\
        --collect-overlap-pairs \\
        --skip-duplicated
    """
}
