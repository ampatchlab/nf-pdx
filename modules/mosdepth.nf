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

params.publish_mosdepth = false


/*
 * Processes
 */

process mosdepth {

    tag { indexed_bam.first().name }

    label 'mosdepth'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_mosdepth,
        mode: params.publish_mode,
    )

    input:
    path indexed_bam
    path bed

    output:
    path "${indexed_bam.first().baseName}.mosdepth.{global,region}.dist.txt", emit: dists
    path "${indexed_bam.first().baseName}.mosdepth.summary.txt", emit: summary
    path "${indexed_bam.first().baseName}.regions.bed.gz{,.csi}", emit: regions_bed

    script:
    def bam = indexed_bam.first()
    def by = bed.name != 'null' ? /-b "${bed}"/ : ''

    """
    mosdepth \\
        -t "${task.cpus - 1}" \\
        ${by} \\
        -n \\
        "${bam.baseName}" \\
        "${bam}"
    """
}
