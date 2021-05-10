/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * ampatchlab/nf-pdx: Patient-derived xenograft Nextflow pipeline
 *
 * Copyright (C) 2020-2021 QIMR Berghofer Medical Research Institute
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

params.publish_split_regions = false

params.num_regions = 1000000
params.suffix_length = 5


/*
 * Processes
 */

process split_regions {

    tag { sample }

    label 'coreutils'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_split_regions,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(bed)

    output:
    tuple val(sample), path("${sample}.${/[0-9]/*params.suffix_length}.bed")

    shell:
    '''
    zcat "!{bed}" | shuf | split \\
        -a "!{params.suffix_length}" \\
        -d \\
        -l "!{params.num_regions}" \\
        --filter='LC_ALL=C sort -k1,1V -k2,2n -k3,3n > ${FILE}.bed' \\
        - \\
        "!{sample}."
    '''
}
