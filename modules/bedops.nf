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

params.publish_bedops_convert2bed = false


/*
 * Processes
 */

process bedops_convert2bed {

    tag { sample }

    label 'bedops'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_bedops_convert2bed,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(indexed_vcf)

    output:
    path "${sample}.bed.gz"

    """
    zcat "${indexed_vcf.first()}" |
        convert2bed -i vcf -d - |
        cut -f -3 |
        gzip \\
        > "${sample}.bed.gz"
    """
}
