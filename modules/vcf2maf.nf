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

params.publish_vcf2maf = false


/*
 * Processes
 */

process vcf2maf {

    tag { input_vcf.name }

    label 'vcf2maf'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_vcf2maf,
        mode: params.publish_mode,
    )

    input:
    tuple path(input_vcf), val(test_sample), val(control_sample)
    path indexed_fasta
    val species
    val ncbi_build

    output:
    path "${input_vcf.baseName}.maf"

    """
    vcf2maf.pl \\
        --input-vcf "${input_vcf}" \\
        --output-maf "${input_vcf.baseName}.maf" \\
        --tmp-dir tmp \\
        --tumor-id "${test_sample}" \\
        --normal-id "${control_sample}" \\
        --inhibit-vep \\
        --ref-fasta "${indexed_fasta.first()}" \\
        --species "${species}" \\
        --ncbi-build "${ncbi_build}"
    """
}
