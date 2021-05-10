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

params.publish_vep = false

params.vep_cache_type = null
params.vcf_info_field = 'CSQ'
params.buffer_size = 5000


/*
 * Processes
 */

process vep {

    tag { sample }

    label 'ensembl_vep'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_vep,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(indexed_vcf)
    path indexed_fasta
    path "cache/*"
    path chrom_synonyms
    val species

    output:
    path "${sample}.${species}.vep.vcf.gz{,.tbi}", emit: annotated_variants
    path "${sample}.${species}.stats.html", emit: stats

    script:
    def vep_cache_type = params.vep_cache_type in ['merged', 'refseq']
                       ? "--${params.vep_cache_type}"
                       : ''

    def synonyms = chrom_synonyms.name != 'null'
                 ? /--synonyms "${chrom_synonyms}"/
                 : ''

    """
    vep \\
        --everything \\
        --species "${species}" \\
        --input_file "${indexed_vcf.first()}" \\
        --format vcf \\
        --vcf_info_field "${params.vcf_info_field}" \\
        --vcf \\
        --output_file "${sample}.${species}.vep.vcf.gz" \\
        --stats_file "${sample}.${species}.stats.html" \\
        --warning_file "${sample}.${species}.warnings.txt" \\
        --fork ${task.cpus - 1} \\
        --dir cache \\
        ${vep_cache_type} \\
        --offline \\
        --fasta "${indexed_fasta.first()}" \\
        ${synonyms} \\
        --compress_output bgzip \\
        --allow_non_variant \\
        --buffer_size "${params.buffer_size}"
    tabix \\
        "${sample}.${species}.vep.vcf.gz"
    """
}
