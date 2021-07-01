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

params.publish_unpack_vep_cache = false
params.publish_vep = false

params.vep_cache_type = null
params.vcf_info_field = 'CSQ'
params.buffer_size = 5000


/*
 * Processes
 */

process unpack_vep_cache {

    tag { vep_cache.name }

    label 'ensembl_vep'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_unpack_vep_cache,
        mode: params.publish_mode,
    )

    input:
    path vep_cache

    output:
    path "cache/*/*/info.txt", emit: cache_info
    path "cache", emit: cache_dir

    """
    mkdir cache
    tar -xf "${vep_cache}" -C cache
    """
}

process vep {

    tag { "${sample}:${assembly}" }

    label 'ensembl_vep'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_vep,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(indexed_vcf)
    path indexed_fasta
    path cache_dir
    path chrom_synonyms
    val species
    val assembly

    output:
    path "${sample}.${assembly}.vep.vcf.gz{,.tbi}", emit: indexed_vcf
    path "${sample}.${assembly}.stats.html", emit: stats_html

    script:
    def vep_cache_type = params.vep_cache_type in ['merged', 'refseq']
                       ? "--${params.vep_cache_type}"
                       : ''

    def synonyms = chrom_synonyms.name != 'null'
                 ? /--synonyms "${chrom_synonyms}"/
                 : ''

    def pick_order = 'canonical,tsl,biotype,rank,ccds,length'

    """
    vep \\
        --everything \\
        --species "${species}" \\
        --assembly "${assembly}" \\
        --input_file "${indexed_vcf.first()}" \\
        --format vcf \\
        --vcf_info_field "${params.vcf_info_field}" \\
        --vcf \\
        --output_file "${sample}.${assembly}.vep.vcf.gz" \\
        --stats_file "${sample}.${assembly}.stats.html" \\
        --warning_file "${sample}.${assembly}.warnings.txt" \\
        --fork ${task.cpus - 1} \\
        --dir "${cache_dir}" \\
        ${vep_cache_type} \\
        --offline \\
        --fasta "${indexed_fasta.first()}" \\
        ${synonyms} \\
        --compress_output bgzip \\
        --allow_non_variant \\
        --buffer_size "${params.buffer_size}" \\
        --check_existing \\
        --total_length \\
        --allele_number \\
        --no_escape \\
        --xref_refseq \\
        --failed 1 \\
        --flag_pick_allele \\
        --pick_order "${pick_order}"
    tabix \\
        "${sample}.${assembly}.vep.vcf.gz"
    """
}
