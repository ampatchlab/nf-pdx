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

params.publish_bcftools_subset_pass = false
params.publish_bcftools_subset_regions = false
params.publish_bcftools_mpileup = false
params.publish_bcftools_concat = false

params.mpileup_max_depth = 20000
params.mpileup_min_mq = 0
params.mpileup_min_bq = 1

params.mpileup_exclude_filters = [
  'LOW_QUAL': 'QUAL<10', // low quality

  'LOW_DP': 'FORMAT/DP<10', // low number of high quality bases
  'LOW_GQ': 'FORMAT/GQ<15', // low genotype quality

  'BQ_BIAS': 'ABS(INFO/BQBZ)>5', // base quality bias
  'MQ_BIAS': 'ABS(INFO/MQBZ)>10', // mapping quality bias
  'MQS_BIAS': 'ABS(INFO/MQSBZ)>5', // mapping quality vs strand bias
  'RP_BIAS': 'ABS(INFO/RPBZ)>5', // read position bias
  'SC_BIAS': 'ABS(INFO/SCBZ)>5', // soft-clip length bias
]

params.mpileup_include_filters = null


/*
 * Processes
 */

process bcftools_subset_pass {

    tag { sample }

    label 'bcftools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_bcftools_subset_pass,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(indexed_vcf)

    output:
    path "${sample}.pass.vcf.gz{,.tbi}"

    """
    bcftools view \\
        --no-version \\
        -f PASS \\
        -o "${sample}.pass.vcf.gz" \\
        -Oz \\
        "${indexed_vcf.first()}"
    bcftools index \\
        -t \\
        "${sample}.pass.vcf.gz"
    """
}

process bcftools_subset_regions {

    tag { sample }

    label 'bcftools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_bcftools_subset_regions,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(indexed_vcf)
    path bed

    output:
    tuple val(sample), path("${sample}.${bed.baseName}.vcf.gz{,.tbi}")

    """
    bcftools view \\
        --no-version \\
        -R "${bed}" \\
        -o "${sample}.${bed.baseName}.vcf.gz" \\
        -Oz \\
        "${indexed_vcf.first()}"
    bcftools index \\
        -t \\
        "${sample}.${bed.baseName}.vcf.gz"
    """
}

process bcftools_mpileup {

    tag { bed.baseName }

    label 'bcftools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_bcftools_mpileup,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(bed), path(indexed_bam_files)
    path indexed_fasta

    output:
    tuple val(sample), path("${bed.baseName}.vcf.gz{,.tbi}")

    script:
    def info_fields = ['AD', 'ADF', 'ADR'].collect { "INFO/$it" }
    def format_fields = ['SP', 'DP', 'AD', 'ADF', 'ADR'].collect { "FORMAT/$it" }

    def mpileup_exclude_filters = params.mpileup_exclude_filters
        .collect { name, expr -> "bcftools filter --no-version -Ou -m+ -s '${name}' -e '${expr}' - |" }
        .join('\n'+' '*4)
    def mpileup_include_filters = params.mpileup_include_filters
        .collect { name, expr -> "bcftools filter --no-version -Ou -m+ -s '${name}' -i '${expr}' - |" }
        .join('\n'+' '*4)

    def bam_file_list = indexed_bam_files
        .findAll { infile -> infile.name.endsWith('.bam') }
        .collect { /"${it}"/ }
        .join(' \\\n'+' '*8)

    """
    bcftools mpileup \\
        --no-version \\
        -Ou \\
        -d "${params.mpileup_max_depth}" \\
        -f "${indexed_fasta.first()}" \\
        -q "${params.mpileup_min_mq}" \\
        -Q "${params.mpileup_min_bq}" \\
        -R "${bed}" \\
        -a "${info_fields.join(',')},${format_fields.join(',')}" \\
        ${bam_file_list} |
    bcftools norm \\
        --no-version \\
        -Ou \\
        -f "${indexed_fasta.first()}" \\
        -m +any \\
        - |
    bcftools call \\
        --no-version \\
        -Ou \\
        -m \\
        -v \\
        -a GQ,GP \\
        - |
    bcftools norm \\
        --no-version \\
        -Ou \\
        -f "${indexed_fasta.first()}" \\
        -m -any \\
        - |
    ${mpileup_exclude_filters ?: '\\'}
    ${mpileup_include_filters ?: '\\'}
    bcftools view \\
        --no-version \\
        -Oz \\
        -o "${bed.baseName}.vcf.gz" \\
        -
    bcftools index \\
        -t \\
        "${bed.baseName}.vcf.gz"
    """
}

process bcftools_concat {

    tag { sample }

    label 'bcftools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_bcftools_concat,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(indexed_vcf_files)

    output:
    path "${sample}.vcf.gz{,.tbi}"

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toGiga()}G" : ''

    def vcf_file_list = indexed_vcf_files
        .findAll { infile ->
            [ '.bcf', '.bcf.gz', '.vcf', '.vcf.gz' ].any { infile.name.endsWith(it) }
        }
        .sort { it.name }
        .collect { /"${it}"/ }
        .join(' \\\n'+' '*8)

    """
    bcftools concat \\
        --no-version \\
        -a \\
        -D \\
        ${vcf_file_list} |
    bcftools sort \\
        ${avail_mem} \\
        -Oz \\
        -o "${sample}.vcf.gz" \\
        -T . \\
        -
    bcftools index \\
        -t \\
        "${sample}.vcf.gz"
    """
}
