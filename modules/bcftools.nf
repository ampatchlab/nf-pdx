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
params.publish_bcftools_call = false
params.publish_bcftools_concat = false

// mpileup input options
params.mpileup_max_depth = 250
params.mpileup_min_mq = 0
params.mpileup_min_bq = 1

// mpileup genotype likelihoods options
params.mpileup_ext_prob = 20
params.mpileup_gap_frac = 0.05
params.mpileup_tandem_qual = 500
params.mpileup_skip_indels = false
params.mpileup_max_idepth = 250
params.mpileup_min_ireads = 2
params.mpileup_max_read_len = 500
params.mpileup_open_prob = 40
params.mpileup_per_sample_mF = false
params.mpileup_indel_bias = 1.00

// variant calling filters
params.exclude_filters = null
params.include_filters = null


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
    tuple val(sample), path("${sample}.pass.vcf.gz{,.tbi}")

    """
    bcftools view \\
        --no-version \\
        --output "${sample}.pass.vcf.gz" \\
        --output-type z \\
        --apply-filters PASS \\
        "${indexed_vcf.first()}"
    bcftools index \\
        --tbi \\
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
        --output "${sample}.${bed.baseName}.vcf.gz" \\
        --output-type z \\
        --regions-file "${bed}" \\
        "${indexed_vcf.first()}"
    bcftools index \\
        --tbi \\
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
    tuple val(sample), path("${bed.baseName}.bcf.gz{,.csi}")

    script:
    def info_fields = ['AD', 'ADF', 'ADR'].collect { "INFO/$it" }
    def format_fields = ['SP', 'DP', 'AD', 'ADF', 'ADR'].collect { "FORMAT/$it" }

    def bam_list = indexed_bam_files.findAll { it.name.endsWith('.bam') }

    """
    cat << __EOF__ > "${sample}.list"
    ${bam_list.join('\n'+' '*4)}
    __EOF__

    bcftools mpileup \\
        --no-version \\
        --output "${bed.baseName}.bcf.gz" \\
        --output-type b \\
        --max-depth "${params.mpileup_max_depth}" \\
        --fasta-ref "${indexed_fasta.first()}" \\
        --min-MQ "${params.mpileup_min_mq}" \\
        --min-BQ "${params.mpileup_min_bq}" \\
        --regions-file "${bed}" \\
        --annotate "${info_fields.join(',')},${format_fields.join(',')}" \\
        --bam-list "${sample}.list" \\
        --ext-prob "${params.mpileup_ext_prob}" \\
        --gap-frac "${params.mpileup_gap_frac}" \\
        --tandem-qual "${params.mpileup_tandem_qual}" \\
        ${params.mpileup_skip_indels ? '--skip-indels' : ''} \\
        --max-idepth "${params.mpileup_max_idepth}" \\
        --min-ireads "${params.mpileup_min_ireads}" \\
        --max-read-len "${params.mpileup_max_read_len}" \\
        --open-prob "${params.mpileup_open_prob}" \\
        ${params.mpileup_per_sample_mF ? '--per-sample-mF' : ''} \\
        --indel-bias "${params.mpileup_indel_bias}"
    bcftools index \\
        --csi \\
        "${bed.baseName}.bcf.gz"
    """
}

process bcftools_call {

    tag { sample }

    label 'bcftools'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_bcftools_call,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(indexed_vcf_files)
    path indexed_fasta

    output:
    tuple val(sample), path("${sample}.vcf.gz{,.tbi}")

    script:
    def avail_mem = task.memory ? "--max-mem ${task.memory.toGiga()}G" : ''

    def exclude_filters = params.exclude_filters.collect { name, expr ->

        """\
        >&2 echo "Applying '${name}' exclude expression..."
        input="\${output}"
        output='filters/exclude/${sample}.${name}.bcf.gz'

        bcftools filter \\
            --no-version \\
            --output "\${output}" \\
            --output-type b \\
            --mode + \\
            --soft-filter '${name}' \\
            --exclude '${expr}' \\
            "\${input}"
        """.stripIndent().replaceAll(/(?m)^/, ' '*4)
    }

    def include_filters = params.include_filters.collect { name, expr ->

        """\
        >&2 echo "Applying '${name}' include expression..."
        input="\${output}"
        output='filters/include/${sample}.${name}.bcf.gz'

        bcftools filter \\
            --no-version \\
            --output "\${output}" \\
            --output-type b \\
            --mode + \\
            --soft-filter '${name}' \\
            --include '${expr}' \\
            "\${input}"
        """.stripIndent().replaceAll(/(?m)^/, ' '*4)
    }

    def vcf_list = indexed_vcf_files
        .findAll { infile ->
            [ '.bcf', '.bcf.gz', '.vcf', '.vcf.gz' ].any { infile.name.endsWith(it) }
        }
        .sort { it.name }

    """
    mkdir -p filters/{exclude,include}

    output="${sample}.list"
    cat << __EOF__ > "\${output}"
    ${vcf_list.join('\n'+' '*4)}
    __EOF__


    >&2 echo "Concatenating ${vcf_list.size()} VCF/BCF files..."
    input="\${output}"
    output="${sample}.concat.bcf.gz"

    bcftools concat \\
        --no-version \\
        --output "\${output}" \\
        --output-type b \\
        --allow-overlaps \\
        --remove-duplicates \\
        --file-list "\${input}"


    >&2 echo "Sorting variants..."
    input="\${output}"
    output="${sample}.sort.bcf.gz"

    bcftools sort \\
        ${avail_mem} \\
        --output "\${output}" \\
        --output-type b \\
        --temp-dir . \\
        "\${input}"


    >&2 echo "Joining biallelic sites into multiallelic records..."
    input="\${output}"
    output="${sample}.join.bcf.gz"

    bcftools norm \\
        --no-version \\
        --output "\${output}" \\
        --output-type b \\
        --fasta-ref "${indexed_fasta.first()}" \\
        --multiallelics +any \\
        "\${input}"


    >&2 echo "Calling SNP/indels using the multiallelic-caller..."
    input="\${output}"
    output="${sample}.call.bcf.gz"

    bcftools call \\
        --no-version \\
        --output "\${output}" \\
        --output-type b \\
        --multiallelic-caller \\
        --variants-only \\
        --annotate GQ,GP \\
        "\${input}"


    >&2 echo "Splitting multiallelic sites into biallelic records..."
    input="\${output}"
    output="${sample}.split.bcf.gz"

    bcftools norm \\
        --no-version \\
        --output "\${output}" \\
        --output-type b \\
        --fasta-ref "${indexed_fasta.first()}" \\
        --multiallelics -any \\
        "\${input}"


    ${exclude_filters.join('\n\n').trim() ?: '# No exclude filters'}


    ${include_filters.join('\n\n').trim() ?: '# No include filters'}


    >&2 echo "Converting BCF to VCF..."
    input="\${output}"
    output="${sample}.vcf.gz"

    bcftools view \\
        --no-version \\
        --output "\${output}" \\
        --output-type z \\
        "\${input}"
    bcftools index \\
        --tbi \\
        "\${output}"
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
    tuple val(sample), path("${sample}.sorted.vcf.gz{,.tbi}")

    script:
    def avail_mem = task.memory ? "--max-mem ${task.memory.toGiga()}G" : ''

    def vcf_list = indexed_vcf_files
        .findAll { infile ->
            [ '.bcf', '.bcf.gz', '.vcf', '.vcf.gz' ].any { infile.name.endsWith(it) }
        }
        .sort { it.name }

    """
    cat << __EOF__ > "${sample}.list"
    ${vcf_list.join('\n'+' '*4)}
    __EOF__

    bcftools concat \\
        --no-version \\
        --output "${sample}.bcf.gz" \\
        --output-type b \\
        --allow-overlaps \\
        --remove-duplicates \\
        --file-list "${sample}.list"
    bcftools sort \\
        ${avail_mem} \\
        --output "${sample}.sorted.vcf.gz" \\
        --output-type z \\
        --temp-dir . \\
        "${sample}.bcf.gz"
    bcftools index \\
        --tbi \\
        "${sample}.sorted.vcf.gz"
    """
}
