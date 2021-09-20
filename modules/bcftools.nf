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

params.mpileup_max_depth = 20000
params.mpileup_min_mq = 0
params.mpileup_min_bq = 1

params.exclude_filters = [
  'LOW_QUAL': 'QUAL<10', // low quality

  'LOW_DP': 'FORMAT/DP<10', // low number of high quality bases
  'LOW_GQ': 'FORMAT/GQ<15', // low genotype quality

  'BQ_BIAS': 'ABS(INFO/BQBZ)>5', // base quality bias
  'MQ_BIAS': 'ABS(INFO/MQBZ)>10', // mapping quality bias
  'MQS_BIAS': 'ABS(INFO/MQSBZ)>5', // mapping quality vs strand bias
  'RP_BIAS': 'ABS(INFO/RPBZ)>5', // read position bias
  'SC_BIAS': 'ABS(INFO/SCBZ)>5', // soft-clip length bias
]

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
        -o "${bed.baseName}.bcf.gz" \\
        -Ob \\
        -d "${params.mpileup_max_depth}" \\
        -f "${indexed_fasta.first()}" \\
        -q "${params.mpileup_min_mq}" \\
        -Q "${params.mpileup_min_bq}" \\
        -R "${bed}" \\
        -a "${info_fields.join(',')},${format_fields.join(',')}" \\
        -b "${sample}.list"
    bcftools index \\
        -c \\
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
    def avail_mem = task.memory ? "-m ${task.memory.toGiga()}G" : ''

    def exclude_filters = params.exclude_filters.collect { name, expr ->

        """\
        >&2 echo "Applying '${name}' exclude expression..."
        input="\${output}"
        output='filters/exclude/${sample}.${name}.bcf.gz'

        bcftools filter \\
            --no-version \\
            -o "\${output}" \\
            -Ob \\
            -m+ \\
            -s '${name}' \\
            -e '${expr}' \\
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
            -o "\${output}" \\
            -Ob \\
            -m+ \\
            -s '${name}' \\
            -i '${expr}' \\
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
        -o "\${output}" \\
        -Ob \\
        -a \\
        -D \\
        -f "\${input}"


    >&2 echo "Sorting variants..."
    input="\${output}"
    output="${sample}.sort.bcf.gz"

    bcftools sort \\
        ${avail_mem} \\
        -o "\${output}" \\
        -Ob \\
        -T . \\
        "\${input}"


    >&2 echo "Joining biallelic sites into multiallelic records..."
    input="\${output}"
    output="${sample}.join.bcf.gz"

    bcftools norm \\
        --no-version \\
        -o "\${output}" \\
        -Ob \\
        -f "${indexed_fasta.first()}" \\
        -m +any \\
        "\${input}"


    >&2 echo "Calling SNP/indels using the multiallelic-caller..."
    input="\${output}"
    output="${sample}.call.bcf.gz"

    bcftools call \\
        --no-version \\
        -o "\${output}" \\
        -Ob \\
        -m \\
        -v \\
        -a GQ,GP \\
        "\${input}"


    >&2 echo "Splitting multiallelic sites into biallelic records..."
    input="\${output}"
    output="${sample}.split.bcf.gz"

    bcftools norm \\
        --no-version \\
        -o "\${output}" \\
        -Ob \\
        -f "${indexed_fasta.first()}" \\
        -m -any \\
        "\${input}"


    ${exclude_filters.join('\n\n').trim() ?: '# No exclude filters'}


    ${include_filters.join('\n\n').trim() ?: '# No include filters'}


    >&2 echo "Converting BCF to VCF..."
    input="\${output}"
    output="${sample}.vcf.gz"

    bcftools view \\
        --no-version \\
        -Oz \\
        -o "\${output}" \\
        "\${input}"
    bcftools index \\
        -t \\
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
    path "${sample}.sorted.vcf.gz{,.tbi}"

    script:
    def avail_mem = task.memory ? "-m ${task.memory.toGiga()}G" : ''

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
        -o "${sample}.bcf.gz" \\
        -Ob \\
        -a \\
        -D \\
        -f "${sample}.list"
    bcftools sort \\
        ${avail_mem} \\
        -Oz \\
        -o "${sample}.sorted.vcf.gz" \\
        -T . \\
        "${sample}.bcf.gz"
    bcftools index \\
        -t \\
        "${sample}.sorted.vcf.gz"
    """
}
