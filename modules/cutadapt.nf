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

params.publish_cutadapt = false

params.min_read_length = 0
params.base_qual_cutoff = null


/*
 * Processes
 */

process cutadapt {

    tag { readgroup }

    label 'cutadapt'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_cutadapt,
        mode: params.publish_mode,
    )

    input:
    tuple val(readgroup), path(input_reads)
    path adapter_files

    output:
    tuple val(readgroup), path(output_reads), emit: trimmed_reads
    path "${readgroup}.log", emit: logs

    script:
    def (in1_fq, in2_fq) = input_reads

    def in1_str = in1_fq ? /"$in1_fq"/ : ''
    def in2_str = in2_fq ? /"$in2_fq"/ : ''

    // adapter params
    def (r1_file, r2_file) = adapter_files

    def r1_adapters = in1_fq && r1_file?.name != 'null-1.fa' ? /-a "file:${r1_file}"/ : ''
    def r2_adapters = in2_fq && r2_file?.name != 'null-2.fa' ? /-A "file:${r2_file}"/ : ''

    // output params
    def (out1_fq, out2_fq) = output_reads = input_reads.collect {
        get_fastq_basename( it ) + '.trimmed.fastq.gz'
    }

    def out1_str = in1_fq && out1_fq ? /-o "${out1_fq}"/ : ''
    def out2_str = in2_fq && out2_fq ? /-p "${out2_fq}"/ : ''

    // filtering and trimming params
    def min_read_length = params.min_read_length ? "-m ${params.min_read_length}" : ''
    def base_qual_cutoff = params.base_qual_cutoff ? "-q ${params.base_qual_cutoff}" : ''

    """
    cutadapt \\
        -j "${task.cpus}" \\
        -Z \\
        ${r1_adapters} \\
        ${r2_adapters} \\
        ${min_read_length} \\
        ${base_qual_cutoff} \\
        --trim-n \\
        ${out1_str} \\
        ${out2_str} \\
        ${in1_str} \\
        ${in2_str} \\
        > "${readgroup}.log"
    """
}


/*
 * Functions
 */

def get_fastq_basename( fastq ) {

    def compression_extns = [ '.gz', '.bz2', '.xz' ]
    def fastq_extns = [ '.fastq', '.fq' ]

    def fn = fastq.name

    for( List<String> extns : [ compression_extns, fastq_extns ] ) {

        def extn = extns.find { fn.endsWith( it ) }

        fn = fn.take( fn.size() - ( extn ? extn.size() : 0 ) )
    }

    return fn
}
