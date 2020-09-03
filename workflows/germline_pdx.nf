#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * ampatchlab/nf-typeHLA: HLA typing Nextflow pipeline
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


include { bwa_index } from '../modules/bwakit.nf' params( params )
include { bwa_mem } from '../modules/bwakit.nf' params( params )
include { concat_ref_genomes } from '../modules/biopython.nf' params( params )
include { cutadapt } from '../modules/cutadapt.nf' params( params )
include { fastqc } from '../modules/fastqc.nf' params( params )
include { mark_duplicates } from '../modules/picard.nf' params( params )
include { multiqc } from '../modules/multiqc.nf' params( params )
include { samtools_faidx } from '../modules/samtools.nf' params( params )
include { samtools_merge } from '../modules/samtools.nf' params( params )
include { samtools_sort } from '../modules/samtools.nf' params( params )


workflow germline_pdx {

    take:

    sample_readgroup_reads_tuples
    cutadapt_adapter_files
    human_genome_assembly
    human_ref_fasta
    human_vep_cache
    mouse_genome_assembly
    mouse_ref_fasta
    mouse_vep_cache
    multiqc_cfg

    main:

    // PRE-PROCESSING - Concatenate the reference genomes
    human_ref_inputs = [ human_genome_assembly, human_ref_fasta ]
    mouse_ref_inputs = [ mouse_genome_assembly, mouse_ref_fasta ]

    combined_ref_genome = concat_ref_genomes( human_ref_inputs, mouse_ref_inputs )


    // PRE-PROCESSING - Index the combined reference genome with BWA and Samtools
    bwa_index( combined_ref_genome )
    samtools_faidx( combined_ref_genome )


    // PRE-PROCESSING - Replace the 'sample' value with a special group key object
    sample_readgroup_reads_tuples \
        | groupTuple() \
        | map { sample, readgroups, reads ->
            tuple( groupKey(sample, readgroups.size()), readgroups, reads )
        } \
        | transpose() \
        | set { sample_readgroups }


    // STEP 1 - Run Cutadapt on the raw reads
    sample_readgroups \
        | map { sample_key, readgroup, reads -> tuple( readgroup, reads ) } \
        | set { cutadapt_inputs }

    cutadapt( cutadapt_inputs, cutadapt_adapter_files )


    // STEP 2 - Run FastQC on the raw and trimmed reads
    sample_readgroups \
        | map { sample_key, readgroup, reads -> tuple( readgroup, reads ) } \
        | mix( cutadapt.out.trimmed_reads ) \
        | flatMap { readgroup, reads -> reads } \
        | fastqc


    // STEP 3 -  Align reads using BWA-MEM
    sample_readgroups \
        | map { sample_key, readgroup, reads -> tuple( readgroup, sample_key.toString() ) } \
        | join( cutadapt.out.trimmed_reads ) \
        | map { readgroup, sample, reads -> tuple( sample, readgroup, reads ) } \
        | set { bwa_mem_inputs }

    bwa_mem( bwa_mem_inputs, bwa_index.out )


    // STEP 4 - Coordinate sort the aligned readgroups
    bwa_mem.out \
        | map { sample, aligned_readgroup -> aligned_readgroup } \
        | samtools_sort \
        | map { bam, bai -> tuple( bam.getBaseName(3), tuple( bam, bai ) ) } \
        | set { csorted_readgroups }


    // STEP 5 - Merge the aligned readgroups by sample and mark duplicates
    sample_readgroups \
        | map { sample_key, readgroup, reads -> tuple( readgroup, sample_key ) } \
        | join( csorted_readgroups ) \
        | map { readgroup, sample_key, bam_tuple -> tuple( sample_key, bam_tuple.first() ) } \
        | groupTuple() \
        | map { sample_key, bams -> tuple( sample_key.toString(), bams ) } \
        | samtools_merge \
        | map { bam, bai -> bam } \
        | mark_duplicates


    // STEP 6 - Create a MultiQC report
    Channel.empty() \
        | mix( cutadapt.out.logs ) \
        | mix( fastqc.out ) \
        | mix( mark_duplicates.out.metrics ) \
        | collect \
        | set { multiqc_inputs }

    multiqc( multiqc_inputs, multiqc_cfg )
}
