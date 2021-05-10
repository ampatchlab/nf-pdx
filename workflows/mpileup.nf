#!/usr/bin/env nextflow

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


include { bcftools_concat } from '../modules/bcftools.nf' params( params )
include { bcftools_subset_pass } from '../modules/bcftools.nf' params( params )
include { bcftools_mpileup } from '../modules/bcftools.nf' params( params )
include { bedops_convert2bed } from '../modules/bedops.nf' params( params )
include { split_regions } from '../modules/coreutils.nf' params( params )



workflow mpileup {

    take:

    analysis_tuples
    vcf_inputs
    indexed_ref_fasta


    main:

    // STEP 1 - Chunkify variants in BED format
    bedops_convert2bed( vcf_inputs ) \
        | map { bed -> tuple( bed.getBaseName(2), bed ) } \
        | split_regions \
        | join( analysis_tuples ) \
        | map { sample, beds, indexed_test_bam, indexed_control_bam = [] ->
            def bed_files = [beds].flatten()
            def group_key = groupKey( sample, bed_files.size() )

            tuple( group_key, bed_files, [ *indexed_test_bam, *indexed_control_bam ] )
        } \
        | transpose( by: 1 ) \
        | set { mpileup_inputs }


    // STEP 2 - Pileup and re-call each chunk of variant positions
    bcftools_mpileup( mpileup_inputs, indexed_ref_fasta ) \
        | groupTuple() \
        | map { sample, vcf_tuples -> tuple( sample.toString(), vcf_tuples.flatten() ) } \
        | bcftools_concat


    emit:

    bcftools_concat.out.map { vcf, tbi -> tuple( vcf.getBaseName(2), tuple( vcf, tbi ) ) }
}
