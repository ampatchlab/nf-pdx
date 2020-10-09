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
include { manta_somatic_wf } from '../modules/manta.nf' params( params )
include { split_regions } from '../modules/coreutils.nf' params( params )
include { strelka_germline_wf } from '../modules/strelka.nf' params( params )
include { strelka_somatic_wf } from '../modules/strelka.nf' params( params )

include { bcftools_concat as vcf_concat } from '../modules/bcftools.nf' params( params )



workflow germline_variant_calling {

    take:

    germline_analysis_tuples
    indexed_ref_fasta
    strelka_call_regions


    main:

    // STEP 1 - Call variants using Strelka
    strelka_germline_wf(
        germline_analysis_tuples,
        indexed_ref_fasta,
        strelka_call_regions,
    )


    // STEP 2 - Chunkify all PASS variants in BED format
    strelka_germline_wf.out.variants \
        | bcftools_subset_pass \
        | map { vcf, tbi -> tuple( vcf.getBaseName(2), tuple( vcf, tbi ) ) } \
        | bedops_convert2bed \
        | map { bed -> tuple( bed.getBaseName(3), bed ) } \
        | split_regions \
        | join( germline_analysis_tuples ) \
        | map { sample, beds, indexed_bam_files ->
            def bed_files = [beds].flatten()
            def group_key = groupKey( sample, bed_files.size() )

            tuple( group_key, bed_files, indexed_bam_files )
        } \
        | transpose( by: 1 ) \
        | set { mpileup_inputs }


    // STEP 3 - Pileup and re-call each chunk of variant positions using BCFtools
    bcftools_mpileup( mpileup_inputs, indexed_ref_fasta ) \
        | groupTuple() \
        | map { sample, vcf_tuples -> tuple( sample.toString(), vcf_tuples.flatten() ) } \
        | bcftools_concat


    emit:

    bcftools_concat.out.map { vcf, tbi -> tuple( vcf.getBaseName(2), tuple( vcf, tbi ) ) }
}


workflow somatic_variant_calling {

    take:

    somatic_analysis_tuples
    indexed_ref_fasta
    manta_call_regions
    strelka_call_regions


    main:

    // STEP 1 - Run Manta to call candidate indels
    manta_somatic_wf(
        somatic_analysis_tuples,
        indexed_ref_fasta,
        manta_call_regions,
    )


    // STEP 2 - Call variants using Strelka
    somatic_analysis_tuples \
        | join( manta_somatic_wf.out.candidate_small_indels ) \
        | set { strelka_inputs }

    strelka_somatic_wf(
        strelka_inputs,
        indexed_ref_fasta,
        strelka_call_regions,
    )


    // STEP 3 - Chunkify all PASS variants in BED format
    strelka_somatic_wf.out.somatic_snvs \
        | join( strelka_somatic_wf.out.somatic_indels ) \
        | map { analysis_id, indexed_snvs_vcf, indexed_indels_vcf ->
            tuple( analysis_id, [indexed_snvs_vcf, indexed_indels_vcf].flatten() )
        } \
        | vcf_concat \
        | map { vcf, tbi -> tuple( vcf.getBaseName(2), tuple( vcf, tbi ) ) } \
        | bcftools_subset_pass \
        | map { vcf, tbi -> tuple( vcf.getBaseName(2), tuple( vcf, tbi ) ) } \
        | bedops_convert2bed \
        | map { bed -> tuple( bed.getBaseName(3), bed ) } \
        | split_regions \
        | join( somatic_analysis_tuples ) \
        | map { analysis_id, beds, indexed_test_bam, indexed_control_bam ->
            def bed_files = [beds].flatten()
            def group_key = groupKey( analysis_id, bed_files.size() )

            tuple( group_key, bed_files, [ *indexed_test_bam, *indexed_control_bam ] )
        } \
        | transpose( by: 1 ) \
        | set { mpileup_inputs }


    // STEP 4 - Pileup and re-call each chunk of variant positions using BCFtools
    bcftools_mpileup( mpileup_inputs, indexed_ref_fasta ) \
        | groupTuple() \
        | map { analysis_id, vcf_tuples -> tuple( analysis_id.toString(), vcf_tuples.flatten() ) } \
        | bcftools_concat


    emit:

    bcftools_concat.out.map { vcf, tbi -> tuple( vcf.getBaseName(2), tuple( vcf, tbi ) ) }
}
