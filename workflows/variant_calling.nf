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


include { bcftools_subset_pass } from '../modules/bcftools.nf' params( params )
include { manta_somatic_wf } from '../modules/manta.nf' params( params )
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


    // STEP 2 - Extract all PASS variants
    bcftools_subset_pass( strelka_germline_wf.out.variants )


    emit:

    bcftools_subset_pass.out.map { vcf, tbi ->
        tuple( vcf.getBaseName(3), tuple( vcf, tbi ) )
    }
}


workflow somatic_variant_calling {

    take:

    somatic_analysis_tuples
    indexed_ref_fasta
    manta_call_regions
    strelka_call_regions
    skip_manta


    main:

    if( !skip_manta ) {

        // STEP 1 - Call candidate small indels using Manta
        manta_somatic_wf(
            somatic_analysis_tuples,
            indexed_ref_fasta,
            manta_call_regions,
        )

        strelka_inputs = somatic_analysis_tuples \
            | join( manta_somatic_wf.out.candidate_small_indels )
    }
    else {

        candidate_small_indels = [
            "${baseDir}/assets/null-1",
            "${baseDir}/assets/null-2",
        ]

        strelka_inputs = somatic_analysis_tuples \
            | map { tuple( *it, candidate_small_indels ) }
    }


    // STEP 2 - Call variants using Strelka
    strelka_somatic_wf(
        strelka_inputs,
        indexed_ref_fasta,
        strelka_call_regions,
    )


    // STEP 3 - Extract all PASS variants
    strelka_somatic_wf.out.somatic_snvs \
        | join( strelka_somatic_wf.out.somatic_indels ) \
        | map { analysis_id, indexed_snvs_vcf, indexed_indels_vcf ->
            tuple( analysis_id, [indexed_snvs_vcf, indexed_indels_vcf].flatten() )
        } \
        | vcf_concat \
        | map { vcf, tbi -> tuple( vcf.getBaseName(2), tuple( vcf, tbi ) ) } \
        | bcftools_subset_pass \


    emit:

    bcftools_subset_pass.out.map { vcf, tbi ->
        tuple( vcf.getBaseName(3), tuple( vcf, tbi ) )
    }
}
