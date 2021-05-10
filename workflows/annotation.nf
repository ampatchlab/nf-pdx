#!/usr/bin/env nextflow

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


include { vep } from '../modules/vep.nf' params( params )
include { vepvcf2tsv } from '../modules/pysam' params( params )



workflow ensembl_vep {

    take:

    indexed_vcf_tuples
    indexed_ref_fasta
    vep_cache
    chrom_synonyms
    species


    main:

    // STEP 1 - Annotate the variants using VEP
    vep(
        indexed_vcf_tuples,
        indexed_ref_fasta,
        vep_cache,
        chrom_synonyms,
        species,
    )


    // STEP 2 - Convert the VEP VCF files to tab-separated values
    vep.out.annotated_variants \
        | map { vcf, tbi -> tuple( vcf.getBaseName(2), tuple( vcf, tbi ) ) } \
        | vepvcf2tsv


    emit:

    vep_vcf = vep.out.annotated_variants

    all_variants_tsv = vepvcf2tsv.out.all_variants
    pass_variants_tsv = vepvcf2tsv.out.pass_variants
}
