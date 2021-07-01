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


nextflow.enable.dsl=2

import nextflow.config.ConfigParser

nextflow_config = file( "${baseDir}/nextflow.config" ).text
parsed_config = new ConfigParser().setIgnoreIncludes( true ).parse( nextflow_config )
defaults = parsed_config.params

check_params()


/*
 * Imports
 */

// functions
include { parse_readgroup_csv } from './functions/input_csv_parsers.nf' params( params )
include { parse_somatic_csv } from './functions/input_csv_parsers.nf' params( params )
include { parse_germline_csv } from './functions/input_csv_parsers.nf' params( params )

include { get_header } from './functions/input_csv_parsers.nf' params( params )

// modules
include { concat_ref_genomes } from './modules/biopython.nf' params( params )
include { bwa_index } from './modules/bwakit.nf' params( params )
include { multiqc } from './modules/multiqc.nf' params( params )
include { samtools_faidx } from './modules/samtools.nf' params( params )

include { gunzip as gunzip_human_somatic_vcf } from './modules/gzip.nf' params( params )
include { gunzip as gunzip_mouse_somatic_vcf } from './modules/gzip.nf' params( params )

include { mosdepth as mosdepth_human } from './modules/mosdepth.nf' params( params )
include { mosdepth as mosdepth_mouse } from './modules/mosdepth.nf' params( params )

include { qualimap as qualimap_human } from './modules/qualimap.nf' params( params )
include { qualimap as qualimap_mouse } from './modules/qualimap.nf' params( params )

include { bcftools_subset_regions as subset_human_germline_variants } from './modules/bcftools.nf' params( params )
include { bcftools_subset_regions as subset_human_somatic_variants } from './modules/bcftools.nf' params( params )
include { bcftools_subset_regions as subset_mouse_germline_variants } from './modules/bcftools.nf' params( params )
include { bcftools_subset_regions as subset_mouse_somatic_variants } from './modules/bcftools.nf' params( params )

include { unpack_vep_cache as unpack_human_vep_cache } from './modules/vep.nf' params( params )
include { unpack_vep_cache as unpack_mouse_vep_cache } from './modules/vep.nf' params( params )

include { vcf2maf as human_vcf2maf } from './modules/vcf2maf.nf' params( params )
include { vcf2maf as mouse_vcf2maf } from './modules/vcf2maf.nf' params( params )

// workflows
include { dna_alignment } from './workflows/alignment.nf' params( params )

include { germline_variant_calling } from './workflows/variant_calling.nf' params( params )
include { somatic_variant_calling } from './workflows/variant_calling.nf' params( params )

include { mpileup as germline_mpileup } from './workflows/mpileup.nf' params( params )
include { mpileup as somatic_mpileup } from './workflows/mpileup.nf' params( params )

include { ensembl_vep as human_germline_annotation } from './workflows/annotation.nf' params( params )
include { ensembl_vep as human_somatic_annotation } from './workflows/annotation.nf' params( params )
include { ensembl_vep as mouse_germline_annotation } from './workflows/annotation.nf' params( params )
include { ensembl_vep as mouse_somatic_annotation } from './workflows/annotation.nf' params( params )


/*
 * Params
 */

// Cutadapt adapter files
params.r1_adapter_file = params.adapters in params.adapter_files
    ? params.adapter_files[ params.adapters ].r1_adapters
    : "${baseDir}/resource-adapters/null-1.fa"
params.r2_adapter_file = params.adapters in params.adapter_files
    ? params.adapter_files[ params.adapters ].r2_adapters
    : "${baseDir}/resource-adapters/null-2.fa"

// Reference genome files
params.human_ref_fasta = params.human_genome in params.human_genomes
    ? params.human_genomes[ params.human_genome ].ref_fasta
    : null
params.human_vep_cache = params.human_genome in params.human_genomes
    ? params.human_genomes[ params.human_genome ].vep_cache
    : null
params.mouse_ref_fasta = params.mouse_genome in params.mouse_genomes
    ? params.mouse_genomes[ params.mouse_genome ].ref_fasta
    : null
params.mouse_vep_cache = params.mouse_genome in params.mouse_genomes
    ? params.mouse_genomes[ params.mouse_genome ].vep_cache
    : null

// Mosdepth BED file
params.human_mosdepth_bed_file = "${baseDir}/assets/null"
params.mouse_mosdepth_bed_file = "${baseDir}/assets/null"

// QualiMap feature files
params.human_qualimap_feature_file = "${baseDir}/assets/null"
params.mouse_qualimap_feature_file = "${baseDir}/assets/null"

// Manta and Strelka call regions
params.call_regions = "${baseDir}/assets/null"

params.manta_call_regions = params.call_regions
params.strelka_call_regions = params.call_regions


/*
 * Workflow
 */

workflow {

    readgroup_inputs = parse_readgroup_csv( params.readgroup_csv )

    List<String> header = get_header( params.readgroup_csv )
    Boolean paired_end = header.contains( 'fastq1' ) && header.contains( 'fastq2' )
    Boolean single_end = header.contains( 'fastq' )

    cutadapt_adapter_files = [ params.r1_adapter_file, params.r2_adapter_file ]

    human_ref_inputs = [
        params.human_genome.replaceAll( /\./, '_' ),
        params.human_ref_fasta,
    ]
    mouse_ref_inputs = [
        params.mouse_genome.replaceAll( /\./, '_' ),
        params.mouse_ref_fasta,
    ]


    // STEP 1 - Create the combined human and mouse reference genome
    concat_ref_genomes( human_ref_inputs, mouse_ref_inputs )

    combined_ref_fasta = concat_ref_genomes.out.fasta

    human_regions_bed = concat_ref_genomes.out.human_regions_bed
    human_chrom_synonyms = concat_ref_genomes.out.human_chrom_synonyms

    mouse_regions_bed = concat_ref_genomes.out.mouse_regions_bed
    mouse_chrom_synonyms = concat_ref_genomes.out.mouse_chrom_synonyms


    // STEP 2 - Index the combined human and mouse reference genome
    samtools_faidx( combined_ref_fasta )
    bwa_index( combined_ref_fasta )

    indexed_ref_fasta = combined_ref_fasta \
        | concat( samtools_faidx.out ) \
        | collect


    // STEP 3 - Perform adapter trimming and alignment to the reference
    dna_alignment(
        readgroup_inputs,
        bwa_index.out,
        cutadapt_adapter_files,
    )


    // STEP 4 - Run Mosdepth
    mosdepth_inputs = ! params.skip_mosdepth
        ? dna_alignment.out.alignments.map { sample, indexed_bam -> indexed_bam }
        : Channel.empty()

    mosdepth_human( mosdepth_inputs, params.human_mosdepth_bed_file ?: human_regions_bed )
    mosdepth_mouse( mosdepth_inputs, params.mouse_mosdepth_bed_file ?: mouse_regions_bed )


    // STEP 5 - Run QualiMap
    qualimap_inputs = ! params.skip_qualimap
        ? dna_alignment.out.alignments.map { sample, indexed_bam -> indexed_bam.first() }
        : Channel.empty()

    qualimap_human( qualimap_inputs, params.human_qualimap_feature_file ?: human_regions_bed )
    qualimap_mouse( qualimap_inputs, params.mouse_qualimap_feature_file ?: mouse_regions_bed )


    // STEP 6 - Extract the indexed Ensembl VEP cache files
    if( params.germline_csv || params.somatic_csv ) {

        // human

        human_vep_cache = unpack_human_vep_cache( params.human_vep_cache )
        human_cache_info = human_vep_cache.cache_info.map {
            it.readLines().collectEntries() { it.split('\t') }
        }
        human_vep_species = human_cache_info.map { it['species'] }
        human_vep_assembly = human_cache_info.map { it['assembly'] }

        // mouse

        mouse_vep_cache = unpack_mouse_vep_cache( params.mouse_vep_cache )
        mouse_cache_info = mouse_vep_cache.cache_info.map {
            it.readLines().collectEntries() { it.split('\t') }
        }
        mouse_vep_species = mouse_cache_info.map { it['species'] }
        mouse_vep_assembly = mouse_cache_info.map { it['assembly'] }
    }


    // STEP 7 - Call and annotate germline variants
    if( params.germline_csv ) {

        parse_germline_csv( params.germline_csv ) \
            | map { analysis, samples ->
                tuple( groupKey(analysis, samples.size()), samples )
            } \
            | transpose() \
            | map { analysis, sample -> tuple( sample, analysis ) } \
            | combine( dna_alignment.out.alignments, by: 0 ) \
            | map { sample, analysis, indexed_bam -> tuple( analysis, indexed_bam ) } \
            | groupTuple() \
            | map { analysis, indexed_bam_files ->
                tuple( analysis.toString(), indexed_bam_files.flatten() )
            } \
            | set { germline_inputs }

        germline_variant_calling(
            germline_inputs,
            indexed_ref_fasta,
            params.strelka_call_regions,
        )

        germline_mpileup(
            germline_inputs,
            germline_variant_calling.out,
            indexed_ref_fasta,
        )

        // human

        subset_human_germline_variants(
            germline_mpileup.out,
            human_regions_bed,
        )

        human_germline_annotation(
            subset_human_germline_variants.out,
            indexed_ref_fasta,
            human_vep_cache.cache_dir,
            human_chrom_synonyms,
            human_vep_species,
            human_vep_assembly,
        )

        // mouse

        subset_mouse_germline_variants(
            germline_mpileup.out,
            mouse_regions_bed,
        )

        mouse_germline_annotation(
            subset_mouse_germline_variants.out,
            indexed_ref_fasta,
            mouse_vep_cache.cache_dir,
            mouse_chrom_synonyms,
            mouse_vep_species,
            mouse_vep_assembly,
        )
    }


    // STEP 8 - Call and annotate somatic variants
    if( params.somatic_csv ) {

        test_control_inputs = parse_somatic_csv( params.somatic_csv )

        test_control_inputs \
            | map { analysis, test, control ->
                tuple( test, tuple( analysis, test, control ) )
            } \
            | combine( dna_alignment.out.alignments, by: 0 ) \
            | map { sample, key_tuple, indexed_bam -> tuple( *key_tuple, indexed_bam ) } \
            | set { test_inputs }

        test_control_inputs \
            | map { analysis, test, control ->
                tuple( control, tuple( analysis, test, control ) )
            } \
            | combine( dna_alignment.out.alignments, by: 0 ) \
            | map { sample, key_tuple, indexed_bam -> tuple( *key_tuple, indexed_bam ) } \
            | set { control_inputs }

        test_control_inputs \
            | join( test_inputs, by: [0,1,2] ) \
            | join( control_inputs, by: [0,1,2] ) \
            | map { analysis, test, control, indexed_test_bam, indexed_control_bam ->
                tuple( analysis, indexed_test_bam, indexed_control_bam )
            } \
            | set { somatic_inputs }

        somatic_variant_calling(
            somatic_inputs,
            indexed_ref_fasta,
            params.manta_call_regions,
            params.strelka_call_regions,
            single_end,
        )

        somatic_mpileup(
            somatic_inputs,
            somatic_variant_calling.out,
            indexed_ref_fasta,
        )

        // human

        subset_human_somatic_variants(
            somatic_mpileup.out,
            human_regions_bed,
        )

        human_somatic_annotation(
            subset_human_somatic_variants.out,
            indexed_ref_fasta,
            human_vep_cache.cache_dir,
            human_chrom_synonyms,
            human_vep_species,
            human_vep_assembly,
        )

        human_somatic_annotation.out.vep_vcf \
            | map { vcf, tbi -> vcf } \
            | gunzip_human_somatic_vcf \
            | map { vcf -> tuple( vcf.getSimpleName(), vcf ) } \
            | join( test_control_inputs ) \
            | map { analysis, vcf, test, control -> tuple( vcf, test, control ) } \
            | set { human_vcf2maf_inputs }

        human_vcf2maf(
            human_vcf2maf_inputs,
            indexed_ref_fasta,
            human_vep_species,
            human_vep_assembly,
        )

        // mouse

        subset_mouse_somatic_variants(
            somatic_mpileup.out,
            mouse_regions_bed,
        )

        mouse_somatic_annotation(
            subset_mouse_somatic_variants.out,
            indexed_ref_fasta,
            mouse_vep_cache.cache_dir,
            mouse_chrom_synonyms,
            mouse_vep_species,
            mouse_vep_assembly,
        )

        mouse_somatic_annotation.out.vep_vcf \
            | map { vcf, tbi -> vcf } \
            | gunzip_mouse_somatic_vcf \
            | map { vcf -> tuple( vcf.getSimpleName(), vcf ) } \
            | join( test_control_inputs ) \
            | map { analysis, vcf, test, control -> tuple( vcf, test, control ) } \
            | set { mouse_vcf2maf_inputs }

        mouse_vcf2maf(
            mouse_vcf2maf_inputs,
            indexed_ref_fasta,
            mouse_vep_species,
            mouse_vep_assembly,
        )
    }


    // STEP 9 - Create a MultiQC report
    human_logs = Channel.empty() \
        | mix( mosdepth_human.out.dists ) \
        | mix( mosdepth_human.out.summary ) \
        | mix( qualimap_human.out ) \
        | collect

    mouse_logs = Channel.empty() \
        | mix( mosdepth_mouse.out.dists ) \
        | mix( mosdepth_mouse.out.summary ) \
        | mix( qualimap_mouse.out ) \
        | collect

    logs = Channel.empty() \
        | mix( dna_alignment.out.logs ) \
        | collect

    multiqc( logs, human_logs, mouse_logs, params.multiqc_config )
}


workflow.onComplete {

    log.info "Workflow completed at: ${workflow.complete}"
    log.info "Time taken: ${workflow.duration}"
    log.info "Execution status: ${workflow.success ? 'success' : 'failed'}"
    log.info "Output directory: ${params.publish_dir}"
}


workflow.onError {

    log.info "Execution halted: ${workflow.errorMessage}"
}


/*
 * Functions
 */

def check_params() {

    if( params.help ) {
        usage()
        exit 0
    }

    if( params.version ) {
        log.info( workflow.manifest.version )
        exit 0
    }

    if( !params.readgroup_csv ) {
        log.info "Readgroup CSV not specified. Please using the `--readgroup_csv` parameter"
        exit 1
    }
}

def usage() {

    log.info"""
    Usage:
        nextflow run -profile <profile> -revision <revision> ampatchlab/nf-pdx [options]


    Nextflow execution options:

        -profile STR
            Nextflow configuration profile to use. Available profiles include:
            'awsbatch', 'conda', 'docker' and 'singularity'

        -revision STR
            Git branch/tag (version) of this workflow to use

        -work-dir DIR
            Directory where intermediate result files are stored

        -help
            Show additional execution options and exit


    Workflow input params:

        --readgroup_csv FILE
            Comma-separated list of sample and readgroup inputs
            Please see: https://github.com/ampatchlab/nf-pdx#readgroup-inputs

        --germline_csv FILE
            Comma-separated list of analysis identifiers and sample inputs
            Please see: https://github.com/ampatchlab/nf-pdx#germline-inputs

        --somatic_csv FILE
            Comma-separated list of analysis identifiers and test and control sample inputs
            Please see: https://github.com/ampatchlab/nf-pdx#somatic-inputs


    Reference genome params:

        --human_genome STR
            Human genome name [Default: ${defaults.human_genome}]

        --human_ref_fasta FILE
            Override the human reference FASTA file with FILE [Default: ${defaults.human_ref_fasta ?: null}]

        --human_vep_cache FILE
            Override the human VEP cache with FILE [Default: ${defaults.human_vep_cache ?: null}]

        --mouse_genome STR
            Human genome name [Default: ${defaults.mouse_genome}]

        --mouse_ref_fasta FILE
            Override the mouse reference FASTA file with FILE [Default: ${defaults.mouse_ref_fasta ?: null}]

        --mouse_vep_cache FILE
            Override the mouse VEP cache with FILE [Default: ${defaults.mouse_vep_cache ?: null}]


    Adapter trimming params:

        --adapters STR
            The adapters to trim [Either: ${defaults.adapter_files.keySet().join(", ")}; Default: ${defaults.adapters}]

        --r1_adapter_file FILE
            Override the R1 adapter file with FILE [Default: ${defaults.r1_adapter_file ?: null}]

        --r2_adapter_file FILE
            Override the R2 adapter file with FILE [Default: ${defaults.r2_adapter_file ?: null}]


    Mosdepth params:

        --human_mosdepth_bed_file FILE
            Override the Mosdepth BED file of human regions [Default: ${defaults.human_mosdepth_bed_file ?: null}]

        --mouse_mosdepth_bed_file FILE
            Override the Mosdepth BED file of mouse regions [Default: ${defaults.mouse_mosdepth_bed_file ?: null}]

        --skip_mosdepth
            Skips Mosdepth process execution


    Qualimap params:

        --human_qualimap_feature_file FILE
            Override the QualiMap feature file of human regions in GFF/GTF or BED format [Default: ${defaults.human_qualimap_feature_file ?: null}]

        --mouse_qualimap_feature_file FILE
            Override the QualiMap feature file of mouse regions in GFF/GTF or BED format [Default: ${defaults.mouse_qualimap_feature_file ?: null}]

        --skip_qualimap
            Skips QualiMap process execution


    Strelka and Manta params:

        --call_regions FILE
            Restrict variant calling to the regions in BED FILE [Default: ${defaults.call_regions ?: null}]

        --exome
            Provide appropriate settings for WES and other regional enrichment analyses


    Ensembl-VEP params:

        --vep_cache_type STR
            Apply when using a merged or refseq VEP cache [Either: 'merged', 'refseq'; Default: ${defaults.vep_cache_type ?: null}]


    Output params:

        --publish_dir DIR
            Path where the results will be published [Default: ${defaults.publish_dir}]

        --publish_mode STR
            The mode used to publish files to the target directory [Default: ${defaults.publish_mode}]


    Standard params:

        --help
            Show this message and exit

        --version
            Show the pipeline version and exit
    """.stripIndent()
}
