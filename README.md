# ampatchlab/nf-pdx

[![Build Status](https://codebuild.ap-southeast-2.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoibVA2TXNnVys3QkQzWHdzOENCbE4vVThmL0dmQWR4TlpXRlM3eE1DaWZrK2dYbkdkQ1A1bVJpWlFCSmZFbmpySU96SVEzU2VjTFBaUkFJMUpPN2JZZmFjPSIsIml2UGFyYW1ldGVyU3BlYyI6IlRLblJJNmpuU01UMnF5TWoiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=master)](https://ap-southeast-2.console.aws.amazon.com/codesuite/codebuild/projects/nf-pdx/history)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Patient-derived xenograft Nextflow pipeline

## Usage

```
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
        Human genome name [Default: GRCh38]

    --human_ref_fasta FILE
        Override the human reference FASTA file with FILE [Default: null]

    --human_vep_cache FILE
        Override the human VEP cache with FILE [Default: null]

    --mouse_genome STR
        Human genome name [Default: GRCm38]

    --mouse_ref_fasta FILE
        Override the mouse reference FASTA file with FILE [Default: null]

    --mouse_vep_cache FILE
        Override the mouse VEP cache with FILE [Default: null]


Adapter trimming params:

    --adapters STR
        The adapters to trim [Either: TruSeq, NexteraTransposase, BGISeq; Default: null]

    --r1_adapter_file FILE
        Override the R1 adapter file with FILE [Default: null]

    --r2_adapter_file FILE
        Override the R2 adapter file with FILE [Default: null]


Mosdepth params:

    --human_mosdepth_bed_file FILE
        Override the Mosdepth BED file of human regions [Default: null]

    --mouse_mosdepth_bed_file FILE
        Override the Mosdepth BED file of mouse regions [Default: null]

    --skip_mosdepth
        Skips Mosdepth process execution


Qualimap params:

    --human_qualimap_feature_file FILE
        Override the QualiMap feature file of human regions in GFF/GTF or BED format [Default: null]

    --mouse_qualimap_feature_file FILE
        Override the QualiMap feature file of mouse regions in GFF/GTF or BED format [Default: null]

    --skip_qualimap
        Skips QualiMap process execution


Strelka and Manta params:

    --call_regions FILE
        Restrict variant calling to the regions in BED FILE [Default: null]

    --exome
        Provide appropriate settings for WES and other regional enrichment analyses


Ensembl-VEP params:

    --vep_cache_type STR
        Apply when using a merged or refseq VEP cache [Either: 'merged', 'refseq'; Default: null]


Output params:

    --publish_dir DIR
        Path where the results will be published [Default: ./results]

    --publish_mode STR
        The mode used to publish files to the target directory [Default: copy]


Standard params:

    --help
        Show this message and exit

    --version
        Show the pipeline version and exit
```

## Readgroup inputs

Readgroup inputs must be supplied using the `--readgroup_csv` parameter. This parameter
expects a simple file containing comma-separated values and an appropriate header line.

For paired-end data, the readgroup CSV must have the following columns:

 * sample: Unique sample name or ID (required)
 * readgroup: Unique readgroup name or ID (optional)
 * fastq1: Absolute path of the 'R1' FASTQ file (required)
 * fastq2: Absolute path of the 'R2' FASTQ file (required)

For single-end data, the readgroup CSV must have the following columns:

 * sample: Unique sample name or ID (required)
 * readgroup: Unique readgroup name or ID (optional)
 * fastq: Absolute path of the FASTQ file (required)

If a particular sample has multiple FASTQ files (or pairs of FASTQ files), then these may
be specified on additional lines with a unique readgroup identifier. All readgroups belonging
to a particular sample will be aligned and merged.


## Germline inputs

To enable germline variant calling, another CSV file must be supplied using the
`--germline_csv` parameter. This parameter expects a simple CSV file and a header line.

The germline CSV file must have the following columns:

 * analysis: Unique analysis label or ID (required)
 * samples: One or more unique sample names or IDs (required)

Note that if multiple samples are specified, these MUST be pipe-delimited. If they are comma-
separated, only the first sample will be analyzed. Joint germline analysis is limited to
trios or family scale (tens of samples) and is not intended to support cohort, case/control
or population analysis.


## Somatic inputs

To enable somatic variant calling, another CSV file must be supplied using the
`--somatic_csv` parameter. This parameter expects a simple CSV file and a header line.

The somatic CSV file must have the following columns:

 * analysis: Unique analysis label or ID (required)
 * test: Unique sample name or ID of the 'test' sample (required)
 * control: Unique sample name or ID of the 'control' sample (required)
