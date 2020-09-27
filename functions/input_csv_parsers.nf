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

import nextflow.splitter.CsvSplitter


def parse_input_csv( csv ) {

    def splitter = new CsvSplitter().options( header:true )
    def reader = new BufferedReader( new FileReader( csv ) )

    splitter.parseHeader(reader)

    List<String> header = splitter.columnsHeader

    def sample_cols = [ 'sample' ] + ( 'readgroup' in header ? [ 'readgroup' ] : [] )

    def pe_fastq_cols = [ 'fastq1', 'fastq2' ]
    def se_fastq_cols = [ 'fastq' ]

    def pe_header_cols = sample_cols + pe_fastq_cols
    def se_header_cols = sample_cols + se_fastq_cols

    Boolean every_pe_header_col = pe_header_cols.every { header.count(it) == 1 }
    Boolean every_se_header_col = se_header_cols.every { header.count(it) == 1 }

    Boolean any_pe_fastq_cols = pe_fastq_cols.any { header.contains( it ) }
    Boolean any_se_fastq_cols = se_fastq_cols.any { header.contains( it ) }

    def header_cols = []
    def fastq_cols = []

    if( every_pe_header_col && !any_se_fastq_cols ) {
        header_cols = pe_header_cols
        fastq_cols = pe_fastq_cols
    }
    else if( every_se_header_col && !any_pe_fastq_cols ) {
        header_cols = se_header_cols
        fastq_cols = se_fastq_cols
    }
    else {
        error "${csv.name} - Invalid column header configuration: ${header}"
    }

    def readgroups = []
    def fastq_filenames = []

    def inputs = []

    Map<String,String> row

    while( row = splitter.fetchRecord( reader ) ) {

        if( !header_cols.every() ) {
            error "${csv.name} - Missing one of ${header_cols} for row: ${row}"
        }

        def sample_fields = sample_cols.collect { row[ it ].replaceAll( /\./, '_' ) }

        def sample = sample_fields.first()
        def readgroup = sample_fields.last()

        if( readgroup in readgroups ) {
            error "${csv.name} - Must not use sample/readgroup labels more than once: ${readgroup}"
        }
        readgroups.add( readgroup )

        List<Path> fastqs = fastq_cols.collect { file( row[ it ], glob: false ) }

        for( Path fastq : fastqs ) {

            def extensions = [ '.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fastq.bz2', '.fq.bz2' ]

            if( !extensions.find { fastq.name.endsWith( it ) } ) {
                error "${csv.name} - FASTQ file has an unsupported extension: ${fastq.name}"
            }
            if( fastq.name in fastq_filenames ) {
                error "${csv.name} - FASTQ file must not be used more than once: ${fastq.name}"
            }
            fastq_filenames.add( fastq.name )
        }

        inputs.add( tuple( sample, readgroup, fastqs ) )
    }

    return Channel.from( inputs )
}
