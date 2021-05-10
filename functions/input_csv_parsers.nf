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


class DictReader implements Iterable<Map<String, String>> {

    private BufferedReader reader
    private CsvSplitter splitter

    private final String filePath

    DictReader(String filePath ) {
        this.reader = new BufferedReader( new FileReader( filePath ) )
        this.splitter = new CsvSplitter().options( header: true )

        splitter.parseHeader(reader)
    }

    public List<String> getHeader() {
        return splitter.columnsHeader
    }

    @Override
    public Iterator<Map<String, String>> iterator() {

        return new Iterator<Map<String, String>>() {

            private Map<String, String> nextLine

            @Override
            public boolean hasNext() {
               if( nextLine == null ) {
                   nextLine = splitter.fetchRecord( reader )
               }

               return nextLine != null
            }

            @Override
            public Map<String, String> next() {
                Map<String, String> nextline = nextLine
                nextLine = null

                return nextline
            }
        }
    }
}


def get_header( csv ) {

    DictReader reader = new DictReader( csv )

    return reader.getHeader()
}


def parse_readgroup_csv( csv ) {

    DictReader reader = new DictReader( csv )

    List<String> header = reader.getHeader()

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
        error "${csv} - Invalid column header configuration: ${header}"
    }

    def readgroups = []
    def fastq_filenames = []

    def inputs = []

    for( Map<String,String> row : reader ) {

        if( !header_cols.every() ) {
            error "${csv} - Missing one of ${header_cols} for row: ${row}"
        }

        def sample_fields = sample_cols.collect { row[ it ].replaceAll( /\./, '_' ) }

        def sample = sample_fields.first()
        def readgroup = sample_fields.last()

        if( readgroup in readgroups ) {
            error "${csv} - Cannot use this sample/readgroup label more than once: ${readgroup}"
        }
        readgroups.add( readgroup )

        List<Path> fastqs = fastq_cols.collect { file( row[ it ], glob: false ) }

        for( Path fastq : fastqs ) {

            def extensions = [ '.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fastq.bz2', '.fq.bz2' ]

            if( !extensions.find { fastq.name.endsWith( it ) } ) {
                error "${csv} - FASTQ file has an unsupported extension: ${fastq.name}"
            }
            if( fastq.name in fastq_filenames ) {
                error "${csv} - FASTQ file must not be used more than once: ${fastq.name}"
            }
            fastq_filenames.add( fastq.name )
        }

        inputs.add( tuple( sample, readgroup, fastqs ) )
    }

    return Channel.from( inputs )
}


def parse_germline_csv( csv ) {

    DictReader reader = new DictReader( csv )

    List<String> header = reader.getHeader()

    def required_cols = [ 'analysis', 'samples' ]

    if( !required_cols.every { header.count(it) == 1 } ) {
        error "${csv} - Invalid column header configuration: ${header}"
    }

    def analyses = []

    def inputs = []

    for( Map<String,String> row : reader ) {

        if( !required_cols.every() ) {
            error "${csv} - Missing one of ${required_cols} for row: ${row}"
        }

        def (analysis, samples) = required_cols.collect { row[ it ].replaceAll( /\./, '_' ) }

        if( analysis in analyses ) {
            error "${csv} - Cannot use this analysis more than once: ${analysis}"
        }
        analyses.add( analysis )

        List<String> sample_list = samples.tokenize( '|' ).unique()

        inputs.add( tuple( analysis, sample_list ) )
    }

    return Channel.from( inputs )
}


def parse_somatic_csv( csv ) {

    DictReader reader = new DictReader( csv )

    List<String> header = reader.getHeader()

    def required_cols = [ 'analysis', 'test', 'control' ]

    if( !required_cols.every { header.count(it) == 1 } ) {
        error "${csv} - Invalid column header configuration: ${header}"
    }

    def analyses = []

    def inputs = []

    for( Map<String,String> row : reader ) {

        if( !required_cols.every() ) {
            error "${csv} - Missing one of ${required_cols} for row: ${row}"
        }

        def (analysis, test, control) = required_cols.collect { row[ it ].replaceAll( /\./, '_' ) }

        if( analysis in analyses ) {
            error "${csv} - Cannot use this analysis more than once: ${analysis}"
        }
        if( test == control ) {
            error "${csv} - Cannot use the same test/control sample for analysis: ${analysis}"
        }
        analyses.add( analysis )

        inputs.add( tuple( analysis, test, control ) )
    }

    return Channel.from( inputs )
}
