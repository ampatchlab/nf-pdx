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

params.publish_vepvcf2tsv = false

params.vcf_info_field = 'CSQ'


/*
 * Processes
 */

process vepvcf2tsv {

    tag { sample }

    label 'pysam'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_vepvcf2tsv,
        mode: params.publish_mode,
    )

    input:
    tuple val(sample), path(indexed_vcf)

    output:
    path "${sample}.tsv.gz{,.tbi}", emit: all_variants
    path "${sample}.pass.tsv.gz{,.tbi}", emit: pass_variants

    """
    #!/usr/bin/env python
    import csv
    import re
    import shlex

    from contextlib import ExitStack

    import pysam

    all_tsv = '${sample}.tsv'
    pass_tsv = '${sample}.pass.tsv'

    tabix = pysam.TabixFile('${indexed_vcf.first()}', parser=pysam.asVCF())

    *comment_lines, header_line = tabix.header

    info_cols = []
    format_cols = []
    csq_cols = []

    for comment_line in comment_lines:
        meta_line = re.match('^##(?P<key>.*)=<(?P<value>.*)>', comment_line)

        if meta_line:
            lexer = shlex.shlex(meta_line.group('value'), posix=True)

            lexer.whitespace_split = True
            lexer.whitespace = ','

            meta = dict(r.split('=', 1) for r in lexer)

            if meta_line.group('key') == 'INFO':
                if meta['ID'] == '${params.vcf_info_field}':
                    csq_cols = re.search('Format: (.*)', meta['Description']).group(1).split('|')
                else:
                    info_cols.append('/'.join(['INFO', meta['ID']]))

            if meta_line.group('key') == 'FORMAT':
                format_cols.append(meta['ID'])

    headers = re.match('^#(.*)', header_line).group(1).split('\\t')
    mandatory_header_cols, remaining_header_cols = headers[:8], headers[8:]

    if remaining_header_cols:
        # first column must be 'FORMAT' if genotype data is present
        assert remaining_header_cols.pop(0) == 'FORMAT'
        # sample names must not match any of the eight fixed mandatory columns
        assert not any(s in mandatory_header_cols for s in remaining_header_cols)

    sample_cols = ['/'.join([s, key]) for s in remaining_header_cols for key in format_cols]

    fieldnames = mandatory_header_cols[:-1] + info_cols + sample_cols + csq_cols

    with ExitStack() as stack:

        all_tsv_fh = stack.enter_context(open(all_tsv, 'w'))
        all_tsv_writer = csv.DictWriter(all_tsv_fh, fieldnames=fieldnames, delimiter='\\t')
        all_tsv_writer.writeheader()

        pass_tsv_fh = stack.enter_context(open(pass_tsv, 'w'))
        pass_tsv_writer = csv.DictWriter(pass_tsv_fh, fieldnames=fieldnames, delimiter='\\t')
        pass_tsv_writer.writeheader()

        for row in (dict(zip(headers, r)) for r in tabix.fetch()):
            result_dict = { c: row[c] for c in mandatory_header_cols[:-1] }
            csq_list = []
            for rec in row['INFO'].split(';'):
                key, sep, value = rec.partition('=')
                if key == '${params.vcf_info_field}':
                    for r in value.split(','):
                        csq_list.append({ k: v.replace('&', ',') for k, v in zip(csq_cols, r.split('|')) })
                else:
                    result_dict['/'.join(['INFO', key])] = value if sep else True

            for sample in remaining_header_cols:
                key_list = ['/'.join([sample, x]) for x in row['FORMAT'].split(':')]
                for k, v in zip(key_list, row[sample].split(':')):
                    result_dict[k] = v

            for csq_dict in csq_list or [{}]:
                all_tsv_writer.writerow({**result_dict, **csq_dict})
                if row['FILTER'] == 'PASS':
                    pass_tsv_writer.writerow({**result_dict, **csq_dict})

    index_kwargs = { 'seq_col': 0, 'start_col': 1, 'end_col': 1, 'line_skip': 1 }

    pysam.tabix_index(all_tsv, **index_kwargs)
    pysam.tabix_index(pass_tsv, **index_kwargs)
    """
}
