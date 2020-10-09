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


/*
 * Params
 */

params.publish_dir = './results'
params.publish_everything = false
params.publish_mode = 'copy'

params.publish_concat_ref_genomes = false


/*
 * Processes
 */

process concat_ref_genomes {

    tag { "combined_${human_prefix}_${mouse_prefix}.fa" }

    label 'biopython'

    publishDir(
        path: "${params.publish_dir}/${task.process.replaceAll(':', '/')}",
        enabled: params.publish_everything || params.publish_concat_ref_genomes,
        mode: params.publish_mode,
    )

    input:
    tuple val(human_prefix), path(human_ref_fasta)
    tuple val(mouse_prefix), path(mouse_ref_fasta)

    output:
    path "combined_${human_prefix}_${mouse_prefix}.fa", emit: fasta
    path "${human_prefix}.tsv", emit: human_chrom_synonyms
    path "${human_prefix}.txt", emit: human_chrom_list
    path "${mouse_prefix}.tsv", emit: mouse_chrom_synonyms
    path "${mouse_prefix}.txt", emit: mouse_chrom_list

    """
    #!/usr/bin/env python
    import csv
    import gzip

    from contextlib import ExitStack
    from Bio import SeqIO

    def openfile(filename, mode='r'):
        if filename.endswith('.gz'):
            if mode == 'r':
                return gzip.open(filename, 'rt')
            if mode == 'w':
                return gzip.open(filename, 'wt')
            return gzip.open(filename, mode)
        else:
            return open(filename, mode)

    def seq_writer(output_handle, format_string):
        while True:
            rec = (yield)
            SeqIO.write(rec, output_handle, format_string)

    human_prefix = '${human_prefix}'
    mouse_prefix = '${mouse_prefix}'

    human_fasta = '${human_ref_fasta}'
    mouse_fasta = '${mouse_ref_fasta}'

    with ExitStack() as stack:

        human_ref_fh = stack.enter_context(openfile(human_fasta))
        mouse_ref_fh = stack.enter_context(openfile(mouse_fasta))

        out_fn = '_'.join(['combined', human_prefix, mouse_prefix]) + '.fa'
        out_fh = stack.enter_context(openfile(out_fn, 'w'))

        writer = seq_writer(out_fh, 'fasta')
        next(writer)

        human_chrom_list = stack.enter_context(openfile(human_prefix + '.txt', 'w'))
        mouse_chrom_list = stack.enter_context(openfile(mouse_prefix + '.txt', 'w'))

        human_chrom_syns = stack.enter_context(openfile(human_prefix + '.tsv', 'w'))
        mouse_chrom_syns = stack.enter_context(openfile(mouse_prefix + '.tsv', 'w'))

        csv_kwargs = { 'delimiter': '\\t', 'lineterminator': '\\n' }

        human_syns_csv = csv.writer(human_chrom_syns, **csv_kwargs)
        mouse_syns_csv = csv.writer(mouse_chrom_syns, **csv_kwargs)

        for rec in SeqIO.parse(human_ref_fh, 'fasta'):
            prefixed_id = '_'.join([human_prefix, rec.id])
            human_chrom_list.write(prefixed_id + '\\n')
            human_syns_csv.writerow([prefixed_id, rec.id])
            rec.id = prefixed_id
            writer.send(rec)

        for rec in SeqIO.parse(mouse_ref_fh, 'fasta'):
            prefixed_id = '_'.join([mouse_prefix, rec.id])
            mouse_chrom_list.write(prefixed_id + '\\n')
            mouse_syns_csv.writerow([prefixed_id, rec.id])
            rec.id = prefixed_id
            writer.send(rec)

        writer.close()
    """
}
