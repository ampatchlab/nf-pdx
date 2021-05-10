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
        saveAs: { fn ->
            if( fn.endsWith('.bed') ) {
                return "regions_bed/${fn}"
            }
            if( fn.endsWith('.tsv') ) {
                return "chrom_synonyms/${fn}"
            }
            return fn
        }
    )

    input:
    tuple val(human_prefix), path(human_ref_fasta)
    tuple val(mouse_prefix), path(mouse_ref_fasta)

    output:
    path "combined_${human_prefix}_${mouse_prefix}.fa", emit: fasta
    path "${human_prefix}.bed", emit: human_regions_bed
    path "${human_prefix}.tsv", emit: human_chrom_synonyms
    path "${mouse_prefix}.bed", emit: mouse_regions_bed
    path "${mouse_prefix}.tsv", emit: mouse_chrom_synonyms

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

    def seqio_writer(output_handle, format_string):
        while True:
            rec = (yield)
            SeqIO.write(rec, output_handle, format_string)

    human_chroms = []
    mouse_chroms = []

    human_prefix = '${human_prefix}'
    mouse_prefix = '${mouse_prefix}'

    human_fasta = '${human_ref_fasta}'
    mouse_fasta = '${mouse_ref_fasta}'

    writer_kwargs = { 'delimiter': '\\t', 'lineterminator': '\\n' }

    with ExitStack() as stack:

        human_ref_fh = stack.enter_context(openfile(human_fasta))
        mouse_ref_fh = stack.enter_context(openfile(mouse_fasta))

        out_fn = '_'.join(['combined', human_prefix, mouse_prefix]) + '.fa'
        out_fh = stack.enter_context(openfile(out_fn, 'w'))

        seq_writer = seqio_writer(out_fh, 'fasta')
        next(seq_writer)

        for rec in SeqIO.parse(human_ref_fh, 'fasta'):
            human_chroms.append((rec.id, len(rec.seq)))
            rec.id = '_'.join([human_prefix, rec.id])
            seq_writer.send(rec)

        for rec in SeqIO.parse(mouse_ref_fh, 'fasta'):
            mouse_chroms.append((rec.id, len(rec.seq)))
            rec.id = '_'.join([mouse_prefix, rec.id])
            seq_writer.send(rec)

        seq_writer.close()

    # BED regions
    with openfile(human_prefix + '.bed', 'w') as bed:
        writer = csv.writer(bed, **writer_kwargs)
        for rec_id, seq_length in human_chroms:
            writer.writerow(['_'.join([human_prefix, rec_id]), 0, seq_length])
    with openfile(mouse_prefix + '.bed', 'w') as bed:
        writer = csv.writer(bed, **writer_kwargs)
        for rec_id, seq_length in mouse_chroms:
            writer.writerow(['_'.join([mouse_prefix, rec_id]), 0, seq_length])

    # VEP synonyms
    with openfile(human_prefix + '.tsv', 'w') as synonyms:
        writer = csv.writer(synonyms, **writer_kwargs)
        for rec_id, seq_length in human_chroms:
            writer.writerow(['_'.join([human_prefix, rec_id]), rec_id])
    with openfile(mouse_prefix + '.tsv', 'w') as synonyms:
        writer = csv.writer(synonyms, **writer_kwargs)
        for rec_id, seq_length in mouse_chroms:
            writer.writerow(['_'.join([mouse_prefix, rec_id]), rec_id])
    """
}
