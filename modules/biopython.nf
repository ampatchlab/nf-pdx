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
    path "combined_${human_prefix}_${mouse_prefix}.fa"

    script:
    """
    #!/usr/bin/env python
    import gzip

    from contextlib import ExitStack
    from itertools import chain

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

    def prefix_chroms(filename, prefix):
        for rec in SeqIO.parse(filename, 'fasta'):
            rec.id = '_'.join([prefix, rec.id])
            yield rec

    with ExitStack() as stack:

        out_fn = 'combined_${human_prefix}_${mouse_prefix}.fa'
        out_fh = stack.enter_context(openfile(out_fn, 'w'))

        human_ref_fh = stack.enter_context(openfile('${human_ref_fasta}'))
        mouse_ref_fh = stack.enter_context(openfile('${mouse_ref_fasta}'))

        human_recs = prefix_chroms(human_ref_fh, prefix='${human_prefix}')
        mouse_recs = prefix_chroms(mouse_ref_fh, prefix='${mouse_prefix}')

        SeqIO.write(chain(human_recs, mouse_recs), out_fh, 'fasta')
    """
}
