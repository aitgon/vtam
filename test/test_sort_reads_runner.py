# -*- coding: utf-8 -*-
import os
import pathlib
from unittest import TestCase

import pandas

from vtam.utils.PathManager import PathManager
from vtam.utils.SortReadsRunner import SortReadsRunner


class TestSortReadsRunner(TestCase):

    def test_sort_reads_runner(self):

        #######################################################################
        #
        # Fasta path
        #
        #######################################################################

        reads_fasta_str = """>M00842:118:000000000-ABGKE:1:1101:18229:3444 1:N:0:5
TCGATCACGATGTTCCACTAATCACAAGGATATTGGTACCCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGAT
CATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATT
GTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATTGGAGGATTTGGAAATTGATTAGTTGC
TGTAGATCGACA
>M00842:118:000000000-ABGKE:1:1101:19104:2756 1:N:0:5
TGATCGATGATCAGTCCACTAATCACAAGGATATTGGTACCCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGA
TCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAAT
TGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATTGGAGGATTTGGAAATTGATTAGTTG
CTGTAGATCGACA
>M00842:118:000000000-ABGKE:1:1101:15187:2443 1:N:0:5
TGATCGATGATCAGTCCACTAATCACAAGGATATTGGTACCCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGA
TCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAAT
TGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATTGGAGGATTTGGTAATTGATTAGTAG
CTGTAGATCGACA"""
        this_tempdir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(this_tempdir).mkdir(exist_ok=True)

        fasta_file_name = 'reads.fasta'
        fasta_path = os.path.join(this_tempdir, fasta_file_name)
        outdir = os.path.join(PathManager.get_module_test_path(), 'output')
        pathlib.Path(outdir).mkdir(exist_ok=True)

        with open(fasta_path, 'w') as fout:
            fout.write(reads_fasta_str)

        alignement_parameters = {'min_id': 0.8, 'minseqlength': 32, 'overhang': 0}

        fasta_information_df = pandas.DataFrame(
            {
                'Run': [1, 1],
                'Marker': [1, 1],
                'Biosample': [2, 5],
                'Replicate': [2, 2],
                'TagFwd': ['tcgatcacgatgt', 'tgatcgatgatcag'],
                'PrimerFwd': ['TCCACTAATCACAARGATATTGGTAC', 'TCCACTAATCACAARGATATTGGTAC'],
                'TagRev': ['tgtcgatctacagc', 'tgtcgatctacagc'],
                'PrimerRev': ['WACTAATCAATTWCCAAATCCTCC', 'WACTAATCAATTWCCAAATCCTCC'],
                'Fasta': [fasta_file_name, fasta_file_name],
             }
        )


        #######################################################################
        #
        # Sort Reads
        #
        #######################################################################

        # Create SortReadsRunner
        sort_reads_tsv_path = os.path.join(this_tempdir, 'sort_reads.tsv')
        sort_reads_runner = SortReadsRunner(fasta_path=fasta_path, alignement_parameters=alignement_parameters,
                                                fasta_information_df=fasta_information_df, outdir=outdir)

        sort_reads_runner.run()

        fasta_info_str_bak = """Run	Marker	Biosample	Replicate	Fasta
1	1	2	2	reads_000.fasta
1	1	5	2	reads_001.fasta
"""

        with open(os.path.join(outdir, 'fasta_info.tsv')) as fin:
            fasta_info_str = fin.read()
        assert fasta_info_str_bak == fasta_info_str

        fasta_str_bak = """>M00842:118:000000000-ABGKE:1:1101:18229:3444
CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATT
"""

        with open(os.path.join(outdir, 'reads_000.fasta')) as fin:
            fasta_str = fin.read()
        assert fasta_str_bak == fasta_str
