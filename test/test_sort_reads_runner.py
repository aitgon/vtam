# -*- coding: utf-8 -*-
import os
import pathlib
import sqlite3
from unittest import TestCase

import pandas

from vtam.utils.PathManager import PathManager
from vtam.utils.ReadTrimmer import ReadTrimmer
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
        fasta_id = 1

        with open(fasta_path, 'w') as fout:
            fout.write(reads_fasta_str)

        alignement_parameters = {'min_id': 0.8, 'minseqlength': 32}

        fasta_information_df = pandas.DataFrame(
            {
                'run_id': [1, 1],
                'marker_id': [1, 1],
                'biosample_id': [2, 5],
                'replicate_id': [2, 2],
                'tag_fwd_sequence': ['tcgatcacgatgt', 'tgatcgatgatcag'],
                'primer_fwd_sequence': ['TCCACTAATCACAARGATATTGGTAC', 'TCCACTAATCACAARGATATTGGTAC'],
                'tag_rev_sequence': ['tgtcgatctacagc', 'tgtcgatctacagc'],
                'primer_rev_sequence': ['WACTAATCAATTWCCAAATCCTCC', 'WACTAATCAATTWCCAAATCCTCC'],
                'fasta_file_name': [fasta_file_name, fasta_file_name],
             }
        )

        #######################################################################
        #
        # Sort Reads
        #
        #######################################################################

        # Create SortReadsRunner
        sort_reads_tsv_path = os.path.join(this_tempdir, 'sort_reads.tsv')
        sort_reads_runner = SortReadsRunner(fasta_path=fasta_path, fasta_id=fasta_id, alignement_parameters=alignement_parameters,
                                                fasta_information_df=fasta_information_df, sort_reads_tsv=sort_reads_tsv_path)

        sort_reads_runner.run()

        sort_reads_tsv_str_bak = """M00842:118:000000000-ABGKE:1:1101:15187:2443\t1\t1\t1\t5\t2\tCCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATT
M00842:118:000000000-ABGKE:1:1101:19104:2756\t1\t1\t1\t5\t2\tCCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATT
M00842:118:000000000-ABGKE:1:1101:18229:3444\t1\t1\t1\t2\t2\tCCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATT
"""
        with open(sort_reads_runner.sort_reads_tsv) as fout:
            sort_reads_tsv_str = fout.read()
        # import pdb; pdb.set_trace()
        assert sorted(sort_reads_tsv_str.split("\n")[:3]) == sorted(sort_reads_tsv_str_bak.split("\n")[:3])
