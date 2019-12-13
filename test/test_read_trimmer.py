# -*- coding: utf-8 -*-
import os
import pathlib
import sqlite3
from unittest import TestCase

import pandas

from vtam.utils.PathManager import PathManager
from vtam.utils.ReadTrimmer import ReadTrimmer


class TestReadTrimmer(TestCase):

    def setUp(self):
        self.this_tempdir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))

    def test_01(self):

        read_fasta_str = """>M00842:118:000000000-ABGKE:1:1101:18229:3444 1:N:0:5
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

        align_parameters = {'min_id': 0.8, 'minseqlength': 32, 'overhang': 0}

        this_tempdir_fwd = os.path.join(self.this_tempdir, 'fwd')
        pathlib.Path(this_tempdir_fwd).mkdir(parents=True, exist_ok=True)
        read_fasta_path = os.path.join(this_tempdir_fwd, "read.fasta")

        fasta_information_df = pandas.DataFrame(
            {
                'tag_fwd': ['tcgatcacgatgt', 'tgatcgatgatcag'],
                'primer_fwd': ['TCCACTAATCACAARGATATTGGTAC', 'TCCACTAATCACAARGATATTGGTAC'],
                'tag_rev': ['tgtcgatctacagc', 'tgtcgatctacagc'],
                'primer_rev': ['WACTAATCAATTWCCAAATCCTCC', 'WACTAATCAATTWCCAAATCCTCC'],
                'run_id': [1, 1],
                'marker_id': [1, 1],
                'biosample_id': [2, 5],
                'replicate_id': [2, 2],
                'reads_fasta_path': [read_fasta_path, read_fasta_path],
             }
        )
        with open(read_fasta_path, 'w') as fout:
            fout.write(read_fasta_str)

        tag_fwd_sequence_list = ['tcgatcacgatgt', 'tgatcgatgatcag']
        primer_fwd_sequence_list = ['TCCACTAATCACAARGATATTGGTAC', 'TCCACTAATCACAARGATATTGGTAC']

        #
        # Create ReadTrimmer
        sort_reads_runner_fwd = ReadTrimmer(reads_fasta_path=read_fasta_path, tag_sequence_list=tag_fwd_sequence_list,
                                            primer_sequence_list=primer_fwd_sequence_list,
                                            align_parameters=align_parameters, tempdir=this_tempdir_fwd)

        #######################################################################
        #
        # Test write_tag_primer_fasta
        #
        #######################################################################

        tag_primer_fasta_str_bak = """>tcgatcacgatgtTCCACTAATCACAARGATATTGGTAC
tcgatcacgatgtTCCACTAATCACAARGATATTGGTAC
>tgatcgatgatcagTCCACTAATCACAARGATATTGGTAC
tgatcgatgatcagTCCACTAATCACAARGATATTGGTAC
"""
        sort_reads_runner_fwd.write_tag_primer_fasta()
        with open(sort_reads_runner_fwd.tag_primer_fasta_path) as fout:
            tag_primer_fasta_str = fout.read()
        assert tag_primer_fasta_str == tag_primer_fasta_str_bak
#
        #######################################################################
        #
        # align_tag_primer_to_reads
        #
        #######################################################################

        vsearch_align_str_bak = """M00842:118:000000000-ABGKE:1:1101:18229:3444	tcgatcacgatgtTCCACTAATCACAARGATATTGGTAC	39	1	39	1	39	TCGATCACGATGTTCCACTAATCACAAGGATATTGGTAC
M00842:118:000000000-ABGKE:1:1101:19104:2756	tgatcgatgatcagTCCACTAATCACAARGATATTGGTAC	40	1	40	1	40	TGATCGATGATCAGTCCACTAATCACAAGGATATTGGTAC
M00842:118:000000000-ABGKE:1:1101:15187:2443	tgatcgatgatcagTCCACTAATCACAARGATATTGGTAC	40	1	40	1	40	TGATCGATGATCAGTCCACTAATCACAAGGATATTGGTAC
"""
        sort_reads_runner_fwd.align_tag_primer_to_reads()

        with open(sort_reads_runner_fwd.alignements_tsv_path) as fout:
            vsearch_align_str = fout.read()
        assert vsearch_align_str_bak.split("\n")[:3].sort() == vsearch_align_str_bak.split("\n")[:3].sort()

        #######################################################################
        #
        # alignements_with_high_quality_tsv_path
        #
        #######################################################################

        sort_reads_runner_fwd.keep_alignements_where_tag_is_well_aligned_to_read(overhang=align_parameters['overhang'])

        with open(sort_reads_runner_fwd.alignements_with_high_quality_tsv_path) as fout:
            vsearch_align_high_quality_str = fout.read()
        assert vsearch_align_str_bak.split("\n")[:3].sort() == vsearch_align_high_quality_str.split("\n")[:3].sort()

        #######################################################################
        #
        # Trim reads
        #
        #######################################################################

        reads_trimmed_5prime_fasta_str_bak = """>M00842:118:000000000-ABGKE:1:1101:19104:2756;tag_sequence=tgatcgatgatcag;primer_sequence=tccactaatcacaargatattggtac
CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATTGGAGGATTTGGAAATTGATTAGTTGCTGTAGATCGACA
>M00842:118:000000000-ABGKE:1:1101:18229:3444;tag_sequence=tcgatcacgatgt;primer_sequence=tccactaatcacaargatattggtac
CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATTGGAGGATTTGGAAATTGATTAGTTGCTGTAGATCGACA
>M00842:118:000000000-ABGKE:1:1101:15187:2443;tag_sequence=tgatcgatgatcag;primer_sequence=tccactaatcacaargatattggtac
CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATTGGAGGATTTGGTAATTGATTAGTAGCTGTAGATCGACA
"""

        sort_reads_runner_fwd.trim_reads()

        with open(sort_reads_runner_fwd.reads_trimmed_fasta_path) as fout:
            reads_trimmed_5prime_fasta_str = fout.read()
        assert sorted(reads_trimmed_5prime_fasta_str.split("\n")[:6]) == sorted(reads_trimmed_5prime_fasta_str_bak.split("\n")[:6])

        #######################################################################
        #
        # Reverse-complement fasta_path read with trimmed 5 prime
        #
        #######################################################################

        reads_reversed_5prime_trimmed_fasta_path = os.path.join(this_tempdir_fwd, "reads_5prime_trimmed_reversed.fasta")
        ReadTrimmer.fasta_file_to_reverse_complement(sort_reads_runner_fwd.reads_trimmed_fasta_path,
                                                     reads_reversed_5prime_trimmed_fasta_path)

        #######################################################################
        #
        # Trim reverse complement in one step and return final trimmed reads
        #
        #######################################################################

        tag_rev_sequence_list = ['tgtcgatctacagc', 'tgtcgatctacagc']
        primer_rev_sequence_list = ['WACTAATCAATTWCCAAATCCTCC', 'WACTAATCAATTWCCAAATCCTCC']

        this_tempdir_rev = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'rev')

        # Create ReadTrimmer with 5 prime trimmed reversed reads
        sort_reads_runner_rev = ReadTrimmer(reads_fasta_path=reads_reversed_5prime_trimmed_fasta_path, tag_sequence_list=tag_rev_sequence_list,
                                            primer_sequence_list=primer_rev_sequence_list,
                                            align_parameters=align_parameters, tempdir=this_tempdir_rev)

        sort_reads_runner_rev.trim_reads()

        reads_trimmed_final_fasta_path = os.path.join(this_tempdir_rev, "reads_trimmed.fasta")
        ReadTrimmer.fasta_file_to_reverse_complement(sort_reads_runner_rev.reads_trimmed_fasta_path,
                                                     reads_trimmed_final_fasta_path)

        reads_trimmed_final_str_bak = """>M00842:118:000000000-ABGKE:1:1101:19104:2756;tag_sequence=tgatcgatgatcag;primer_sequence=tccactaatcacaargatattggtac;tag_sequence=tgtcgatctacagc;primer_sequence=wactaatcaattwccaaatcctcc
CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATT
>M00842:118:000000000-ABGKE:1:1101:15187:2443;tag_sequence=tgatcgatgatcag;primer_sequence=tccactaatcacaargatattggtac;tag_sequence=tgtcgatctacagc;primer_sequence=wactaatcaattwccaaatcctcc
CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATT
>M00842:118:000000000-ABGKE:1:1101:18229:3444;tag_sequence=tcgatcacgatgt;primer_sequence=tccactaatcacaargatattggtac;tag_sequence=tgtcgatctacagc;primer_sequence=wactaatcaattwccaaatcctcc
CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCACTTAGATTTATTATTCGAATAGAATTAAGAACTCCTGGTAGATTTATTGGCAACGACCAAATTTATAACGTAATTGTTACATCTCATGCATTTATTATAATTTTTTTTATAGTTATACCAATCATAATT
"""

        with open(reads_trimmed_final_fasta_path) as fout:
            reads_trimmed_final_str = fout.read()
        assert sorted(reads_trimmed_final_str.split("\n")[:6]) == sorted(reads_trimmed_final_str_bak.split("\n")[:6])

        # #######################################################################
        # #
        # # Sort reads based on final fasta_path with trimmed sequence
        # #
        # #######################################################################
        #
        # sorted_reads_tsv_path = os.path.join(this_tempdir_rev, "sorted_reads.tsv")
        # ReadTrimmer.annotate_reads(reads_trimmed_final_fasta_path=reads_trimmed_final_fasta_path,
        #                            fasta_information=fasta_information_df, sort_reads_tsv_path=sorted_reads_tsv_path)
        # import pdb; pdb.set_trace()

