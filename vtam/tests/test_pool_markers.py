import io
import os
import sys

import pandas
import pathlib
import unittest

from vtam.CommandPoolRunMarkers import CommandPoolRunMarkers
from vtam.utils.PathManager import PathManager
from vtam.utils.RunnerVSearch import RunnerVSearch
from vtam.utils.DataframeVariant import DataframeVariant


class TestPoolMarkers(unittest.TestCase):

    def setUp(self):
        asv_table_str = """variant_id	marker_name	run_name	sequence_length	read_count	sample1	sample2	sample3	chimera_borderline	sequence
3	MFZR	prerun	176	9713	9712	1	0	FALSE	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
33	MFZR	prerun	174	9713	9703	10	0	FALSE	CTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
333	ZFZR	prerun	157	10000	9900	10	0	FALSE	TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
836	MFZR	prerun	176	11588	123	56	0	FALSE	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
8368	ZFZR	prerun	157	545	500	0	45	FALSE	TGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
83683	MFZR	prerun	175	484	0	28	456	FALSE	TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
"""
        asv_table_df = pandas.read_csv(
            io.StringIO(asv_table_str), sep="\t", header=0)
        self.asv_table_df = asv_table_df
        # Create this_tempdir
        this_tempdir = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(this_tempdir).mkdir(exist_ok=True)
        # Define fasta_path tsv_path
        fasta_path = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__),
            'variants.fa')
        # Create variant variant_read_count_input_df
        variant_df = asv_table_df[[
            'variant_id', 'sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        # Create fasta_path file from asv_table_df
        variant_df_utils = DataframeVariant(variant_df)
        variant_df_utils.to_fasta(fasta_path, add_column='size')
        # Define vsearch output tsv_path
        vsearch_output_path = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__),
            'centroid_out.fa')
        # Define cluster output tsv_path
        vsearch_cluster_output_path = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__),
            'cluster.fa')
        #
        # Create object and run_name vsearch
        os.environ["VTAM_THREADS"] = "1"
        vsearch_parameters = {'--cluster_size': fasta_path,
                              '--clusters': vsearch_cluster_output_path,
                              '--id': 1, '--sizein': None,
                              '--centroids': vsearch_output_path,
                              "--threads": int(os.getenv('VTAM_THREADS')),
                              }
        vsearch_cluster = RunnerVSearch(parameters=vsearch_parameters)
        vsearch_cluster.run()

    def test_cluster_sequences_with_vsearch(self):
        ####################################################################
        #
        # tests run_vsearch_to_cluster_sequences
        #
        ####################################################################
        vsearch_output_centroid_fasta_content_bak = """>836;size=11588
TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACC
CGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGAttttttttATAGTG
ATACCTATTATAATT
>333;size=10000
TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGG
ACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGAttttttttATAGTGATACCTATTATAATT
>83683;size=484
TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACC
CGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGAttttttttATAGTG
ATACCTATTATAATT
"""
        pool_marker_runner = CommandPoolRunMarkers(self.asv_table_df, readcounts=False)
        vsearch_output_centroid_fasta, vsearch_output_cluster_path = pool_marker_runner.run_vsearch_to_cluster_sequences()
        with open(vsearch_output_centroid_fasta) as fin:
            vsearch_output_centroid_fasta_content = fin.read()
        assert vsearch_output_centroid_fasta_content == vsearch_output_centroid_fasta_content_bak

        ####################################################################
        #
        # tests get_vsearch_clusters_to_df
        #
        ####################################################################

        cluster_df = pool_marker_runner.get_vsearch_clusters_to_df()
        cluster_str_bak = """centroid_variant_id	variant_id
836	836
836	8368
333	333
333	33
333	3
83683	83683"""
        if sys.platform.startswith('win'):
            cluster_df_bak = pandas.read_csv(
                io.StringIO(cluster_str_bak), sep="\t", header=0, dtype={'centroid_variant_id': 'int32', 'variant_id': 'int32'})
        else:
            cluster_df_bak = pandas.read_csv(
                io.StringIO(cluster_str_bak), sep="\t", header=0, dtype={'centroid_variant_id': 'int64', 'variant_id': 'int64'})
        pandas.testing.assert_frame_equal(cluster_df, cluster_df_bak)

        ####################################################################
        #
        # tests get_pooled_marker_df
        #
        ####################################################################
        pooled_marker_bak_str = """variant	pooled_variants	run	marker	sample1	sample2	sample3	pooled_sequences	sequence
333	3,33,333	prerun	MFZR,ZFZR	1	1	0	CTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT,TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT,TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT	TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
836	836,8368	prerun	MFZR,ZFZR	1	1	1	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT,TGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
83683	83683	prerun	MFZR	0	1	1	TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT	TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT"""
        pooled_marker_bak_df = pandas.read_csv(
            io.StringIO(pooled_marker_bak_str), sep="\t", header=0)

        pooled_marker_df = pool_marker_runner.get_pooled_marker_df()

        pandas.testing.assert_frame_equal(
            pooled_marker_df, pooled_marker_bak_df)
