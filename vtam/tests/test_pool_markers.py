import io
import os
import pathlib
from unittest import TestCase

import pandas

from vtam.utils.PathManager import PathManager
from vtam.CommandPoolRunMarkers import CommandPoolRunMarkers
from vtam.utils.VSearch import VSearch
from vtam.utils.VariantDFutils import VariantDFutils


class TestPoolMarkers(TestCase):

    def setUp(self):
        asv_table_str = """variant_id	marker	run	sequence_length	read_count	sample1	sample2	sample3	chimera_borderline	sequence
3	MFZR	prerun	176	9713	9712	1	0	FALSE	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
33	MFZR	prerun	174	9713	9703	10	0	FALSE	CTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
333	ZFZR	prerun	157	10000	9900	10	0	FALSE	TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
836	MFZR	prerun	176	11588	123	56	0	FALSE	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
8368	ZFZR	prerun	157	545	500	0	45	FALSE	TGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
83683	MFZR	prerun	175	484	0	28	456	FALSE	TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
"""
        asv_table_df = pandas.read_csv(io.StringIO(asv_table_str), sep="\t", header=0)
        self.asv_table_df = asv_table_df
        # Create this_tempdir
        this_tempdir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(this_tempdir).mkdir(exist_ok=True)
        # Define fasta_path fastqinfo_tsv_path
        fasta_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'variants.fa')
        # Create variant variant_read_count_input_df
        variant_df = asv_table_df[['variant_id', 'sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        # Create fasta_path file from asv_table_df
        variant_df_utils = VariantDFutils(variant_df)
        variant_df_utils.to_fasta(fasta_path, add_column='size')
        # Define vsearch output fastqinfo_tsv_path
        vsearch_output_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'centroid_out.fa')
        # Define cluster output fastqinfo_tsv_path
        vsearch_cluster_output_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'cluster.fa')
        #
        # Create object and run vsearch
        os.environ["VTAM_THREADS"] = "1"
        vsearch_parameters = {'--cluster_size': fasta_path,
                              '--clusters':  vsearch_cluster_output_path,
                              '--id': 1, '--sizein': None,
                              '--centroids': vsearch_output_path,
                              "--threads": int(os.getenv('VTAM_THREADS')),
                              }
        vsearch_cluster = VSearch(parameters = vsearch_parameters)
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
        pool_marker_runner = CommandPoolRunMarkers(self.asv_table_df)
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
        cluster_df_bak = pandas.read_csv(io.StringIO(cluster_str_bak), sep="\t", header=0)
        pandas.testing.assert_frame_equal(cluster_df, cluster_df_bak)

        ####################################################################
        #
        # tests get_pooled_marker_df
        #
        ####################################################################
        pooled_marker_bak_str = """centroid_variant_id	variant_id	run	marker	sample1	sample2	sample3	sequence
333	3,33,333	prerun	MFZR,ZFZR	1	1	0	TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
836	836,8368	prerun	MFZR,ZFZR	1	1	1	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
83683	83683	prerun	MFZR	0	1	1	TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
"""
        pooled_marker_bak_df = pandas.read_csv(io.StringIO(pooled_marker_bak_str), sep="\t", header=0)

        pooled_marker_df = pool_marker_runner.get_pooled_marker_df()
        pandas.testing.assert_frame_equal(pooled_marker_df, pooled_marker_bak_df)


