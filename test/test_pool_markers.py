import io
import os
from unittest import TestCase

import pandas

from vtam.utils.PathManager import PathManager
from vtam.utils.PoolMarkerRunner import PoolMarkerRunner
from vtam.utils.VSearch import VSearch
from vtam.utils.VariantDFutils import VariantDFutils


class TestPoolMarkers(TestCase):

    def setUp(self):
        asv_table_str = """variant_id	marker_name	run_name	sequence_length	read_count	sample1	sample2	sample3	phylum	class	order	family	genus	species	ltg_tax_id	ltg_tax_name	identity	ltg_rank	chimera_borderline	sequence
3	MFZR	prerun	176	9713	9712	1	0	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
33	MFZR	prerun	174	9713	9703	10	0	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	CTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
333	ZFZR	prerun	157	10000	9900	10	0	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
836	MFZR	prerun	176	11588	123	56	0	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
8368	ZFZR	prerun	157	545	500	0	45	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
83683	MFZR	prerun	175	484	0	28	456	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT"""
        asv_table_df = pandas.read_csv(io.StringIO(asv_table_str), sep="\t", header=0)
        self.asv_table_df = asv_table_df
        # Create tempdir
        PathManager.mkdir_p(os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__)))
        # Define fasta_path path
        fasta_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'variants.fa')
        # Create variant df
        variant_df = asv_table_df[['variant_id', 'sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        # Create fasta_path file from asv_table_df
        variant_df_utils = VariantDFutils(variant_df)
        variant_df_utils.to_fasta(fasta_path, add_column='size')
        # Define vsearch output path
        vsearch_output_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'centroid_out.fa')
        # Define cluster output path
        vsearch_cluster_output_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'cluster.fa')
        #
        # Create object and run vsearch
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
        # test run_vsearch_to_cluster_sequences
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
        pool_marker_runner = PoolMarkerRunner(self.asv_table_df)
        vsearch_output_centroid_fasta, vsearch_output_cluster_path = pool_marker_runner.run_vsearch_to_cluster_sequences()
        with open(vsearch_output_centroid_fasta) as fin:
            vsearch_output_centroid_fasta_content = fin.read()
        assert vsearch_output_centroid_fasta_content == vsearch_output_centroid_fasta_content_bak

        ####################################################################
        #
        # test get_vsearch_clusters_to_df
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
        # test get_pooled_marker_df
        #
        ####################################################################
        pooled_marker_bak_str = """centroid_variant_id	variant_id	run_name	marker_name	sample1	sample2	sample3	phylum	class	order	family	genus	species	ltg_tax_id	ltg_tax_name	ltg_rank
333	3,33,333	prerun	MFZR,ZFZR	1	1	0	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	species
836	836,8368	prerun	MFZR,ZFZR	1	1	1	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	species
83683	83683	prerun	MFZR	0	1	1	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	species
"""
        pooled_marker_bak_df = pandas.read_csv(io.StringIO(pooled_marker_bak_str), sep="\t", header=0)

        pooled_marker_df = pool_marker_runner.get_pooled_marker_df()
        pandas.testing.assert_frame_equal(pooled_marker_df, pooled_marker_bak_df)


