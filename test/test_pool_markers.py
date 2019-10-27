# -*- coding: utf-8 -*-
import io
import os
import pandas

from vtam.utils.VSearch import VSearch2
from vtam.utils.VSearchClusterOutput import VSearchClusterOutput
from vtam.utils.FastaIO import FastaIO
from vtam.utils.PathManager import PathManager
from unittest import TestCase




class TestPoolMarkers(TestCase):

    def setUp(self):
        otu_table = """variant_id	marker_name	run_name	sequence_length	read_count	Tpos1_prerun	phylum	class	order	family	genus	species	ltg_tax_id	ltg_tax_name	identity	ltg_rank	chimera_borderline	variant_sequence
3	MFZR	prerun	176	9713	9713	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
33	MFZR	prerun	174	9713	9713	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	CTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
333	ZFZR	prerun	157	10000	10000	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
836	MFZR	prerun	176	11588	11588	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
8368	ZFZR	prerun	157	1000	1000	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
83683	MFZR	prerun	175	100000	100000	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT"""
        otu_table_df = pandas.read_csv(io.StringIO(otu_table), sep="\t", header=0)
        # Create tempdir
        PathManager.mkdir_p(os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__)))
        # Define fasta path
        fasta_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'variants.fa')
        # Create variant df
        variant_df = otu_table_df[['variant_id', 'variant_sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        # Create fasta file from otu_table_df
        FastaIO.variant_df_to_fasta_file(variant_df, fasta_path, add_column='size')
        # Define vsearch output path
        vsearch_output_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'centroid_out.fa')
        # Define cluster output path
        vsearch_cluster_output_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'cluster.fa')
        #
        # Create object and run vsearch
        vsearch_parameters = {'--cluster_size': fasta_path, '--clusters':  vsearch_cluster_output_path,
                              '--id': 1, '--sizein': None,
                              '--centroids': vsearch_output_path}
        vsearch_cluster = VSearch2(parameters = vsearch_parameters)
        vsearch_cluster.run()
        #
        # Analyze vsearch cluster output
        vsearch_cluster = VSearchClusterOutput(clusters_path = vsearch_cluster_output_path)
        cluster_df = vsearch_cluster.clusters_to_df()
        #
        # Merge cluster_df with otu_table
        cluster_df = cluster_df.merge(otu_table_df, on='variant_id')
        cluster_df.sort_values(by=cluster_df.columns.tolist(), inplace=True)
        import pdb; pdb.set_trace()


    def test_01(self):
        pass

