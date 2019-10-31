# -*- coding: utf-8 -*-
import io
import os
import pandas

from vtam.utils.VSearch import VSearch
from vtam.utils.VSearchClusterOutput import VSearchClusterOutput
from vtam.utils.VariantDFutils import VariantDFutils
from vtam.utils.PathManager import PathManager
from unittest import TestCase




class TestPoolMarkers(TestCase):

    def setUp(self):
        otu_table = """variant_id	marker_name	run_name	sequence_length	read_count	sample1	sample2	sample3	phylum	class	order	family	genus	species	ltg_tax_id	ltg_tax_name	identity	ltg_rank	chimera_borderline	variant_sequence
3	MFZR	prerun	176	9713	9712	1	0	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
33	MFZR	prerun	174	9713	9703	10	0	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	CTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
333	ZFZR	prerun	157	10000	9900	10	0	Arthropoda	Insecta	Ephemeroptera	Baetidae	Baetis	Baetis rhodani	189839	Baetis rhodani	100	species	False	TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
836	MFZR	prerun	176	11588	123	56	0	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
8368	ZFZR	prerun	157	545	500	0	45	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
83683	MFZR	prerun	175	484	0	28	456	Chordata	Actinopteri	Cypriniformes	Cyprinidae	Phoxinus	Phoxinus phoxinus	58324	Phoxinus phoxinus	100	species	False	TCTAAATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCCCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT"""
        otu_table_df = pandas.read_csv(io.StringIO(otu_table), sep="\t", header=0)
        # Create tempdir
        PathManager.mkdir_p(os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__)))
        # Define fasta path
        fasta_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'variants.fa')
        # Create variant df
        variant_df = otu_table_df[['variant_id', 'variant_sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        # Create fasta file from otu_table_df
        variant_df_utils = VariantDFutils(variant_df)
        variant_df_utils.to_fasta(fasta_path, add_column='size')
        # Define vsearch output path
        vsearch_output_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'centroid_out.fa')
        # Define cluster output path
        vsearch_cluster_output_path = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__), 'cluster.fa')
        #
        # Create object and run vsearch
        vsearch_parameters = {'--cluster_size': fasta_path, '--clusters':  vsearch_cluster_output_path,
                              '--id': 1, '--sizein': None,
                              '--centroids': vsearch_output_path}
        vsearch_cluster = VSearch(parameters = vsearch_parameters)
        vsearch_cluster.run()
        #
        # Analyze vsearch cluster output
        vsearch_cluster = VSearchClusterOutput(clusters_path = vsearch_cluster_output_path)
        cluster_df = vsearch_cluster.clusters_to_df()
        #
        # Merge cluster_df with otu_table
        cluster_df = cluster_df.merge(otu_table_df, on='variant_id')
        cluster_df.sort_values(by=cluster_df.columns.tolist(), inplace=True)

        # ONGOING
        # Information to keep about each centroid
        centroid_df = cluster_df.loc[cluster_df.centroid_variant_id == cluster_df.variant_id,
                                     ['centroid_variant_id', 'phylum',
                                      'class', 'order', 'family', 'genus', 'species', 'ltg_tax_id', 'ltg_tax_name',
                                      'ltg_rank']]
        centroid_df.drop_duplicates(inplace=True)
        #
        # Centroid to aggregated variants
        centroid_to_variant_id_df = cluster_df[['centroid_variant_id', 'variant_id']].drop_duplicates(inplace=False)
        centroid_to_variant_id_df = centroid_to_variant_id_df.groupby('centroid_variant_id')['variant_id'].apply(lambda x: ','.join(map(str, sorted(x))))
        #
        # Centroid to aggregated variants and biosamples
        cluster_df_col = cluster_df.columns
        biosample_columns = list(set(cluster_df_col[6:]).difference(set(list(cluster_df_col[-12:]))))
        # Drop some columns
        are_reads = lambda x: int(sum(x)>0)
        # agg_dic = {'variant_id': lambda x: ','.join(map(str, x)), 'sample1': are_reads, 'sample2': are_reads, 'sample3': are_reads}
        agg_dic = {}
        for k in ['variant_id', 'run_name', 'marker_name']:
            agg_dic[k] = lambda x: ','.join(map(str, list(set(x))))
        for k in biosample_columns:
            agg_dic[k] = are_reads
        output_df = cluster_df[['centroid_variant_id', 'variant_id', 'run_name', 'marker_name'] + biosample_columns]
        output_df = output_df.groupby('centroid_variant_id').agg(agg_dic).reset_index()
        output_df = output_df.merge(centroid_df, on='centroid_variant_id')
        output_df_cols = list(output_df.columns)
        # output_df_cols.insert(0, output_df_cols.pop(output_df_cols.index('marker'))) #
        # output_df = output_df[]
        output_df.to_csv("otu_clusters.tsv", sep="\t", index=False)
        # import pdb; pdb.set_trace()


    def test_01(self):
        pass

