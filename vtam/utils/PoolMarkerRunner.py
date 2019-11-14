import os
import pathlib
import sys

import pandas
from Bio import SeqIO

from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.VSearch import VSearch
from vtam.utils.VariantDFutils import VariantDFutils
from vtam.utils.VTAMexception import VTAMexception


class PoolMarkerRunner(object):
    """Class for the Pool Marker wrapper"""

    def __init__(self, asv_table_df):

        try:
            assert asv_table_df.columns.tolist()[:5] == ['variant_id', 'marker_name', 'run_name', 'sequence_length',
                                                         'read_count']
            assert asv_table_df.columns.tolist()[-12:] == ['phylum', 'class', 'order', 'family', 'genus', 'species', 'ltg_tax_id', 'ltg_tax_name', 'identity',
             'ltg_rank', 'chimera_borderline', 'sequence']
        except:
            Logger.instance().error(VTAMexception("The ASV table structure is wrong. It is expected to start with these columns:"
                                                  "'variant_id', 'marker_name', 'run_name', 'sequence_length', 'read_count'"
                                                  " followed by biosample names and ending with"
                                                  "'phylum', 'class', 'order', 'family', 'genus', 'species', 'ltg_tax_id', "
                                                  "'ltg_tax_name', 'identity', 'ltg_rank', 'chimera_borderline', 'sequence'."))
            sys.exit(1)

        self.biosample_names = asv_table_df.columns.tolist()[5:-12]

        # self.asv_table_df = pandas.read_csv(asv_table_tsv, sep="\t", header=0)
        self.asv_table_df = asv_table_df
        self.tmp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        PathManager.mkdir_p(self.tmp_dir)

        self.cluster_path = None # returned by cluster_sequences_with_vsearch

        self.cluster_df = None # returned by get_vsearch_clusters_to_df

    def cluster_sequences_with_vsearch(self):
        # Define fasta_path path
        fasta_path = os.path.join(self.tmp_dir, 'variants.fa')
        # Create variant df
        variant_df = self.asv_table_df[['variant_id', 'sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        variant_df.set_index('id', inplace=True)
        # Create fasta_path file from asv_table_df
        variant_df_utils = VariantDFutils(variant_df)
        variant_df_utils.to_fasta(fasta_path, add_column='size')
        # Define vsearch output path
        vsearch_output_centroid_fasta = os.path.join(self.tmp_dir, 'centroid.fa')
        # Define cluster output path
        vsearch_output_cluster_path = os.path.join(self.tmp_dir, 'cluster.fa')

        #
        # Create object and run vsearch
        vsearch_parameters = {'--cluster_size': fasta_path, '--clusters':  vsearch_output_cluster_path,
                              '--id': 1, '--sizein': None,
                              '--centroids': vsearch_output_centroid_fasta}
        vsearch_cluster = VSearch(parameters = vsearch_parameters)
        vsearch_cluster.run()
        self.cluster_path = vsearch_output_cluster_path
        return vsearch_output_centroid_fasta, vsearch_output_cluster_path

    def get_vsearch_clusters_to_df(self):

        """
        Analysis vsearch cluster output, which a path that corresponds to the same path with ticker 0, 1, 2

        For instance, if self.cluster_path=/tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa,
        then there are /tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa0, ...1, ...2, etc with the different clusters

        :return: pandas.DataFrame with columns: variant_id_centroid and variant_id
        """
        if self.cluster_path is None:
            self.cluster_sequences_with_vsearch()

        clusters_dir = os.path.dirname(self.cluster_path)
        clusters_basename = os.path.basename(self.cluster_path)
        cluster_i = 0
        cluster_i_path = os.path.join(clusters_dir, clusters_basename + str(cluster_i))
        cluster_instance_list = []
        centroid = None
        while pathlib.Path(cluster_i_path).exists():
            seq_j = 0
            for seq_record in SeqIO.parse(cluster_i_path, "fasta"):
                cluster_instance = {}
                variant_id = seq_record.id.split(';')[0]
                if seq_j == 0:
                    centroid = variant_id
                cluster_instance['centroid_variant_id'] = centroid
                cluster_instance['variant_id'] = variant_id
                seq_j += 1
                cluster_instance_list.append(cluster_instance)
            cluster_i += 1
            cluster_i_path = os.path.join(clusters_dir, clusters_basename + str(cluster_i))
        cluster_df = pandas.DataFrame(data=cluster_instance_list)
        cluster_df = cluster_df.astype(int)
        self.cluster_df = cluster_df
        return cluster_df

    def get_pooled_marker_df(self):

        if self.cluster_df is None:
            self.get_vsearch_clusters_to_df()

        # Merge cluster_df with asv_table
        self.cluster_df = self.cluster_df.merge(self.asv_table_df, on='variant_id')
        self.cluster_df.sort_values(by=self.cluster_df.columns.tolist(), inplace=True)
        #
        # Information to keep about each centroid
        centroid_df = self.cluster_df.loc[self.cluster_df.centroid_variant_id == self.cluster_df.variant_id,
                                     ['centroid_variant_id', 'phylum',
                                      'class', 'order', 'family', 'genus', 'species', 'ltg_tax_id', 'ltg_tax_name',
                                      'ltg_rank']]
        centroid_df.drop_duplicates(inplace=True)
        #
        # Centroid to aggregated variants
        # centroid_to_variant_id_df = self.cluster_df[['centroid_variant_id', 'variant_id']].drop_duplicates(inplace=False)
        # centroid_to_variant_id_df = centroid_to_variant_id_df.groupby('centroid_variant_id')['variant_id'].apply(lambda x: ','.join(map(str, sorted(x))))
        #
        # Centroid to aggregated variants and biosamples
        self.cluster_df_col = self.cluster_df.columns
        # biosample_columns = list(set(self.cluster_df_col[6:]).difference(set(list(self.cluster_df_col[-12:]))))
        # Drop some columns
        are_reads = lambda x: int(sum(x)>0)
        agg_dic = {}
        for k in ['variant_id', 'run_name', 'marker_name']:
            agg_dic[k] = lambda x: ','.join(map(str, sorted(list(set(x)))))
        for k in self.biosample_names:
            agg_dic[k] = are_reads
        pooled_marker_df = self.cluster_df[['centroid_variant_id', 'variant_id', 'run_name', 'marker_name'] + self.biosample_names]
        pooled_marker_df = pooled_marker_df.groupby('centroid_variant_id').agg(agg_dic).reset_index()
        pooled_marker_df = pooled_marker_df.merge(centroid_df, on='centroid_variant_id')
        return pooled_marker_df