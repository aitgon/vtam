import os
import pathlib

import pandas
from Bio import SeqIO

class VSearchClusterOutput(object):
    """Analyzes clusters parameter of vsearch"""

    def __init__(self, clusters_path):
        self.clusters_path =clusters_path




    def clusters_to_df(self):
        """
        Analysis vsearch cluster output, which a path that corresponds to the same path with ticker 0, 1, 2

        For instance, if self.clusters_path=/tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa,
        then there are /tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa0, ...1, ...2, etc with the different clusters

        :return: pandas.DataFrame with columns: variant_id_centroid and variant_id
        """
        clusters_dir = os.path.dirname(self.clusters_path)
        clusters_basename = os.path.basename(self.clusters_path)
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
        return cluster_df


