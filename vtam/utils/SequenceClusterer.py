import os
import shlex
import subprocess
import sys

import pandas
from vtam.utils.PathManager import PathManager


class SequenceClusterer(object):

    def __init__(self, variant_info_df, cluster_identity):

        """Takes as input df with at least these columns: variant_id, read_cout, sequence"""

        self.variant_info_df = variant_info_df
        self.cluster_identity = cluster_identity

    def compute_clusters(self):

        tempcluster_dir = PathManager.instance().get_tempdir()

        i_fas = os.path.join(tempcluster_dir, 'cluster_input.fas')
        with open(i_fas, 'w') as fout:
            for idx, row in self.variant_info_df.iterrows():
                valdict = {}
                valdict['variant_id'] = row.variant_id
                valdict['read_count'] = row.read_count
                valdict['sequence'] = row.sequence
                fout.write(
                    ">{variant_id};size={read_count}\n{sequence}\n".format(
                        **valdict))
        cmd = "vsearch --cluster_size cluster_input.fas --id {} --otutabout otutabout.txt --clusters test".format(self.cluster_identity)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=tempcluster_dir)

        otutabout_path = os.path.join(tempcluster_dir, "otutabout.txt")
        otutabout_df = pandas.read_csv(otutabout_path, sep="\t")
        otutabout_df.rename({'#OTU ID': 'centroid'}, axis=1, inplace=True)

        otutabout_long_df = pandas.melt(otutabout_df, id_vars=['centroid'],
                                        var_name='variant_id',
                                        value_name='read_count')
        otutabout_long_df.rename({'centroid': 'clusterid'}, axis=1,
                                 inplace=True)
        otutabout_long_df = otutabout_long_df.loc[
            otutabout_long_df.read_count > 0]
        otutabout_long_df.variant_id = otutabout_long_df.variant_id.astype(
            'int')

        cluster_count_df = otutabout_long_df[
            ['clusterid', 'variant_id']].groupby('clusterid').count()
        cluster_count_df.rename({'variant_id': 'clustersize'}, axis=1,
                                inplace=True)
        cluster_count_df = otutabout_long_df[
            ['clusterid', 'variant_id']].merge(cluster_count_df,
                                               on='clusterid')

        return cluster_count_df
