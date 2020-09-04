import os
import shlex
import subprocess
import sys

import pandas

from vtam.utils.PathManager import PathManager

from vtam.models.Biosample import Biosample
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.utils.NameIdConverter import NameIdConverter
from vtam.utils.SequenceClusterer import SequenceClusterer
from vtam.utils.VariantReadCountLikeDF import VariantReadCountLikeDF


class AsvTableRunner(object):

    def __init__(self, variant_read_count_df, engine, biosample_list):

        self.variant_read_count_df = variant_read_count_df
        self.engine = engine
        self.biosample_list = biosample_list

    def to_tsv(self, asvtable_path):

        asvtable_df = self.create_asvtable_df()
        asvtable_df.to_csv(asvtable_path, sep="\t", header=True, index=False)

    def create_asvtable_df(self):

        """Merge asvtable information and reorder columns"""

        asvtable_variant_info_df = self.get_asvtable_variants()
        asvtable_biosamples_df = self.get_asvtable_biosamples()

        ############################################################################################
        #
        # Merge variant and biosample sides of asvtables
        #
        ############################################################################################

        asvtable_df = asvtable_variant_info_df.merge(asvtable_biosamples_df, on=['run_id', 'marker_id', 'variant_id'])
        asvtable_df.run_id = NameIdConverter(id_name_or_sequence_list=asvtable_df.run_id.tolist(), engine=self.engine) \
            .to_names(Run)
        asvtable_df.marker_id = NameIdConverter(id_name_or_sequence_list=asvtable_df.marker_id.tolist(), engine=self.engine) \
            .to_names(Marker)
        asvtable_df.rename({'run_id': 'run', 'marker_id': 'marker', 'variant_id': 'variant'}, axis=1, inplace=True)

        ############################################################################################
        #
        # Reorder columns
        #
        ############################################################################################

        column_list = asvtable_df.columns.tolist()
        column_list.remove("clusterid")
        column_list.remove("clustersize")
        column_list.remove("chimera_borderline")
        column_list.remove("sequence")
        column_list = column_list + ['clusterid', 'clustersize', 'chimera_borderline', 'sequence']

        column_list.remove("sequence_length")
        column_list.remove("read_count")
        column_list.insert(3, "sequence_length")
        column_list.insert(4, "read_count")

        asvtable_df = asvtable_df[column_list]

        return asvtable_df

    def get_asvtable_variants(self):

        """This function gets variant-related data such sequence, sequence length and chimera
        borderline status

        Example:
    run_id  marker_id  variant_id  read_count                                           sequence  sequence_length  chimera_borderline  clusterid  clustersize
0        1          1          95         563  ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCT...              175               False       1651            2
1        1          1         271        3471  ACTTTATTTTATTTTTGGTGCTTGATCAGGAATAGTAGGAACTTCT...              175               False       3679            2
2        1          1         603         489  CCTTTATCTTGTATTTGGTGCCTGGGCCGGAATGGTAGGGACCGCC...              175               False        603            1
3        1          1         954       34970  CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCA...              175               False        954            1
4        1          1        1309         755  CTTATATTTTATTTTTGGTGCTTGATCAGGGATAGTGGGAACTTCT...              175               False       1309            2
            """

        asvtable_variant_info_df = VariantReadCountLikeDF(self.variant_read_count_df).get_N_i_df()
        asvtable_variant_info_df.rename({'N_i': 'read_count'}, axis=1, inplace=True)

        asvtable_variant_info_df['sequence'] = NameIdConverter(
            id_name_or_sequence_list=asvtable_variant_info_df.variant_id.tolist(), engine=self.engine)\
            .variant_id_to_sequence()

        asvtable_variant_info_df['sequence_length'] = asvtable_variant_info_df.sequence.apply(len)

        asvtable_variant_info_df['chimera_borderline'] = NameIdConverter(
            id_name_or_sequence_list=asvtable_variant_info_df.variant_id.tolist(), engine=self.engine)\
            .variant_id_is_chimera_borderline()

        ############################################################################################
        #
        # Cluster sequences
        #
        ############################################################################################

        # tempcluster_dir = PathManager.instance().get_tempdir()
        #
        # i_fas = os.path.join(tempcluster_dir, 'cluster_input.fas')
        # with open(i_fas, 'w') as fout:
        #     for idx,row in asvtable_variant_info_df.iterrows():
        #         valdict = {}
        #         valdict['variant_id'] = row.variant_id
        #         valdict['read_count'] = row.read_count
        #         valdict['sequence'] = row.sequence
        #         fout.write(">{variant_id};size={read_count}\n{sequence}\n".format(**valdict))
        # cmd = "vsearch --cluster_size cluster_input.fas --id 0.97 --otutabout otutabout.txt --clusters test"
        # if sys.platform.startswith("win"):
        #     args = cmd
        # else:
        #     args = shlex.split(cmd)
        # subprocess.run(args=args, check=True, cwd=tempcluster_dir)
        #
        # otutabout_path = os.path.join(tempcluster_dir, "otutabout.txt")
        # otutabout_df = pandas.read_csv(otutabout_path, sep="\t")
        # otutabout_df.rename({'#OTU ID': 'centroid'}, axis=1, inplace=True)
        #
        # otutabout_long_df = pandas.melt(otutabout_df, id_vars=['centroid'], var_name='variant_id', value_name='read_count')
        # otutabout_long_df.rename({'centroid': 'clusterid'}, axis=1, inplace=True)
        # otutabout_long_df = otutabout_long_df.loc[otutabout_long_df.read_count > 0]
        # otutabout_long_df.variant_id = otutabout_long_df.variant_id.astype('int')
        #
        # cluster_count_df = otutabout_long_df[['clusterid', 'variant_id']].groupby('clusterid').count()
        # cluster_count_df.rename({'variant_id': 'clustersize'}, axis=1, inplace=True)
        # cluster_count_df = otutabout_long_df[['clusterid', 'variant_id']].merge(cluster_count_df, on='clusterid')

        seq_clusterer_obj = SequenceClusterer(asvtable_variant_info_df)
        cluster_count_df = seq_clusterer_obj.compute_clusters()

        asvtable_variant_info_df = asvtable_variant_info_df.merge(cluster_count_df, on='variant_id')

        ############################################################################################
        #
        # End of variant sequence clusters
        #
        ############################################################################################

        return asvtable_variant_info_df

    def get_asvtable_biosamples(self):

        """This function gets biosample-related data such as run, marker, variants and read_count per biosample

        Example:
biosample  run_id  marker_id  variant_id  tpos1_run1  tnegtag_run1  14ben01  14ben02
0               1          1          95         563             0        0        0
1               1          1         271        3471             0        0        0
2               1          1         603         489             0        0        0
3               1          1         954       34970             0        0        0
4               1          1        1309         755             0        0        0
5               1          1        2371        9124             0        0        0
6               1          1        2684           0             0        0     5343
7               1          1        4134           0             0        0      990
            """

        asvtable_2nd_df = VariantReadCountLikeDF(self.variant_read_count_df).get_N_ij_df()
        asvtable_2nd_df.biosample_id = NameIdConverter(id_name_or_sequence_list=asvtable_2nd_df.biosample_id.tolist(), engine=self.engine)\
            .to_names(Biosample)
        asvtable_2nd_df.rename({'biosample_id': 'biosample'}, axis=1, inplace=True)
        asvtable_2nd_df = asvtable_2nd_df.pivot_table(index=['run_id', 'marker_id', 'variant_id'], columns='biosample',
                            values='N_ij', fill_value=0).reset_index()

        ############################################################################################
        #
        # Set order and fill 0-read count biosamples with Zeros
        #
        ############################################################################################

        asvtable_2nd_2_df = (asvtable_2nd_df.copy())[['run_id', 'marker_id', 'variant_id']]
        for biosample_name in self.biosample_list:
            if biosample_name in asvtable_2nd_df.columns:  # biosample with read counts
                asvtable_2nd_2_df[biosample_name] = asvtable_2nd_df[biosample_name].tolist()
            else:  # biosample without read counts
                asvtable_2nd_2_df[biosample_name] = [0] * asvtable_2nd_df.shape[0]

        return asvtable_2nd_2_df
