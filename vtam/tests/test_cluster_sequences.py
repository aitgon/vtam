import os
import shlex
import subprocess
import sys

import pandas
import unittest

from vtam.utils.PathManager import PathManager


class TestArgParser(unittest.TestCase):

    def setUp(self):

        self.package_path = os.path.join(PathManager.get_package_path())
        self.test_path = os.path.join(PathManager.get_test_path())
        self.outdir_path = os.path.join(self.test_path, 'outdir')

    def test01(self):

        asvtable_default_tsv = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default.tsv")

        # create fasta
        tempdir = PathManager.instance().get_tempdir()
        asvtable_df = pandas.read_csv(asvtable_default_tsv, sep='\t')

        i_fas = os.path.join(tempdir, 'cluster_input.fas')
        with open(i_fas, 'w') as fout:
            for idx,row in asvtable_df.iterrows():
                valdict = {}
                valdict['variant'] = row.variant
                valdict['read_count'] = row.read_count
                valdict['sequence'] = row.sequence
                fout.write(">{variant};size={read_count}\n{sequence}\n".format(**valdict))
        cmd = "vsearch --cluster_size cluster_input.fas --id 0.97 --otutabout otutabout.txt --clusters test"
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=tempdir)

        otutabout_path = os.path.join(tempdir, "otutabout.txt")
        otutabout_df = pandas.read_csv(otutabout_path, sep="\t")
        otutabout_df.rename({'#OTU ID': 'centroid'}, axis=1, inplace=True)

        otutabout_long_df = pandas.melt(otutabout_df, id_vars=['centroid'], var_name='variant', value_name='read_count')
        otutabout_long_df.rename({'centroid': 'clusterid'}, axis=1, inplace=True)
        otutabout_long_df = otutabout_long_df.loc[otutabout_long_df.read_count > 0]
        otutabout_long_df.variant = otutabout_long_df.variant.astype('int')

        cluster_count_df = otutabout_long_df[['clusterid', 'variant']].groupby('clusterid').count()
        cluster_count_df.rename({'variant': 'clustersize'}, axis=1, inplace=True)
        cluster_count_df = otutabout_long_df[['clusterid', 'variant']].merge(cluster_count_df, on='clusterid')

        asvtable_out_df = asvtable_df.merge(cluster_count_df, on='variant')

        cols = list(asvtable_out_df.columns.values)  # Make a list of all of the columns in the df
        cols.pop(cols.index('chimera_borderline'))  # Remove b from list
        cols.pop(cols.index('sequence'))  # Remove x from list
        asvtable_out_df = asvtable_out_df[cols + ['chimera_borderline', 'sequence']]  # Create new dataframe with columns in the order you want
        asvtable_out_df.to_csv('asvtable.tsv', index=False)



