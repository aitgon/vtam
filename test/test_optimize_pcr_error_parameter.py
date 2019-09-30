import os
import sqlite3
from unittest import TestCase

import pandas

from wopmetabarcoding.utils.PathManager import PathManager


class TestOptimizePcrErrorParameter(TestCase):

    def setUp(self):
        self.__testdir_path = os.path.join(PathManager.get_module_test_path())
        self.optimize_variant_path = os.path.join(PathManager.get_module_test_path(), self.__testdir_path, "test_files", "optimize_variants.tsv")
        self.this_step_tmp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        PathManager.mkdir_p(self.this_step_tmp_dir)
        #

    def test_01(self):
        ####
        #
        #  tpos var df
        #
        ####
        variant_df_user = pandas.read_csv(self.optimize_variant_path, sep='\t', header=0)
        # import pdb; pdb.set_trace()
        variant_df_user = variant_df_user[['variant_id', 'biosample']]

        ####
        #
        # Select all variants from the positive samples
        #
        ####

        variant_df = pandas.DataFrame({'id':[], 'sequence':[]})
        for biosample_name in variant_df_user.biosample.unique().tolist():
            con = sqlite3.connect(self.db_opt_path)
            # sql = """select * from Variant as v,   """
            sql = """select distinct v.id, v.sequence from VariantReadCount as vrc, 
            BioSample as b, Variant as v where b.name="{}" and b.id=vrc.biosample_id and v.id=vrc.variant_id order by v.id""".format(biosample_name)
            variant_df = pandas.concat([variant_df, pandas.read_sql(sql=sql, con=con)], axis=0)
            con.close()

        variant_read_count_df = pandas.DataFrame({'marker_id':[], 'run_id':[], 'variant_id':[], 'biosample_id':[], 'replicate_id':[], 'read_count':[]})
        for biosample_name in variant_df_user.biosample.unique().tolist():
            con = sqlite3.connect(self.db_opt_path)
            # sql = """select * from Variant as v,   """
            sql = """select distinct vrc.marker_id, vrc.run_id, vrc.variant_id, vrc.biosample_id, vrc.replicate_id, vrc.read_count from VariantReadCount as vrc, 
            BioSample as b, Variant as v where b.name="{}" and b.id=vrc.biosample_id and v.id=vrc.variant_id order by v.id""".format(biosample_name)
            variant_read_count_df = pandas.concat([variant_read_count_df, pandas.read_sql(sql=sql, con=con)], axis=0)
            con.close()

        variant_vsearch_db_df = variant_df.loc[variant_df.id.isin(variant_df_user.variant_id.unique().tolist())][
            ['id', 'sequence']].drop_duplicates()
        variant_vsearch_db_df.rename(columns={'variant_id': 'id', 'variant_sequence': 'sequence'}, inplace=True)
