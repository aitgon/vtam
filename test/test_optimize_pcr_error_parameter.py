import os
import sqlite3
from unittest import TestCase

import pandas

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.utilities import create_step_tmp_dir
from wopmetabarcoding.wrapper.FilterPCRError import f10_pcr_error_run_vsearch, f10_get_maximal_pcr_error_value


class TestOptimizePcrErrorParameter(TestCase):

    def setUp(self):
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.optimize_variant_path = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "optimize_variants.tsv")
        self.db_opt_path = os.path.join(os.environ['DIR_DATA_NON_GIT'], 'wopmetabarcodin/test/test_files/db_opt.sqlite')
        self.this_step_tmp_dir = create_step_tmp_dir(__file__)
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
        #
        # variant_vsearch_query_df = read_count_df[['variant_id', 'variant_sequence']].drop_duplicates()
        # variant_vsearch_query_df.rename(columns={'variant_id': 'id', 'variant_sequence': 'sequence'}, inplace=True)
        #
        vsearch_output_df = f10_pcr_error_run_vsearch(variant_db_df=variant_vsearch_db_df, variant_usearch_global_df=variant_df, tmp_dir=self.this_step_tmp_dir)
        #
        #
        #
        read_count_unexpected_expected_ratio_max = f10_get_maximal_pcr_error_value(variant_read_count_df, vsearch_output_df)
        #
        #
        #
        self.assertTrue(read_count_unexpected_expected_ratio_max == 0.03182775567516967)
