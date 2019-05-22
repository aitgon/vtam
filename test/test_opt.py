import os
import sqlite3
from math import floor
from unittest import TestCase

import pandas

from wopmetabarcoding.utils.VSearch import VSearch1

from wopmetabarcoding.utils.utilities import create_step_tmp_dir


class TestOpt(TestCase):

    def setUp(self):
        self.db_opt_path = os.path.join(os.environ['DIR_DATA_NON_GIT'], 'wopmetabarcodin/test/test_files/db_opt.sqlite')
        self.tpos_path = os.path.join(os.environ['DIR_DATA_NON_GIT'], 'wopmetabarcodin/test/test_files/tpos.csv')
        self.this_step_tmp_dir = create_step_tmp_dir(__file__)
        #

    def test_01_opt(self):

        ####
        #
        #  Variant df
        #
        ####

        con = sqlite3.connect(self.db_opt_path)
        sql = """select * from Variant  """
        variant_df = pandas.read_sql(sql=sql, con=con)
        con.close()


        ####
        #
        #  tpos var df
        #
        ####
        tpos_variant_df = pandas.read_csv(self.tpos_path, sep=';', header=0)
        tpos_var_df = tpos_variant_df[['id', 'Sequence']]

        # Vsearch output

        Vsearch_output_df = run_vsearch(variant_df, tpos_var_df, self.this_step_tmp_dir)
        Vsearch_output_df['sum_mism_gaps'] = Vsearch_output_df.mism + Vsearch_output_df.gaps

        import pdb;
        pdb.set_trace()

def run_vsearch(variant_df, tpos_var_df,this_step_tmp_dir):
    # length of smallest sequence
    length_min = min(variant_df.sequence.apply(len).tolist())
    # calcul identity
    identity = floor((length_min - 1) / length_min * 100) / 100
    # identity = 0.99

    #
    ###################################################################
    # 1. Make a fasta file with all variants and other one with tpos variant
    ###################################################################

    variant_all_fasta = os.path.join(this_step_tmp_dir, '{}.fasta'.format("pcr_error"))
    with open(variant_all_fasta, 'w') as fout:
        for row in variant_df.itertuples():
            id = row.id
            sequence = row.sequence
            fout.write(">{}\n{}\n".format(id, sequence))

    variant_tpos_fasta = os.path.join(this_step_tmp_dir, '{}.fasta'.format("pcr_error"))
    with open(variant_tpos_fasta, 'w') as fout:
        for row in tpos_var_df.itertuples():
            id = row.id
            sequence = row.Sequence
            fout.write(">{}\n{}\n".format(id, sequence))

        ###################################################################
    # import pdb; pdb.set_trace()
    sample_tsv = os.path.join(this_step_tmp_dir, '{}.tsv'.format("pcr_error"))
    vsearch_usearch_global_args = {'db': variant_tpos_fasta,
                                   'usearch_global': variant_tpos_fasta,
                                   'id': str(identity),
                                   'maxrejects': 0,
                                   'maxaccepts': 0,
                                   'userout': sample_tsv,
                                   'userfields': "query+target+alnlen+ids+mism+gaps",
                                   }

    vsearch_usearch_global = VSearch1(**vsearch_usearch_global_args)
    vsearch_usearch_global.run()
    column_names = ['query', 'target', 'alnlen', 'ids', 'mism', 'gaps']
    vsearch_output_df = pandas.read_csv(sample_tsv, sep='\t', names=column_names)
    return vsearch_output_df


