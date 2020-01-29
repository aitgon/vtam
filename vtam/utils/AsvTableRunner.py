import math

import pandas
import sqlalchemy

from vtam.utils.TaxAssignRunner import f04_1_tax_id_to_taxonomy_lineage
from vtam.utils.VariantReadCountDF import VariantReadCountDF
from vtam.utils.constants import rank_hierarchy_asv_table


class AsvTableRunner(object):


    def __init__(self, engine, variant_read_count_df, variant_df, run_df, marker_df, biosample_df, variant_to_chimera_borderline_df):

        self.engine = engine
        self.variant_read_count_df = variant_read_count_df
        self.variant_df = variant_df
        self.run_df = run_df
        self.marker_df = marker_df
        self.biosample_df = biosample_df
        self.variant_to_chimera_borderline_df = variant_to_chimera_borderline_df

    def run(self):

        # Aggregate replicates
        variant_read_count_obj = VariantReadCountDF(self.variant_read_count_df)
        N_ij_df = variant_read_count_obj.get_N_ij_df()

        asv_df = N_ij_df.merge(self.biosample_df, left_on='biosample_id', right_index=True)
        asv_df.rename({'name': 'biosample_name'}, axis=1, inplace=True)
        asv_df.drop('biosample_id', axis=1, inplace=True)
        asv_df = asv_df.pivot_table(index=['run_id', 'marker_id', 'variant_id'], columns='biosample_name', values='N_ij',
                                     fill_value=0).reset_index()

        ################################################################################################################
        #
        # Asv_df2: biosamples
        #
        ################################################################################################################

        asv_df2 = asv_df
        biosample_name_list = self.biosample_df.name.tolist()
        asv_df2_columns = ['variant_id', 'marker_id', 'run_id'] + [col for col in biosample_name_list if col in asv_df2.iloc[:, 3:].columns.tolist()]
        asv_df2 = asv_df2[asv_df2_columns]

        ################################################################################################################
        #
        # Asv_df1: First part
        #
        ################################################################################################################

        asv_df1 = asv_df
        asv_df1['read_count'] = asv_df1.iloc[:, 3:].apply(sum, axis=1)

        # Add marker_name
        asv_df1 = asv_df1.merge(self.marker_df, left_on = 'marker_id', right_index=True)
        asv_df1.rename({'name': 'marker_name'}, axis=1, inplace=True)
        # asv_df1.drop('marker_id', axis=1, inplace=True)

        # Add run_name
        asv_df1 = asv_df1.merge(self.run_df, left_on='run_id', right_index=True)
        asv_df1.rename({'name': 'run_name'}, axis=1, inplace=True)
        # asv_df1.drop('run_id', axis=1, inplace=True)

        # Add sequence
        asv_df1 = asv_df1.merge(self.variant_df, left_on='variant_id', right_index=True, validate='one_to_one')
        asv_df1['sequence_length'] = asv_df1.sequence.apply(lambda x: len(x))

        asv_df1 = asv_df1[['variant_id', 'marker_id', 'run_id', 'marker_name', 'run_name', 'sequence_length', 'read_count']]

        ################################################################################################################
        #
        # Asv_df3: Last part
        #
        ################################################################################################################

        asv_df3 = asv_df[['variant_id', 'marker_id', 'run_id']]

        # variant_to_chimera_borderline_df = self.get_chimera_borderline_df()
        asv_df3 = asv_df3.merge(self.variant_to_chimera_borderline_df, on=['run_id', 'marker_id', 'variant_id'])

        # Add sequence
        asv_df3 = asv_df3.merge(self.variant_df, left_on='variant_id', right_index=True)

        # Column order
        asv_df3 = asv_df3[['variant_id', 'marker_id', 'run_id', 'chimera_borderline', 'sequence']]

        ################################################################################################################
        #
        # Asv_df: Final merge
        #
        ################################################################################################################

        asv_df_final = asv_df1.merge(asv_df2, on=['variant_id', 'marker_id', 'run_id'])
        asv_df_final = asv_df_final.merge(asv_df3, on=['variant_id', 'marker_id', 'run_id'])
        asv_df_final.drop('marker_id', axis=1, inplace=True)
        asv_df_final.drop('run_id', axis=1, inplace=True)

        return asv_df_final
