import math

import pandas
import sqlalchemy

from vtam.utils.TaxAssignRunner import f04_1_tax_id_to_taxonomy_lineage
from vtam.utils.VariantReadCountDF import VariantReadCountDF
from vtam.utils.constants import rank_hierarchy_asv_table
from vtam.models.TaxAssign import TaxAssign as tax_assign_declarative


class AsvTableRunner(object):


    def __init__(self, engine, variant_read_count_df, variant_df, run_df, marker_df, biosample_df, variant_to_chimera_borderline_df,
                 taxonomy_tsv=None):

        self.engine = engine
        self.variant_read_count_df = variant_read_count_df
        self.variant_df = variant_df
        self.run_df = run_df
        self.marker_df = marker_df
        self.biosample_df = biosample_df
        self.variant_to_chimera_borderline_df = variant_to_chimera_borderline_df
        # self.tax_assign_model = tax_assign_model
        self.taxonomy_tsv = taxonomy_tsv

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

        ################################################################################################################
        #
        # Asv_df3: If self.taxonomy_tsv, then there is taxonomic assignation
        #
        ################################################################################################################

        if not (self.taxonomy_tsv is None):

            #####
            #
            # taxonomy_db to df
            #
            #####
            taxonomy_tsv_path = self.taxonomy_tsv
            taxonomy_df = pandas.read_csv(taxonomy_tsv_path, sep="\t", header=0,
                                             dtype={'tax_id': 'int', 'parent_tax_id': 'int', 'old_tax_id': 'float'})
            # taxonomy_df = f01_taxonomy_tsv_to_df(taxonomy_tsv_path)

            #####
            #
            # ltg_tax_assign
            #
            #####

            # Select to DataFrame
            tax_assign_list = []
            # tax_assign_model_table = self.tax_assign_model.__table__
            tax_assign_model_table = tax_assign_declarative.__table__
            with self.engine.connect() as conn:
                for df_row in self.variant_df.itertuples():
                    variant_id = df_row.Index
                    stmt_ltg_tax_assign = sqlalchemy.select([tax_assign_model_table.c.variant_id,
                                                             tax_assign_model_table.c.identity,
                                                             tax_assign_model_table.c.ltg_rank,
                                                             tax_assign_model_table.c.ltg_tax_id,
                                                             tax_assign_model_table.c.blast_db])\
                        .where(tax_assign_model_table.c.variant_id == variant_id)

                    try:
                        variant_id, identity, ltg_rank, ltg_tax_id, blast_db = conn.execute(stmt_ltg_tax_assign).first()
                        tax_assign_list.append({'variant_id': variant_id, 'identity': identity, 'ltg_rank': ltg_rank, 'ltg_tax_id': ltg_tax_id, 'blast_db': blast_db})
                    except TypeError: # no result
                        pass
            ltg_tax_assign_df = pandas.DataFrame.from_records(tax_assign_list, index='variant_id')
            #
            ltg_tax_assign_df = ltg_tax_assign_df.reset_index().merge(taxonomy_df,
                                                                      left_on='ltg_tax_id', right_on='tax_id', how="left").set_index('variant_id')
            ltg_tax_assign_df.drop(['tax_id', 'parent_tax_id', 'rank', 'old_tax_id'], axis=1, inplace=True)
            ltg_tax_assign_df = ltg_tax_assign_df.rename(columns={'name_txt': 'ltg_tax_name'})

            #
            # Merge ltg tax assign results
            asv_df = asv_df.merge(ltg_tax_assign_df, left_on='variant_id', right_index=True, how='left').drop_duplicates(inplace=False)
            list_lineage = []
            for tax_id in asv_df['ltg_tax_id'].unique().tolist():
                if not math.isnan(tax_id):
                    dic_lineage = f04_1_tax_id_to_taxonomy_lineage(tax_id, taxonomy_df, give_tax_name=True)
                    list_lineage.append(dic_lineage)
            lineage_df = pandas.DataFrame(data=list_lineage)
            lineage_list_df_columns_sorted = list(
                filter(lambda x: x in lineage_df.columns.tolist(), rank_hierarchy_asv_table))
            lineage_list_df_columns_sorted = lineage_list_df_columns_sorted + ['tax_id']
            lineage_df = lineage_df[lineage_list_df_columns_sorted]
            asv_df3 = asv_df3.merge(ltg_tax_assign_df, left_on='variant_id', right_index=True, how='left').drop_duplicates(inplace=False)
            asv_df3 = asv_df3.merge(lineage_df, left_on='ltg_tax_id', right_on='tax_id', how='left').drop_duplicates(inplace=False)
            asv_df3.drop('tax_id', axis=1, inplace=True)

        # Add sequence
        asv_df3 = asv_df3.merge(self.variant_df, left_on='variant_id', right_index=True)

        # Column order
        if not (self.taxonomy_tsv is None):
            asv_df3_columns = ['variant_id', 'marker_id', 'run_id', 'phylum', 'class', 'order', 'family', 'genus',
                               'species', 'ltg_tax_id', 'ltg_tax_name', 'identity', 'blast_db', 'ltg_rank',
                               'chimera_borderline', 'sequence']
        else:
            asv_df3_columns = ['variant_id', 'marker_id', 'run_id', 'chimera_borderline', 'sequence']

        asv_df3 = asv_df3[asv_df3_columns]

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
