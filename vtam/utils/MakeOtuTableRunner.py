import inspect

import pandas
import sqlalchemy

from vtam.utils.Logger import Logger
from vtam.model.Biosample import Biosample
from vtam.utils.FastaInformation import FastaInformation
from vtam.utils.TaxAssignUtilities import f01_taxonomy_sqlite_to_df, f04_1_tax_id_to_taxonomy_lineage
from vtam.utils.VariantReadCountDF import VariantReadCountDF
from vtam.utils.constants import rank_hierarchy_otu_table


class MakeOtuTableRunner(object):

    def __init__(self, engine, input_file_fastainfo, run_model, marker_model, biosample_model, replicate_model,
                 filter_chimera_borderline_model, filter_codon_stop_model, variant_model, tax_assign_model, input_file_taxonomy):

        self.engine = engine
        self.input_file_fastainfo = input_file_fastainfo
        self.run_model = run_model
        self.marker_model = marker_model
        self.biosample_model = biosample_model
        self.replicate_model = replicate_model
        self.filter_chimera_borderline_model = filter_chimera_borderline_model
        self.filter_codon_stop_model = filter_codon_stop_model
        self.variant_model = variant_model
        self.tax_assign_model = tax_assign_model
        self.input_file_taxonomy = input_file_taxonomy

    def run(self):

        fasta_info_obj = FastaInformation(self.input_file_fastainfo, self.engine, self.run_model, self.marker_model, self.biosample_model, self.replicate_model)
        variant_read_count_df = fasta_info_obj.get_variant_read_count_df(self.filter_codon_stop_model)
        variant_df = fasta_info_obj.get_variant_df(variant_read_count_like_model=self.filter_codon_stop_model,
                                               variant_model=self.variant_model)

        biosample_df = fasta_info_obj.get_biosample_df(variant_read_count_like_model=self.filter_codon_stop_model)
        marker_df = fasta_info_obj.get_marker_df(variant_read_count_like_model=self.filter_codon_stop_model)
        run_df = fasta_info_obj.get_run_df(variant_read_count_like_model=self.filter_codon_stop_model)

        # Aggregate replicates
        variant_read_count_obj = VariantReadCountDF(variant_read_count_df)
        N_ij_df = variant_read_count_obj.get_N_ij_df()

        otu_df = N_ij_df.merge(biosample_df, left_on='biosample_id', right_on='id')
        otu_df.rename({'name': 'biosample_name'}, axis=1, inplace=True)
        otu_df.drop('biosample_id', axis=1, inplace=True)
        otu_df = otu_df.pivot_table(index=['run_id', 'marker_id', 'variant_id'], columns='biosample_name', values='N_ij',
                                     fill_value=0).reset_index()

        ################################################################################################################
        #
        # Otu_df2: biosamples
        #
        ################################################################################################################

        otu_df2 = otu_df
        biosample_name_list = fasta_info_obj.df.biosample_name.drop_duplicates(keep='first').tolist()
        otu_df2_columns = ['variant_id', 'marker_id', 'run_id'] + [col for col in biosample_name_list if col in otu_df2.iloc[:, 3:].columns.tolist()]
        otu_df2 = otu_df2[otu_df2_columns]

        ################################################################################################################
        #
        # Otu_df1: First part
        #
        ################################################################################################################

        otu_df1 = otu_df
        otu_df1['read_count'] = otu_df1.iloc[:, 3:].apply(sum, axis=1)

        # Add marker_name
        otu_df1 = otu_df1.merge(marker_df, left_on = 'marker_id', right_index=True)
        otu_df1.rename({'name': 'marker_name'}, axis=1, inplace=True)
        # otu_df1.drop('marker_id', axis=1, inplace=True)

        # Add run_name
        otu_df1 = otu_df1.merge(run_df, left_on='run_id', right_index=True)
        otu_df1.rename({'name': 'run_name'}, axis=1, inplace=True)
        # otu_df1.drop('run_id', axis=1, inplace=True)

        # Add sequence
        otu_df1 = otu_df1.merge(variant_df, left_on='variant_id', right_index=True, validate='one_to_one')
        otu_df1['sequence_length'] = otu_df1.sequence.apply(lambda x: len(x))

        otu_df1 = otu_df1[['variant_id', 'marker_id', 'run_id', 'marker_name', 'run_name', 'sequence_length', 'read_count']]

        ################################################################################################################
        #
        # Otu_df3: Last part
        #
        ################################################################################################################

        otu_df3 = otu_df[['variant_id', 'marker_id', 'run_id']]

        # Chimera borderline
        variant_to_chimera_borderline_df = self.get_chimera_borderline_df()
        otu_df3 = otu_df3.merge(variant_to_chimera_borderline_df, on=['variant_id', 'marker_id', 'run_id'])

        #####
        #
        # ltg_tax_assign
        #
        #####

        tax_assign_model_table = self.tax_assign_model.__table__
        stmt_ltg_tax_assign = sqlalchemy.select([tax_assign_model_table.c.variant_id,
                               tax_assign_model_table.c.identity,
                               tax_assign_model_table.c.ltg_rank,
                               tax_assign_model_table.c.ltg_tax_id])

        #####
        #
        # taxonomy_db to df
        #
        #####
        taxonomy_sqlite_path = self.input_file_taxonomy
        taxonomy_db_df = f01_taxonomy_sqlite_to_df(taxonomy_sqlite_path)


        # Select to DataFrame
        tax_assign_list = []
        with self.engine.connect() as conn:
            for row in conn.execute(stmt_ltg_tax_assign).fetchall():
                tax_assign_list.append(row)
        ltg_tax_assign_df = pandas.DataFrame.from_records(tax_assign_list,
                                                          columns=['variant_id', 'identity', 'ltg_rank', 'ltg_tax_id'])
        #
        ltg_tax_assign_df = ltg_tax_assign_df.merge(taxonomy_db_df, left_on='ltg_tax_id', right_on='tax_id')
        ltg_tax_assign_df.drop(['tax_id', 'parent_tax_id', 'rank', 'old_tax_id'], axis=1, inplace=True)
        ltg_tax_assign_df = ltg_tax_assign_df.rename(columns={'name_txt': 'ltg_tax_name'})

        #
        # Merge ltg tax assign results
        otu_df = otu_df.merge(ltg_tax_assign_df, on='variant_id')
        list_lineage = []
        for tax_id in otu_df['ltg_tax_id'].unique().tolist():
            dic_lineage = f04_1_tax_id_to_taxonomy_lineage(tax_id, taxonomy_db_df, give_tax_name=True)
            list_lineage.append(dic_lineage)
        lineage_df = pandas.DataFrame(data=list_lineage)
        lineage_list_df_columns_sorted = list(
            filter(lambda x: x in lineage_df.columns.tolist(), rank_hierarchy_otu_table))
        lineage_list_df_columns_sorted = lineage_list_df_columns_sorted + ['tax_id']
        lineage_df = lineage_df[lineage_list_df_columns_sorted]

        otu_df3 = otu_df3.merge(ltg_tax_assign_df, on='variant_id')
        otu_df3 = otu_df3.merge(lineage_df, left_on='ltg_tax_id', right_on='tax_id')
        otu_df3.drop('tax_id', axis=1, inplace=True)

        # Add sequence
        otu_df3 = otu_df3.merge(variant_df, left_on='variant_id', right_index=True, validate='one_to_one')

        # Column order
        otu_df3 = otu_df3[['variant_id', 'marker_id', 'run_id', 'phylum', 'class', 'order', 'family', 'genus',
                           'species', 'ltg_tax_id', 'ltg_tax_name', 'identity', 'ltg_rank', 'chimera_borderline', 'sequence']]

        ################################################################################################################
        #
        # Otu_df: Final merge
        #
        ################################################################################################################

        otu_df_final = otu_df1.merge(otu_df2, on=['variant_id', 'marker_id', 'run_id'])
        otu_df_final = otu_df_final.merge(otu_df3, on=['variant_id', 'marker_id', 'run_id'])
        otu_df_final.drop('marker_id', axis=1, inplace=True)
        otu_df_final.drop('run_id', axis=1, inplace=True)

        return otu_df_final

    def get_chimera_borderline_df(self):
        # #########################################################
        #
        # Get variants that passed the filter
        # Get also chimera borderline information
        #
        # #########################################################
        Logger.instance().debug(
            "file: {}; line: {}; Get variants and sequences that passed the filters".format(__file__, inspect.currentframe().f_lineno,'TaxAssign'))

        filter_codon_stop_model_table = self.filter_codon_stop_model.__table__
        filter_chimera_borderline_model_table = self.filter_chimera_borderline_model.__table__
        variant_model_table = self.variant_model.__table__
        stmt_filter_codon_stop = sqlalchemy.select([
                            filter_codon_stop_model_table.c.variant_id,
                            filter_codon_stop_model_table.c.marker_id,
                            filter_codon_stop_model_table.c.run_id,
                            filter_chimera_borderline_model_table.c.filter_delete, ]) \
            .where(filter_chimera_borderline_model_table.c.variant_id == variant_model_table.c.id) \
            .where(filter_codon_stop_model_table.c.filter_delete == 0).distinct()
        # Select to DataFrame
        variant_to_chimera_borderline_list = []
        with self.engine.connect() as conn:
            for row in conn.execute(stmt_filter_codon_stop).fetchall():
                variant_to_chimera_borderline_list.append(row)
        variant_to_chimera_borderline_df = pandas.DataFrame.from_records(variant_to_chimera_borderline_list,
                                                     columns=['variant_id', 'marker_id', 'run_id', 'chimera_borderline'])
        return variant_to_chimera_borderline_df
