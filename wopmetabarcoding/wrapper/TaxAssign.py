import errno
import inspect
import os
import sqlite3

from Bio.Blast import NCBIWWW
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from sqlalchemy import select
import pandas
from wopmetabarcoding.utils.constants import rank_hierarchy
from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.wrapper.TaxAssignUtilities import f02_variant_df_to_fasta

from wopmetabarcoding.utils.constants import wop_dir

from wopmetabarcoding.wrapper.TaxAssignUtilities import f04_1_tax_id_to_taxonomy_lineage

from wopmetabarcoding.wrapper.TaxAssignUtilities import f06_select_ltg


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.TaxAssign"}

    # Input file
    __input_file_taxonomy = "taxonomy"
    __input_file_accession2taxid = "accession2taxid"
    # Input table
    __input_table_filter_codon_stop = "FilterCodonStop"
    __input_table_Variant = "Variant"

    # Output table
    __output_table_tax_assign = "TaxAssign"

    def specify_input_file(self):
        return[
            TaxAssign.__input_file_taxonomy,
            TaxAssign.__input_file_accession2taxid,
        ]

    def specify_input_table(self):
        return [
            TaxAssign.__input_table_Variant,
            TaxAssign.__input_table_filter_codon_stop,
        ]

    def specify_output_table(self):
        return [
        ]

    def specify_output_table(self):
        return[
            TaxAssign.__output_table_tax_assign,

        ]

    def specify_params(self):
        return {
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind

        ##########################################################
        #
        # 1. Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file path
        taxonomy_sqlite_path = self.input_file(TaxAssign.__input_file_taxonomy)
        accession2taxid_sqlite_path = self.input_file(TaxAssign.__input_file_accession2taxid)
        #
        # Input table models
        filter_codon_stop_model = self.input_table(TaxAssign.__input_table_filter_codon_stop)
        variant_model = self.input_table(TaxAssign.__input_table_Variant)
        # Output table models
        tax_assign_model = self.output_table(TaxAssign.__output_table_tax_assign)
        #
        # Options
        identity_threshold = 97 # percentage
        include_prop = 90 # percentage
        min_number_of_taxa = 3 # count
        #
        ##########################################################
        #
        # Get variants that passed the filter
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Get variants and sequences that passed the filters".format(__file__, inspect.currentframe().f_lineno,'TaxAssign'))

        filter_codon_stop_model_table = filter_codon_stop_model.__table__
        variant_model_table = variant_model.__table__
        stmt_variant = select([filter_codon_stop_model_table.c.variant_id, variant_model_table.c.sequence]) \
            .where(filter_codon_stop_model_table.c.variant_id == variant_model_table.c.id) \
            .where(filter_codon_stop_model_table.c.filter_delete == 0).distinct().order_by("variant_id")
        # Select to DataFrame
        variant_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_list.append(row)
        variant_df = pandas.DataFrame.from_records(variant_list, columns=['variant_id', 'variant_sequence'])

        # creation one fasta file containing all the variant
        #
        ##########################################################
        #
        # 2 Create Fasta from Variants
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Create Fasta from Variants".format(__file__, inspect.currentframe().f_lineno ,'TaxAssign'))
        this_tempdir = os.path.join(tempdir, os.path.basename(__file__))
        try:
            os.makedirs(this_tempdir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        variant_fasta = os.path.join(this_tempdir, 'variant.fasta')
        f02_variant_df_to_fasta(variant_df, variant_fasta)
        #

        ##########################################################
        #
        # 3 Run qblast: test_f06_2_run_qblast
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Run qblast".format(__file__, inspect.currentframe().f_lineno,'TaxAssign'))
        # Run and read blast result
        with open(variant_fasta) as fin:
            result_handle = NCBIWWW.qblast("blastn", "nt", fin.read(), format_type='Tabular')
            blast_result_tsv = os.path.join(this_tempdir, "tax_assign_blast.tsv")
            with open(blast_result_tsv, 'w') as out_handle:
                out_handle.write(result_handle.read())
            result_handle.close()
        blast_result_df = pandas.read_csv(blast_result_tsv, sep="\t", skiprows=13, usecols=[0, 1, 2],
                                          header=None, names=['variant_id', 'gb_accession', 'identity'])

        # let only the output coantaining non null values
        blast_variant_result_df = blast_result_df[blast_result_df['identity'].notnull()].copy()

        ##########################################################
        #
        # 4 test_f06_3_annotate_blast_output_with_tax_id & test_f03_import_blast_output_into_df
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Annotation blast output".format(__file__, inspect.currentframe().f_lineno, 'TaxAssign'))

        # import pdb; pdb.set_trace()
        con = sqlite3.connect(accession2taxid_sqlite_path)
        sql = """SELECT gb_accession, tax_id FROM nucl_gb_accession2taxid WHERE gb_accession IN {}""".format(tuple(blast_variant_result_df.gb_accession.tolist()))
        gb_accession_to_tax_id_df = pandas.read_sql(sql=sql, con=con)
        con.close()
        #
        # final result blast contaning all the info that we gonna need
        blast_result_tax_id_df = blast_variant_result_df.merge(gb_accession_to_tax_id_df, on='gb_accession')

        ##########################################################
        #
        # 5 test_f03_1_tax_id_to_taxonomy_lineage for many gb accessions and many variant_ids
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Taxonomy Lineage creation".format(__file__, inspect.currentframe().f_lineno,'TaxAssign'))
        # getting the taxonomy_db to df
        con = sqlite3.connect(taxonomy_sqlite_path)
        sql = """SELECT *  FROM taxonomy """
        taxonomy_db_df = pandas.read_sql(sql=sql, con=con)
        con.close()
        #lineage from tax_id
        lineage_list = []
        for target_tax_id in blast_result_tax_id_df.tax_id.unique().tolist():
            lineage_list.append(f04_1_tax_id_to_taxonomy_lineage(target_tax_id,taxonomy_db_df))
        #lineage to df
        tax_lineage_df = pandas.DataFrame(lineage_list)
        tax_lineage_df = blast_result_tax_id_df.merge(tax_lineage_df, left_on='tax_id', right_on='tax_id')

        ##########################################################
        #
        #  6 test_f05_select_ltg_identity
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Ltg,".format(__file__, inspect.currentframe().f_lineno,  'TaxAssign'))
        variant_id_list = tax_lineage_df.variant_id.unique().tolist()
        #
        identity_list = [100, 99, 97, 95, 90, 85, 80, 75, 70]
        variant_id = variant_id_list[0]
        list_ltg = []
        #
        for variant_id in variant_id_list:
            for identity in identity_list:
                tax_lineage_by_variant_id_df = tax_lineage_df.loc[((tax_lineage_df['variant_id'] == str(variant_id)) & (tax_lineage_df['identity'] >= identity))].copy()

                if tax_lineage_by_variant_id_df.shape[0] > 0:

                    if identity < identity_threshold: # Case 2: min_number_of_taxa

                        if tax_lineage_by_variant_id_df.shape[0] >= min_number_of_taxa: # More/equal rows than min_number_of_taxa
                            # select column that are intersect with the rank_hirarchy
                            lineage_list_df_columns_sorted = [value for value in tax_lineage_by_variant_id_df if value in rank_hierarchy]
                            tax_lineage_by_variant_id_df =tax_lineage_by_variant_id_df[lineage_list_df_columns_sorted]
                            # drop column with NaN
                            tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df.dropna(axis='columns', how='all')
                            #
                            ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_by_variant_id_df, identity, identity_threshold=identity_threshold,
                                                                  include_prop=include_prop, min_number_of_taxa=min_number_of_taxa)
                            lineage_dic = {}
                            lineage_dic['variant_id'] = variant_id
                            lineage_dic['identity'] = identity
                            lineage_dic['ltg_tax_id'] = ltg_tax_id
                            lineage_dic['ltg_rank'] = ltg_rank
                            lineage_dic['ltg_lineage'] = None  #  Todo: How?
                            # dictionnary to list
                            list_ltg.append(lineage_dic)
                    else: # Case 1: include_prop
                        # select column that are intersect with the rank_hirarchy
                        lineage_list_df_columns_sorted = [value for value in tax_lineage_by_variant_id_df if
                                                          value in rank_hierarchy]
                        tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df[lineage_list_df_columns_sorted]
                        # drop column with NaN
                        tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df.dropna(axis='columns', how='all')
                        #
                        ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_by_variant_id_df, identity, identity_threshold=identity_threshold,
                                                              include_prop=include_prop, min_number_of_taxa=min_number_of_taxa)
                        lineage_dic = {}
                        lineage_dic['variant_id'] = variant_id
                        lineage_dic['identity'] = identity
                        lineage_dic['ltg_tax_id'] = ltg_tax_id
                        lineage_dic['ltg_rank'] = ltg_rank
                        lineage_dic['ltg_lineage'] = None  #  Todo: How?
                        # dictionnary to list
                        list_ltg.append(lineage_dic)

        tax_lineage_vars_df = pandas.DataFrame(list_ltg)

        #let only the variant with the higher identity
        tax_lineage_df = tax_lineage_vars_df.groupby("variant_id").max()
        tax_lineage_df = tax_lineage_df.reset_index()
        import pdb;
        pdb.set_trace()

        tax_lineage_df.to_sql(name='TaxAssig', con=engine.connect(), if_exists='replace')
        #

            ##########################################################
        #
        # Previous preparation -------------------
        # bin/create_db_taxonomy (sqlite) / done
        # bin/create_db_accession_to_tax_id (sqlite) / done  (do result sqlite file  need to stell in tmp directory?)
        #
        # This wrapper -------------------
        # 1 Create Fasta from Variants that passed the Filters: One FASTA for all variants
        # 2 Run qblast: test_f06_2_run_qblast
        # 3 test_f06_3_annotate_blast_output_with_tax_id & test_f03_import_blast_output_into_df
        # 4 test_f03_1_tax_id_to_taxonomy_lineage
        # 5 test_f05_select_ltg_identity
        #
        ##########################################################




