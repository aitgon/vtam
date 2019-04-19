import errno
import inspect
import os
import sqlite3

from Bio.Blast import NCBIWWW
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from sqlalchemy import select
import pandas

from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.wrapper.TaxAssignUtilities import f02_variant_df_to_fasta


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.TaxAssign"}

    # Input file
    # Input table
    __input_table_filter_codon_stop = "FilterCodonStop"
    __input_table_Variant = "Variant"

    # Output table
    __output_table_tax_assign = "TaxAssign"

    def specify_input_file(self):
        return[
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
            TaxAssign.__output_table_tax_assign
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
        #
        # Input table models
        filter_codon_stop_model = self.input_table(TaxAssign.__input_table_filter_codon_stop)
        variant_model = self.input_table(TaxAssign.__input_table_Variant)
        # Output table models
        #
        # Input file path
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
        # 3- test_f06_3_annotate_blast_output_with_tax_id & test_f03_import_blast_output_into_df
        #
        ##########################################################
        # path to the  nuclgb accession2taxid db
        nucl_gb_accession2taxid_sqlite = os.path.join(os.getcwd(), "db_accession2taxid.sqlite")

        # import pdb; pdb.set_trace()
        con = sqlite3.connect(nucl_gb_accession2taxid_sqlite)
        sql = """SELECT gb_accession, tax_id FROM nucl_gb_accession2taxid WHERE gb_accession IN {}""".format(
            tuple(blast_variant_result_df.gb_accession.tolist()))
        gb_accession_to_tax_id_df = pandas.read_sql(sql=sql, con=con)
        con.close()
        #
        # final result blast contaning all the info that we gonna need
        blast_result_tax_id_df = blast_variant_result_df.merge(gb_accession_to_tax_id_df, on='gb_accession')

        import pdb;
        pdb.set_trace()

        ##########################################################
        #
        # 4- test_f03_1_tax_id_to_taxonomy_lineage
        ##########################################################







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




