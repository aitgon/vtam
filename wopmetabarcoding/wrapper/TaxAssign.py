import errno
import inspect
import os
import sqlite3

import sqlalchemy
from Bio.Blast import NCBIWWW
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.utils.constants import data_dir
from sqlalchemy import select
import pandas

from wopmetabarcoding.utils import logger, PathFinder
from wopmetabarcoding.utils.constants import tempdir


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.TaxAssign"}

    # Input file
    # Input table
    __input_table_filter_codon_stop = "FilterCodonStop"
    __input_table_Variant = "Variant"

    # Output table


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
        ##########################################################
        #
        # Get variants that passed the filter
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Get variants and sequences that passed the filters".format(__file__, inspect.currentframe().f_lineno))

        filter_codon_stop_model_table = filter_codon_stop_model.__table__
        stmt_variant_filter_lfn = select([filter_codon_stop_model_table.c.marker_id,
                                          filter_codon_stop_model_table.c.run_id,
                                          filter_codon_stop_model_table.c.variant_id,
                                          filter_codon_stop_model_table.c.biosample_id,
                                          filter_codon_stop_model_table.c.replicate_id,
                                          filter_codon_stop_model_table.c.read_count]) \
            .where(filter_codon_stop_model_table.c.filter_id == 14) \
            .where(filter_codon_stop_model_table.c.filter_delete == 1)
        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        codon_stop_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                                                              columns=['marker_id', 'run_id', 'variant_id',
                                                                       'biosample_id', 'replicate_id', 'read_count'])

        # Select id/ sequence from variant table
        variant_model_table = variant_model.__table__
        stmt_variant = select([variant_model_table.c.id,
                               variant_model_table.c.sequence])

        # Select variant table to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                                                   columns=['id', 'sequence'])

        # extract the id and sequence from variant_df that passed based on the id
        variant_passed_df = variant_df.loc[codon_stop_df['variant_id']]

        # creation one fasta file containing all the variant
        #
        ##########################################################
        #
        # Create Fasta from Variants
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Create Fasta from Variants".format(__file__, inspect.currentframe().f_lineno))
        this_tempdir = os.path.join(tempdir, __file__)
        try:
            os.makedirs(this_tempdir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        variant_passed_fasta = os.path.join(this_tempdir, 'variant_passed.fasta')
        ofile = open(variant_passed_fasta, "w")
        for i in range(0, len(variant_df)):
            df = variant_passed_df.iloc[i]
            ofile.write(">" + str(df['id']) + "\n" + str(df['sequence']) + "\n")
        #
        ofile.close()
        #

        ##########################################################
        #
        # 2- Run qblast: test_f06_2_run_qblast
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Run qblast".format(__file__, inspect.currentframe().f_lineno))
        # Run and read blast result
        with open(variant_passed_fasta) as fin:
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
        import pdb;
        pdb.set_trace()
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




