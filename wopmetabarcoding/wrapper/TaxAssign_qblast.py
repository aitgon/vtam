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


from wopmetabarcoding.wrapper.TaxAssignUtilities import f04_1_tax_id_to_taxonomy_lineage

from wopmetabarcoding.wrapper.TaxAssignUtilities import f07_blast_result_to_ltg_tax_id


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.TaxAssign"}

    # Input file
    __input_file_sample2fasta = "sample2fasta"
    __input_file_taxonomy = "taxonomy"
    __input_file_accession2taxid = "accession2taxid"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_filter_codon_stop = "FilterCodonStop"
    __input_table_variant = "Variant"

    # Output table
    __output_table_tax_assign = "TaxAssign"

    def specify_input_file(self):
        return[
            TaxAssign.__input_file_sample2fasta,
            TaxAssign.__input_file_taxonomy,
            TaxAssign.__input_file_accession2taxid,
        ]

    def specify_input_table(self):
        return [
            TaxAssign.__input_table_marker,
            TaxAssign.__input_table_run,
            TaxAssign.__input_table_biosample,
            TaxAssign.__input_table_replicate,
            TaxAssign.__input_table_variant,
            TaxAssign.__input_table_filter_codon_stop,
        ]

    def specify_output_table(self):
        return[
            TaxAssign.__output_table_tax_assign,

        ]

    def specify_params(self):
        return {
        "identity_threshold":"float",  #  percentage
        "include_prop" : "float",  # percentage
        "min_number_of_taxa" : "int",  #  count
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
        input_file_sample2fasta = self.input_file(TaxAssign.__input_file_sample2fasta)
        #
        # Input table models
        marker_model = self.input_table(TaxAssign.__input_table_marker)
        run_model = self.input_table(TaxAssign.__input_table_run)
        biosample_model = self.input_table(TaxAssign.__input_table_biosample)
        replicate_model = self.input_table(TaxAssign.__input_table_replicate)
        filter_codon_stop_model = self.input_table(TaxAssign.__input_table_filter_codon_stop)
        variant_model = self.input_table(TaxAssign.__input_table_variant)
        # Output table models
        tax_assign_model = self.output_table(TaxAssign.__output_table_tax_assign)
        #
        # Options
        identity_threshold = float(self.option("identity_threshold")) # percentage
        include_prop = float(self.option("include_prop")) # percentage
        min_number_of_taxa = int(self.option("min_number_of_taxa")) # count


        ##########################################################
        #
        # 2. Read sample2fasta to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        sample2fasta_df = pandas.read_csv(input_file_sample2fasta, sep="\t", header=None,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
            'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        sample_instance_list = []
        for row in sample2fasta_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([replicate_model.__table__.c.id]).where(replicate_model.__table__.c.name==replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id, 'replicate_id':replicate_id})

        ##########################################################
        #
        # 3a. Select variants from this run/markerbiosample/replicate combination
        # 3b. Delete variants from this run/markerbiosample/replicate combination
        #
        ##########################################################
        #
        # 3a. Select variants from this run/markerbiosample/replicate combination
        #
        variant_id_delete_list = []
        for sample_instance in sample_instance_list:
            stmt_select_variant_id_delete = select([filter_codon_stop_model.__table__.c.variant_id])\
                .where(filter_codon_stop_model.__table__.c.run_id == sample_instance['run_id'])\
                .where(filter_codon_stop_model.__table__.c.marker_id == sample_instance['marker_id'])\
                .where(filter_codon_stop_model.__table__.c.biosample_id == sample_instance['biosample_id'])\
                .where(filter_codon_stop_model.__table__.c.replicate_id == sample_instance['replicate_id'])
            # Select to DataFrame
            with engine.connect() as conn:
                for row in conn.execute(stmt_select_variant_id_delete).fetchall():
                    variant_id = row[0]
                    if not variant_id in variant_id_delete_list:
                        variant_id_delete_list.append(row)
        #
        # 3b. Delete variants from this run/markerbiosample/replicate combination
        #
        variant_instance_list = [{'variant_id': variant_id} for variant_id in variant_id_delete_list]
        with engine.connect() as conn:
            conn.execute(tax_assign_model.__table__.delete(), variant_instance_list)

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
            "file: {}; line: {}; Create Fasta from Variants".format(__file__, inspect.currentframe().f_lineno))
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
            "file: {}; line: {}; Running qblast with FASTA input {}".format(__file__, inspect.currentframe().f_lineno, variant_fasta))
        # Run and read qblast result
        with open(variant_fasta) as fin:
            result_handle = NCBIWWW.qblast("qblastn", "nt", fin.read(), format_type='Tabular')
            qblast_result_tsv = os.path.join(this_tempdir, "tax_assign_qblast.tsv")
            with open(qblast_result_tsv, 'w') as out_handle:
                out_handle.write(result_handle.read())
            result_handle.close()
        logger.debug(
            "file: {}; line: {}; TSV output from qblast: {}".format(__file__, inspect.currentframe().f_lineno, qblast_result_tsv))
        logger.debug(
            "file: {}; line: {}; Finished qblast".format(__file__, inspect.currentframe().f_lineno, ))
        qblast_result_df = pandas.read_csv(qblast_result_tsv, sep="\t", skiprows=13, usecols=[0, 1, 2],
                                          header=None, names=['variant_id', 'gb_accession', 'identity'])
        # let only the output coantaining non null values
        qblast_variant_result_df = qblast_result_df[qblast_result_df['identity'].notnull()].copy()

        ##########################################################
        #
        # 4 test_f06_3_annotate_qblast_output_with_tax_id & test_f03_import_qblast_output_into_df
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Annotation qblast output".format(__file__, inspect.currentframe().f_lineno))
        logger.debug(
            "file: {}; line: {}; Connect to accession2taxid_sqlite: {}".format(__file__, inspect.currentframe().f_lineno,
                                                                               accession2taxid_sqlite_path))
        con = sqlite3.connect(accession2taxid_sqlite_path)
        sql = """SELECT gb_accession, tax_id FROM nucl_gb_accession2taxid WHERE gb_accession IN {}""".format(tuple(qblast_variant_result_df.gb_accession.tolist()))
        gb_accession_to_tax_id_df = pandas.read_sql(sql=sql, con=con)
        con.close()
        #
        # final result qblast contaning all the info that we gonna need
        qblast_result_tax_id_df = qblast_variant_result_df.merge(gb_accession_to_tax_id_df, on='gb_accession')

        ##########################################################
        #
        # 5 test_f03_1_tax_id_to_taxonomy_lineage for many gb accessions and many variant_ids
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Taxonomy Lineage creation".format(__file__, inspect.currentframe().f_lineno))
        logger.debug(
            "file: {}; line: {}; Connect to taxonomy_sqlite: {}".format(__file__, inspect.currentframe().f_lineno,
                                                                               taxonomy_sqlite_path))
        # getting the taxonomy_db to df
        con = sqlite3.connect(taxonomy_sqlite_path)
        sql = """SELECT *  FROM taxonomy """
        taxonomy_db_df = pandas.read_sql(sql=sql, con=con)
        con.close()
        lineage_list = []
        for target_tax_id in qblast_result_tax_id_df.tax_id.unique().tolist():
            lineage_list.append(f04_1_tax_id_to_taxonomy_lineage(target_tax_id,taxonomy_db_df))
        #lineage to df
        tax_lineage_df = pandas.DataFrame(lineage_list)
        tax_lineage_df = qblast_result_tax_id_df.merge(tax_lineage_df, left_on='tax_id', right_on='tax_id')

        ##########################################################
        #
        #  6 test_f05_select_ltg_identity
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Ltg,".format(__file__, inspect.currentframe().f_lineno))
        #
        # f07_blast_result_to_ltg_tax_id(tax_lineage_df,identity_threshold=97, include_prop=90, min_number_of_taxa=3):
        # this function return a data frame containing the Ltg rank and Ltg Tax_id for each variant
        ltg_df = f07_blast_result_to_ltg_tax_id(tax_lineage_df, identity_threshold=identity_threshold, include_prop=include_prop, min_number_of_taxa=min_number_of_taxa)

        ##########################################################
        #
        # 6. Insert Filter data
        #
        ##########################################################
        with engine.connect() as conn:
                conn.execute(tax_assign_model.__table__.insert(), ltg_df.to_dict('records'))







