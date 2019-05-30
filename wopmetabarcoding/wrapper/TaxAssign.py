import errno
import inspect
import os
import sqlite3

from Bio.Blast.Applications import NcbiblastnCommandline
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from sqlalchemy import select
import pandas
from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.utils.constants import tempdir, download_coi_db, download_taxonomy_sqlite
from wopmetabarcoding.wrapper.TaxAssignUtilities import f02_variant_df_to_fasta, f01_taxonomy_sqlite_to_df

from wopmetabarcoding.wrapper.TaxAssignUtilities import f04_1_tax_id_to_taxonomy_lineage

from wopmetabarcoding.wrapper.TaxAssignUtilities import f07_blast_result_to_ltg_tax_id


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.TaxAssign"}

    # Input file
    __input_file_sample2fasta = "sample2fasta"
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


        #################        taxonomy_sqlite_path = download_taxonomy_sqlite()
#########################################
        #
        # 2. Read sample2f  asta to get run_id, marker_id, biosample_id, replicate_id for current analysis
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
        # 2 Create FASTA file with Variants
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
        # 3 Run local blast
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Running local blast with FASTA input {}".format(__file__, inspect.currentframe().f_lineno, variant_fasta))
        #
        # Run and read local blast result
        blast_output_tsv = os.path.join(this_tempdir, 'blast_output.tsv')
        map_taxids_tsv_path, coi_blast_db_dir = download_coi_db()
        os.environ['BLASTDB'] = coi_blast_db_dir
        #
        # Uncomment to run blast with full NCBI blast db)
        #
        # blastn_cline = NcbiblastnCommandline(query=variant_fasta, db="nt", evalue=1e-5,
        #                                      outfmt='"6 qseqid sacc pident evalue qcovhsp staxids"', dust='yes',
        #                                      qcov_hsp_perc=80, num_threads=1, out=blast_output_tsv)
        #
        # Comment to run blast with full NCBI blast db)
        #
        blastn_cline = NcbiblastnCommandline(query=variant_fasta, db="coi_db", evalue=1e-5,
                                             outfmt='"6 qseqid sacc pident evalue qcovhsp"', dust='yes',
                                             qcov_hsp_perc=80, num_threads=1, out=blast_output_tsv)
        stdout, stderr = blastn_cline()
        # import pdb; pdb.set_trace()
        #

        ##########################################################
        #
        # 3 Merge blast result with tax ids (Based on COI blast db)
        #
        ##########################################################
        blast_result_df = pandas.read_csv(blast_output_tsv, header=None, sep="\t",
                                          names=['variant_id', 'target_id', 'identity', 'evalue', 'coverage'])
        map_tax_id_df = pandas.read_csv(map_taxids_tsv_path, header=None, sep="\t", names=['target_id', 'target_tax_id'])
        lblast_output_df = blast_result_df.merge(map_tax_id_df, on='target_id')
        #
        ##########################################################
        #
        # Process local blast TSV result (Uncomment to run blast with full NCBI blast db)
        #
        ##########################################################
        # blast_output_tsv = "/home/gonzalez/tmp/vtam/tmpr95f12vh/TaxAssign.py/blast_output.tsv"
        # logger.debug(
        #     "file: {}; line: {}; Reading TSV output from local blast: {}".format(__file__, inspect.currentframe().f_lineno, blast_output_tsv))
        # lblast_output_df = pandas.read_csv(blast_output_tsv, sep='\t', header=None,
        #                                   names=['variant_id', 'target_id', 'identity', 'evalue', 'coverage',
        #                                          'target_tax_id'])
        #
        # remove multiple target_tax_ids, rename and convert target_tax_id to numeric
        # lblast_output_df = (pandas.concat([lblast_output_df, lblast_output_df.target_tax_id.str.split(pat=';', n=1, expand=True)],axis=1))[['variant_id', 'identity', 0]]
        # lblast_output_df = lblast_output_df.rename(columns={0: 'target_tax_id'})
        # import pdb; pdb.set_trace()
        #
        ##########################################################
        #
        # Read target_tax_id
        # Compute lineages for each unique target_tax_id
        # Create a DF with these columns: tax_id and its lineage in wide format
        # Merge to the blast result
        #
        ##########################################################
        #
        logger.debug(
            "file: {}; line: {}; Open taxonomy.sqlite DB".format(__file__, inspect.currentframe().f_lineno))
        lblast_output_df.target_tax_id = pandas.to_numeric(lblast_output_df.target_tax_id)
        # getting the taxonomy_db to df
        taxonomy_sqlite_path = download_taxonomy_sqlite()
        taxonomy_db_df = f01_taxonomy_sqlite_to_df(taxonomy_sqlite_path)
        #
        logger.debug(
            "file: {}; line: {}; Annotate each target_tax_id with its lineage as columns in wide format".format(__file__, inspect.currentframe().f_lineno))
        lineage_list = []
        for target_tax_id in lblast_output_df.target_tax_id.unique().tolist():
            lineage_list.append(f04_1_tax_id_to_taxonomy_lineage(target_tax_id, taxonomy_db_df))
        tax_id_to_lineage_df = pandas.DataFrame(lineage_list)
        #
        logger.debug(
            "file: {}; line: {}; Merge blast result including tax_id with their lineages".format(__file__, inspect.currentframe().f_lineno))
        # Merge local blast output with tax_id_to_lineage_df
        variantid_identity_lineage_df = lblast_output_df.merge(tax_id_to_lineage_df, left_on='target_tax_id',
                                                               right_on='tax_id')
        variantid_identity_lineage_df.drop('tax_id', axis=1, inplace=True)
        variantid_identity_lineage_tsv = os.path.join(this_tempdir, 'variantid_identity_lineage.tsv')
        variantid_identity_lineage_df.to_csv(variantid_identity_lineage_tsv, sep="\t", header=True)

        ##########################################################
        #
        #  6 test_f05_select_ltg_identity
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Main loop over variant and identity to"
            "compute the whole set of ltg_tax_id and ltg_rank for each variant_id"
            "to a dataframe".format(__file__, inspect.currentframe().f_lineno))
        #
        # f07_blast_result_to_ltg_tax_id(tax_lineage_df,identity_threshold=97, include_prop=90, min_number_of_taxa=3):
        # this function return a data frame containing the Ltg rank and Ltg Tax_id for each variant
        #
        ltg_df = f07_blast_result_to_ltg_tax_id(variantid_identity_lineage_df, identity_threshold=identity_threshold,
                                                include_prop=include_prop, min_number_of_taxa=min_number_of_taxa)
        ##########################################################
        #
        # 6. Insert Filter data
        #
        ##########################################################
        logger.debug("file: {}; line: {}; Insert variant_id, ltg_tax_id, ltg_rank to DB".format(__file__, inspect.currentframe().f_lineno))
        with engine.connect() as conn:
                conn.execute(tax_assign_model.__table__.insert(), ltg_df.to_dict('records'))






