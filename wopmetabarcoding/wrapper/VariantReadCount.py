import inspect

import pandas
import sqlalchemy
from sqlalchemy import select
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.wrapper.FilterLFNutilities import f1_lfn_delete_singleton


class VariantReadCount(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.VariantReadCount"
    }
    # Input
    # Input file
    __input_file_sort_reads = 'sortreads'
    __input_file_sample2fasta = "sample2fasta"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    # Output
    # Output file
    # Output table
    __output_table_variant = 'Variant'
    __output_table_variant_read_count = 'VariantReadCount'

    def specify_input_table(self):
        return [
            VariantReadCount.__input_table_run,
            VariantReadCount.__input_table_marker,
            VariantReadCount.__input_table_biosample,
            VariantReadCount.__input_table_replicate,
        ]

    def specify_input_file(self):
        return[
            VariantReadCount.__input_file_sort_reads,
            VariantReadCount.__input_file_sample2fasta,
        ]

    def specify_output_table(self):
        return [
            VariantReadCount.__output_table_variant,
            VariantReadCount.__output_table_variant_read_count,
        ]

    def specify_params(self):
        """

        :return:
        """
        return {
            "min_id": "float",
            "minseqlength": "int",
            "overhang": "int"
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # conn = engine.connect()
        #
        # Input file
        sort_reads_tsv = self.input_file(VariantReadCount.__input_file_sort_reads)
        input_file_sample2fasta = self.input_file(VariantReadCount.__input_file_sample2fasta)
        #
        # Input table models
        run_model = self.input_table(VariantReadCount.__input_table_run)
        marker_model = self.input_table(VariantReadCount.__input_table_marker)
        biosample_model = self.input_table(VariantReadCount.__input_table_biosample)
        replicate_model = self.input_table(VariantReadCount.__input_table_replicate)
        #
        # Output
        # Output table
        variant_model = self.output_table(VariantReadCount.__output_table_variant)
        variant_read_count_model = self.output_table(VariantReadCount.__output_table_variant_read_count)

        ################################
        #
        # 1. Read sample2fasta to get run_id, marker_id, biosample_id, replicate_id for current analysis
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        # 3. Read tsv file with sorted reads
        # 4. Group by read sequence
        # 5. Delete singleton
        # 6. Insert into Variant and VariantReadCount tables
        #
        ################################

        ##########################################################
        #
        # 1. Read sample information to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        logger.debug("file: {}; line: {}; Read sample information".format(__file__, inspect.currentframe().f_lineno))
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
        # 2. Remove marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################
        logger.debug("file: {}; line: {}; Remove marker/run/biosample/replicate".format(__file__, inspect.currentframe().f_lineno))
        with engine.connect() as conn:
            conn.execute(variant_read_count_model.__table__.delete(), sample_instance_list)

        ##########################################################
        #
        # 3. Read tsv file with sorted reads
        #
        ##########################################################
        logger.debug("file: {}; line: {}; Read tsv file with sorted reads".format(__file__, inspect.currentframe().f_lineno))
        read_annotation_df = pandas.read_csv(sort_reads_tsv, sep='\t',
                             header=None,
                             names=['read_id', 'fasta_id', 'run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_sequence'])
        read_annotation_df.drop('read_id', axis=1, inplace=True)
        ##########################################################
        #
        # 4. Group by read sequence
        #
        ##########################################################
        logger.debug("file: {}; line: {}; Group by read sequence".format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_df = read_annotation_df.groupby(['fasta_id', 'run_id', 'marker_id', 'biosample_id',
                                    'replicate_id', 'variant_sequence']).size().reset_index(name='read_count')

        ##########################################################
        #
        # 5. Remove singleton
        #
        ##########################################################
        logger.debug("file: {}; line: {}; Remove singletons".format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_df = f1_lfn_delete_singleton(variant_read_count_df)

        ################################
        #
        # 6. Insert into Variant and VariantReadCount tables
        #
        ################################
        logger.debug(
            "file: {}; line: {}; Insert variants".format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_instance_list = []
        sample_instance_list = []
        for row in variant_read_count_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            replicate_id = row.replicate_id
            variant_sequence = row.variant_sequence
            read_count = row.read_count
            try:
                stmt_ins_var = variant_model.__table__.insert().values(sequence=variant_sequence)
                with engine.connect() as conn:
                    stmt_result_var = conn.execute(stmt_ins_var)
                variant_id = stmt_result_var.inserted_primary_key[0]
            except sqlalchemy.exc.IntegrityError:
                stmt_select_var = select([variant_model.__table__.c.id]).where(variant_model.__table__.c.sequence==variant_sequence)
                with engine.connect() as conn:
                    variant_id = conn.execute(stmt_select_var).first()[0]
            variant_read_count_instance_list.append({'run_id': run_id, 'marker_id': marker_id,
                'variant_id':variant_id, 'biosample_id':biosample_id, 'replicate_id':replicate_id, 'read_count':read_count})
            sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id, 'replicate_id':replicate_id})
            #
        ############################################
        # Write variant_read_count table
        ############################################
        logger.debug(
            "file: {}; line: {};  Insert variant read count".format(__file__, inspect.currentframe().f_lineno))
        with engine.connect() as conn:
            conn.execute(variant_read_count_model.__table__.delete(), sample_instance_list)
            conn.execute(variant_read_count_model.__table__.insert(), variant_read_count_instance_list)
