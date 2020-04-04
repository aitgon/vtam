import inspect
import os

import pandas
import sqlalchemy
from Bio import SeqIO
from sqlalchemy import select, bindparam, func
from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.Logger import Logger
from vtam.utils.SampleInformationUtils import FastaInformationTSV
from vtam.utils.VariantReadCountDF import VariantReadCountDF


# from vtam.utils.FilterLFNrunner import f1_lfn_delete_singleton


class VariantReadCount(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.VariantReadCount"
    }
    # Input
    # Input file
    # __input_file_sort_reads = 'sortreads'
    __input_file_readinfo = "readinfo"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
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
        ]

    def specify_input_file(self):
        return[
            # VariantReadCount.__input_file_sort_reads,
            VariantReadCount.__input_file_readinfo,
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
            "read_dir": "str",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        ################################################################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ################################################################################################################

        # Input file
        # sort_reads_tsv = self.input_file(VariantReadCount.__input_file_sort_reads)
        input_file_readinfo = self.input_file(VariantReadCount.__input_file_readinfo)
        #
        # Input table models
        run_model = self.input_table(VariantReadCount.__input_table_run)
        marker_model = self.input_table(VariantReadCount.__input_table_marker)
        biosample_model = self.input_table(VariantReadCount.__input_table_biosample)
        #
        # Output
        # Output table
        variant_model = self.output_table(VariantReadCount.__output_table_variant)
        variant_read_count_model = self.output_table(VariantReadCount.__output_table_variant_read_count)
        # Options
        read_dir = self.option("read_dir")

        ################################################################################################################
        #
        # 1. Read readinfo to get run_id, marker_id, biosample_id, replicate for current analysis
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        # 3. Read tsv file with sorted reads
        # 4. Group by read sequence
        # 5. Delete singleton
        # 6. Insert into Variant and VariantReadCountDF tables
        #
        ################################################################################################################

        ################################################################################################################
        #
        # 1. Read sample information to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ################################################################################################################

        Logger.instance().debug("file: {}; line: {}; Read sample information".format(__file__, inspect.currentframe().f_lineno))
        readinfo_df = pandas.read_csv(input_file_readinfo, sep="\t", header=0)
        sample_instance_list  = []
        readinfo_df.columns = readinfo_df.columns.str.lower()

        for row in readinfo_df.itertuples():
            Logger.instance().debug(row)
            marker_name = row.marker
            run_name = row.run
            biosample_name = row.biosample
            replicate = row.replicate
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id,
                                             'replicate': replicate})

        ################################################################################################################
        #
        # 2. Remove marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        Logger.instance().debug("file: {}; line: {}; Remove marker/run/biosample/replicate".format(__file__, inspect.currentframe().f_lineno))
        with engine.connect() as conn:
            stmt = variant_read_count_model.__table__.delete()
            stmt = stmt.where(variant_read_count_model.__table__.c.run_id == bindparam('run_id'))
            stmt = stmt.where(variant_read_count_model.__table__.c.marker_id == bindparam('marker_id'))
            stmt = stmt.where(variant_read_count_model.__table__.c.biosample_id == bindparam('biosample_id'))
            stmt = stmt.where(variant_read_count_model.__table__.c.replicate == bindparam('replicate'))
            conn.execute(stmt, sample_instance_list)

        ##########################################################
        #
        # 3. Read tsv file with sorted reads
        #
        ##########################################################

        fasta_info_obj = FastaInformationTSV(input_file_readinfo, engine=engine)
        fasta_info_ids_df = fasta_info_obj.get_ids_df()

        Logger.instance().debug("file: {}; line: {}; Read demultiplexed FASTA files".format(__file__, inspect.currentframe().f_lineno))

        variant_read_count_df = pandas.DataFrame()

        for row in fasta_info_ids_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            replicate = row.replicate
            read_fasta = row.sortedfasta

            Logger.instance().debug(
                "file: {}; line: {}; Read FASTA: {}".format(__file__, inspect.currentframe().f_lineno, read_fasta))

            read_fasta_path = os.path.join(read_dir, read_fasta)

            if os.path.exists(read_fasta_path):
                with open(read_fasta_path, "r") as fin:
                    sorted_read_list = [x.upper() for x in fin.read().split("\n") if (not x.startswith('>') and x != '')]
                variant_read_count_df_sorted_i = pandas.DataFrame({'run_id': [run_id] * len(sorted_read_list),
                                                                   'marker_id': [marker_id] * len(sorted_read_list),
                                                                   'biosample_id': [biosample_id] * len(
                                                                       sorted_read_list),
                                                                   'replicate': [replicate] * len(sorted_read_list),
                                                                   'read_sequence': sorted_read_list,
                                                                   'read_count': [1] * len(sorted_read_list)})
                #  Compute read count
                variant_read_count_df_sorted_i = variant_read_count_df_sorted_i.groupby(
                    ['run_id', 'marker_id', 'biosample_id', 'replicate', 'read_sequence']).sum().reset_index()

                variant_read_count_df = variant_read_count_df.append(variant_read_count_df_sorted_i)

            else:
                Logger.instance().warning('This file {} doest not exists'.format(read_fasta_path))

        ################################################################################################################
        #
        # 4. Group by read sequence to variant_read_count with run_id, marker, ...
        #
        ################################################################################################################

        Logger.instance().debug("file: {}; line: {}; Group by read sequence".format(__file__, inspect.currentframe().f_lineno))
        # variant_read_count_input_df = variant_read_count_input_df.groupby(['run_id', 'marker_id', 'biosample_id', 'replicate',
        #                                                        'read_sequence']).size().reset_index(name='read_count')
        variant_read_count_df = variant_read_count_df.groupby(['run_id', 'marker_id', 'biosample_id', 'replicate', 'read_sequence'])\
            .sum().reset_index()
        variant_read_count_df.rename(columns={'read_sequence': 'variant_id'}, inplace=True)
        variant_read_count_df.sort_values(by=variant_read_count_df.columns.tolist())

        ################################################################################################################
        #
        # 5. Remove variants with absolute read count lower than lfn_read_count_threshold
        #
        ################################################################################################################

        variant_read_count_lfn = VariantReadCountDF(variant_read_count_df)
        Logger.instance().debug("file: {}; line: {}; Remove singletons".format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_df = variant_read_count_lfn.filter_out_singletons() # returns variant_read_count wout singletons
        variant_read_count_df.rename(columns={'variant_id': 'variant_sequence'}, inplace=True)

        ################################################################################################################
        #
        # 6. Insert into Variant and VariantReadCountDF tables
        #
        ################################################################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Insert variants".format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_instance_list = []
        sample_instance_list = []
        variant_read_count_df.sort_values(
            by=['variant_sequence', 'run_id', 'marker_id', 'biosample_id', 'replicate'], inplace=True)
        variant_new_set = set()
        variant_new_instance_list = []
        with engine.connect() as conn:
            select_variant_id_max = conn.execute(sqlalchemy.select([func.max(variant_model.__table__.c.id)])).first()[0]
            if select_variant_id_max is None:
                select_variant_id_max = 0
            for row in variant_read_count_df.itertuples():
                run_id = row.run_id
                marker_id = row.marker_id
                biosample_id = row.biosample_id
                replicate = row.replicate
                variant_sequence = row.variant_sequence
                read_count = row.read_count
                select_row = conn.execute(sqlalchemy.select([variant_model.__table__.c.id])
                                          .where(variant_model.__table__.c.sequence == variant_sequence)).first()
                if select_row is None:  # variant_sequence IS NOT in the database, so will INSERT it
                    if not (variant_sequence in variant_new_set):
                        variant_id = select_variant_id_max + len(variant_new_instance_list) + 1
                        variant_new_set.add(variant_sequence)
                        variant_new_instance_list.append({'id': variant_id,
                                                          'sequence': variant_sequence})

                else:  # variant_sequence IS in the database
                    variant_id = select_row[0]
                variant_read_count_instance_list.append({'run_id': run_id, 'marker_id': marker_id,
                    'variant_id':variant_id, 'biosample_id': biosample_id, 'replicate': replicate, 'read_count': read_count})
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id,
                                             'replicate': replicate})

        ################################################################################################################
        #
        # Write variant_read_count table
        #
        ################################################################################################################

        Logger.instance().debug(
            "file: {}; line: {};  Insert variant read count".format(__file__, inspect.currentframe().f_lineno))
        # import pdb; pdb.set_trace()
        with engine.connect() as conn:
            if len(variant_new_instance_list) > 0:
                conn.execute(variant_model.__table__.insert(), variant_new_instance_list)
            conn.execute(variant_read_count_model.__table__.delete(), sample_instance_list)
            conn.execute(variant_read_count_model.__table__.insert(), variant_read_count_instance_list)

        # Touch variant table to update modification date
        if len(sample_instance_list) > 0:
            with engine.connect() as conn:
                variant_id, variant_sequence = conn.execute(sqlalchemy.select([variant_model.__table__])).first()
                stmt_update = variant_model.__table__.update()\
                    .where(variant_model.__table__.c.id == variant_id).values(sequence=variant_sequence)
                conn.execute(stmt_update)
