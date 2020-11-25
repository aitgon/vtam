from Bio import SeqIO
from sqlalchemy import select, bindparam, func
from vtam.utils.Logger import Logger
from vtam.utils.FileSampleInformation import FileSampleInformation
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.DataframeVariantReadCountLike import DataframeVariantReadCountLike
from wopmars.models.ToolWrapper import ToolWrapper
import inspect
import os
import pandas
import sqlalchemy
import sys

# Compatible with both pre- and post Biopython 1.78:
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None


class VariantReadCount(ToolWrapper):

    __mapper_args__ = {
        "polymorphic_identity": __module__
    }
    # Input
    # Input file
    __input_file_sortedinfo = "sortedinfo"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_sample = "Sample"
    # Output
    # Output file
    # Output table
    __output_table_variant = 'Variant'
    __output_table_variant_read_count = 'VariantReadCount'

    def specify_input_table(self):
        return [
            VariantReadCount.__input_table_run,
            VariantReadCount.__input_table_marker,
            VariantReadCount.__input_table_sample,
        ]

    def specify_input_file(self):
        return[
            # VariantReadCount.__input_file_sort_reads,
            VariantReadCount.__input_file_sortedinfo,
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
            "global_read_count_cutoff": "int",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        #######################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        #######################################################################

        # Input file
        # sort_reads_tsv = self.input_file(VariantReadCount.__input_file_sort_reads)
        input_file_sortedinfo = self.input_file(
            VariantReadCount.__input_file_sortedinfo)
        #
        # Input table models
        run_model = self.input_table(VariantReadCount.__input_table_run)
        marker_model = self.input_table(VariantReadCount.__input_table_marker)
        sample_model = self.input_table(
            VariantReadCount.__input_table_sample)
        #
        # Output
        # Output table
        variant_model = self.output_table(
            VariantReadCount.__output_table_variant)
        variant_read_count_model = self.output_table(
            VariantReadCount.__output_table_variant_read_count)
        # Options
        read_dir = self.option("read_dir")
        global_read_count_cutoff = self.option(
            "global_read_count_cutoff")

        #######################################################################
        #
        # 1. Read sortedinfo to get run_id, marker_id, sample_id, replicate for current analysis
        # 2. Delete marker_name/run_name/sample/replicate from variant_read_count_model
        # 3. Read tsv file with sorted reads
        # 4. Group by read sequence
        # 5. Delete variants if below global_read_count_cutoff
        # 6. Insert into Variant and DataframeVariantReadCountLike tables
        #
        #######################################################################

        #######################################################################
        #
        # 1. Read sample information to get run_id, marker_id, sample_id, replicate for current analysis
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Read sample information".format(
                __file__, inspect.currentframe().f_lineno))
        sortedinfo_df = pandas.read_csv(input_file_sortedinfo, sep="\t", header=0)
        sample_instance_list = []
        sortedinfo_df.columns = sortedinfo_df.columns.str.lower()

        for row in sortedinfo_df.itertuples():
            Logger.instance().debug(row)
            marker_name = row.marker
            run_name = row.run
            sample_name = row.sample
            replicate = row.replicate
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(
                    run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(
                    marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get sample_id ###########
                stmt_select_sample_id = select([sample_model.__table__.c.id]).where(
                    sample_model.__table__.c.name == sample_name)
                sample_id = conn.execute(
                    stmt_select_sample_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id,
                                             'marker_id': marker_id,
                                             'sample_id': sample_id,
                                             'replicate': replicate})

        #######################################################################
        #
        # 2. Delete marker_name/run_name/sample/replicate from variant_read_count_model
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Delete marker_name/run_name/sample/replicate".format(
                __file__, inspect.currentframe().f_lineno))

        with engine.connect() as conn:
            stmt_del = variant_read_count_model.__table__.delete()
            stmt_del = stmt_del.where(
                variant_read_count_model.__table__.c.run_id == bindparam('run_id'))
            stmt_del = stmt_del.where(
                variant_read_count_model.__table__.c.marker_id == bindparam('marker_id'))
            stmt_del = stmt_del.where(
                variant_read_count_model.__table__.c.sample_id == bindparam('sample_id'))
            stmt_del = stmt_del.where(
                variant_read_count_model.__table__.c.replicate == bindparam('replicate'))
            conn.execute(stmt_del, sample_instance_list)

        #######################################################################
        #
        # 3. Read tsv file with sorted reads
        #
        #######################################################################

        # fasta_info_obj = FastaInformationTSV(input_file_sortedinfo, engine=engine)
        # sample_info_ids_df = fasta_info_obj.get_ids_df()
        sample_info_tsv_obj = FileSampleInformation(
            tsv_path=input_file_sortedinfo)
        sample_info_ids_df = sample_info_tsv_obj.to_identifier_df(engine=engine)

        Logger.instance().debug(
            "file: {}; line: {}; Read demultiplexed FASTA files".format(
                __file__, inspect.currentframe().f_lineno))

        variant_read_count_df = pandas.DataFrame()

        for row in sample_info_ids_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            sample_id = row.sample_id
            replicate = row.replicate
            read_fasta = row.sortedfasta

            Logger.instance().debug(
                "file: {}; line: {}; Read FASTA: {}".format(
                    __file__, inspect.currentframe().f_lineno, read_fasta))

            read_fasta_path = os.path.join(read_dir, read_fasta)

            if os.path.exists(read_fasta_path):

                ####################################################################################
                #
                # Read FASTA
                #
                ####################################################################################

                if generic_dna:  # Biopython <1.78
                    sorted_read_list = [str(seq_record.seq).upper() for
                                        seq_record in
                                        SeqIO.parse(read_fasta_path,
                                                    format="fasta",
                                                    alphabet=generic_dna)]
                else:  # Biopython =>1.78
                    sorted_read_list = [str(seq_record.seq).upper() for
                                        seq_record in
                                        SeqIO.parse(read_fasta_path,
                                                    format="fasta",
                                                    alphabet=generic_dna)]

                variant_read_count_df_sorted_i = pandas.DataFrame(
                    {
                        'run_id': [run_id] *
                        len(sorted_read_list),
                        'marker_id': [marker_id] *
                        len(sorted_read_list),
                        'sample_id': [sample_id] *
                        len(sorted_read_list),
                        'replicate': [replicate] *
                        len(sorted_read_list),
                        'read_sequence': sorted_read_list,
                        'read_count': [1] *
                        len(sorted_read_list)})
                #  Compute read count
                variant_read_count_df_sorted_i = variant_read_count_df_sorted_i.groupby(
                    ['run_id', 'marker_id', 'sample_id', 'replicate', 'read_sequence']).sum().reset_index()

                variant_read_count_df = variant_read_count_df.append(
                    variant_read_count_df_sorted_i)

            else:
                Logger.instance().warning('This file {} doest not exists'.format(read_fasta_path))

        #######################################################################
        #
        # 4. Group by read sequence to variant_read_count with run_id, marker_name, ...
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Group by read sequence".format(
                __file__, inspect.currentframe().f_lineno))
        variant_read_count_df = variant_read_count_df.groupby(
            ['run_id', 'marker_id', 'sample_id', 'replicate', 'read_sequence']) .sum().reset_index()
        variant_read_count_df.rename(
            columns={
                'read_sequence': 'variant_id'},
            inplace=True)
        variant_read_count_df.sort_values(
            by=variant_read_count_df.columns.tolist())

        #######################################################################
        #
        # 5. Remove variants with read count across all run_name, markers, samples and replicates lower than
        # global_read_count_cutoff parameter
        #
        #######################################################################

        variant_read_count_like_df_obj = DataframeVariantReadCountLike(variant_read_count_df)
        Logger.instance().debug(
            "file: {}; line: {}; Remove variants with global read count lower than parameter 'global_read_count_cutoff'".format(
                __file__, inspect.currentframe().f_lineno))
        variant_read_count_df = variant_read_count_like_df_obj.filter_out_below_global_read_count_cutoff(
            global_read_count_cutoff=global_read_count_cutoff)
        variant_read_count_df.rename(
            columns={
                'variant_id': 'variant_sequence'},
            inplace=True)

        #######################################################################
        #
        # 6. Insert into Variant and VariantReadCount tables
        #
        #######################################################################

        Logger.instance().debug("file: {}; line: {}; Insert variants".format(
                __file__, inspect.currentframe().f_lineno))
        variant_read_count_instance_list = []
        variant_read_count_df.sort_values(
            by=['variant_sequence', 'run_id', 'marker_id', 'sample_id', 'replicate'], inplace=True)
        variant_new_set = set()
        variant_new_instance_list = []
        with engine.connect() as conn:
            # Retrieve maximal variant id if possible
            select_variant_id_max = conn.execute(sqlalchemy.select(
                [func.max(variant_model.__table__.c.id)])).first()[0]
            if select_variant_id_max is None:
                select_variant_id_max = 0  # If no variants, then maximal variant id is 0
            for row in variant_read_count_df.itertuples():
                run_id = row.run_id
                marker_id = row.marker_id
                sample_id = row.sample_id
                replicate = row.replicate
                variant_sequence = row.variant_sequence
                read_count = row.read_count
                select_row = conn.execute(sqlalchemy.select([variant_model.__table__.c.id]) .where(
                    variant_model.__table__.c.sequence == variant_sequence)).first()
                if select_row is None:  # variant_sequence IS NOT in the database, so will INSERT it
                    if not (variant_sequence in variant_new_set):
                        variant_id = select_variant_id_max + \
                            len(variant_new_instance_list) + 1
                        variant_new_set.add(variant_sequence)
                        variant_new_instance_list.append(
                            {'id': variant_id, 'sequence': variant_sequence})
                else:  # variant_sequence IS in the database
                    variant_id = select_row[0]
                variant_read_count_instance_list.append(
                    {
                        'run_id': run_id,
                        'marker_id': marker_id,
                        'variant_id': variant_id,
                        'sample_id': sample_id,
                        'replicate': replicate,
                        'read_count': read_count})

        #######################################################################
        #
        # Exit if variant_read_count_instance_list empty
        #
        #######################################################################

        if not len(variant_read_count_instance_list):
            Logger.instance().warning(
                VTAMexception(
                    "No new variants in these samples. Maybe singletons? The analysis will stop here.".format(
                        self.__class__.__name__)))
            sys.exit(0)

        #######################################################################
        #
        # Write variant_read_count table
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {};  Insert variant read count".format(
                __file__, inspect.currentframe().f_lineno))

        with engine.connect() as conn:

            # Insert if there some new variants
            if len(variant_new_instance_list) > 0:
                conn.execute(
                    variant_model.__table__.insert(),
                    variant_new_instance_list)

            # Insert new variant_read_count_instances
            conn.execute(
                variant_read_count_model.__table__.insert(),
                variant_read_count_instance_list)

        #######################################################################
        #
        # Touch output tables, to update modification date
        #
        #######################################################################

        for output_table_i in self.specify_output_table():
            declarative_meta_i = self.output_table(output_table_i)
            obj = session.query(declarative_meta_i).order_by(
                declarative_meta_i.id.desc()).first()
            session.query(declarative_meta_i).filter_by(
                id=obj.id).update({'id': obj.id})
            session.commit()
