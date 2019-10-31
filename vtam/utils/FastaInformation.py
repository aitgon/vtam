import pandas
import sqlalchemy
import sys

from vtam import Logger, VTAMexception


class FastaInformation(object):
    """This class defines various methods to retrieve data from the DB based on the Fasta Information TSV"""

    def __init__(self, fasta_info_tsv, engine, run_model, marker_model, biosample_model, replicate_model):
        """A new instance needs the path to the fasta information TSV file as wel as DB information to interact with the DB"""
        self.__engine = engine
        self.__run_model = run_model
        self.__marker_model = marker_model
        self.__biosample_model = biosample_model
        self.__replicate_model = replicate_model
        #
        self.__fastainfo_df = pandas.read_csv(fasta_info_tsv, sep="\t", header=0, \
                                              names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
                                                    'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        #
        # self.fastainfo_instance_list = self.get_fasta_info_record_list()

    def get_fasta_info_record_list(self):
        """Based on the Fasta information TSV, returns a list of dictionnaries with run_id, marker_id, biosample_id
        and replicate_id entries (See return)

        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate_id': 1}, {'run_id': 1, ...
        """
        fasta_info_instance_list = []
        for row in self.__fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with self.__engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select([self.__run_model.__table__.c.id]).where(self.__run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([self.__marker_model.__table__.c.id]).where(self.__marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = sqlalchemy.select([self.__biosample_model.__table__.c.id]).where(self.__biosample_model.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = sqlalchemy.select([self.__replicate_model.__table__.c.id]).where(self.__replicate_model.__table__.c.name == replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                fasta_info_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id,
                                             'replicate_id': replicate_id})
        return fasta_info_instance_list

    def get_variant_read_count_df(self, variant_read_count_like_model, filter_id=None):
        """Based on the Fasta samples and the variant_read_count_model, returns the variant_read_count_df

        :param variant_read_count_like_model: SQLalchemy model with columns: run_id, marker_id, biosample_id, replicate_id, variant_id, read_count
        :return: DataFrame with columns: run_id, marker_id, biosample_id, replicate_id, variant_id, read_count
        """

        variant_read_count_like_table = variant_read_count_like_model.__table__

        variant_read_count_list = []
        for sample_instance in self.get_fasta_info_record_list():
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate_id = sample_instance['replicate_id']
            stmt_select = sqlalchemy.select([variant_read_count_like_table.c.run_id,
                                  variant_read_count_like_table.c.marker_id,
                                  variant_read_count_like_table.c.biosample_id,
                                  variant_read_count_like_table.c.replicate_id,
                                  variant_read_count_like_table.c.variant_id,
                                  variant_read_count_like_table.c.read_count]).distinct()\
                                    .where(variant_read_count_like_table.c.run_id == run_id)\
                                    .where(variant_read_count_like_table.c.marker_id == marker_id)\
                                    .where(variant_read_count_like_table.c.biosample_id == biosample_id)\
                                    .where(variant_read_count_like_table.c.replicate_id == replicate_id)
            # Used for filters tables where filter_delete attribute exists
            if 'filter_delete' in [column.key for column in variant_read_count_like_table.columns]:
                stmt_select = stmt_select.where(variant_read_count_like_table.c.filter_delete == 0)
            # used for filter lfn where filter_id = 8 is necessary (do not pass all filters)
            if not filter_id is None:
                stmt_select = stmt_select.where(variant_read_count_like_table.c.filter_id == filter_id)

            with self.__engine.connect() as conn:
                for row in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row)
        #
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count'])

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            Logger.instance().warning(VTAMexception("No variants available for this Filter: {}. "
                                                    "The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)
        return variant_read_count_df

    def get_variant_df(self, variant_read_count_like_model, variant_model):
        """Based on the Fasta information TSV and variant_model, returns the variant_df

        :return: DataFrame with columns: id, sequence
        """

        variant_read_count_like_df = self.get_variant_read_count_df(variant_read_count_like_model)

        variant_record_list = []
        variant_model_table = variant_model.__table__
        for variant_id in variant_read_count_like_df.variant_id.unique().tolist():
            stmt_select = sqlalchemy.select([variant_model_table.c.sequence]).where(variant_model_table.c.id == variant_id)
            with self.__engine.connect() as conn:
                for sequence in conn.execute(stmt_select).fetchone():
                    variant_record_list.append({'id': variant_id, 'sequence': sequence})

        # # Select to DataFrame
        # variant_list = []
        # with self.__engine.connect() as conn:
        #     for row in conn.execute(stmt_variant).fetchall():
        #         variant_list.append(row)
        variant_df = pandas.DataFrame.from_records(variant_record_list)

        return variant_df

