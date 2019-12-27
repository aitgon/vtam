import pandas
import sqlalchemy
import sys

from vtam import Logger, VTAMexception


class FastaInformation(object):
    """This class defines various methods to retrieve data from the DB based on the Fasta Information TSV"""

    def __init__(self, fasta_info_tsv, engine, run_model, marker_model, biosample_model):
        """A new instance needs the path to the fasta information path information TSV file as wel as DB information to interact with the DB"""
        self.__engine = engine
        self.__run_model = run_model
        self.__marker_model = marker_model
        self.__biosample_model = biosample_model
        #
        self.df = pandas.read_csv(fasta_info_tsv, sep="\t", header=0, \
                                  names=['tag_fwd_sequence', 'primer_fwd_sequence', 'tag_rev_sequence', 'primer_rev_sequence', 'marker_name', 'biosample_name',\
                                                    'replicate', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta_file_name'])
        #
        # self.fastainfo_instance_list = self.get_fasta_information_record_list()

    def get_sample_information_df(self, add_tag_primer_fasta=False):
        """Based on the Fasta information TSV, returns a list of dictionnaries with run_id, marker_id, biosample_id
        and replicate entries (See return)

        :param tag_primer_fasta_information: Boolean. Default=False. If True, will also return tag and primer sequences and fasta file name
        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate': 1}, {'run_id': 1, ...
        """
        fasta_info_instance_list = []
        for row in self.df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            # replicate_name = row.replicate_name
            replicate = row.replicate
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
                # get replicate ###########
                # stmt_select_replicate = sqlalchemy.select([self.__replicate_model.__table__.c.id]).where(self.__replicate_model.__table__.c.name == replicate_name)
                # replicate = conn.execute(stmt_select_replicate).first()[0]
            fasta_information_obj = {'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id, 'replicate': replicate}
            if add_tag_primer_fasta:
                fasta_information_obj['tag_fwd_sequence'] = row.tag_fwd_sequence
                fasta_information_obj['primer_fwd_sequence'] = row.primer_fwd_sequence
                fasta_information_obj['tag_rev_sequence'] = row.tag_rev_sequence
                fasta_information_obj['primer_rev_sequence'] = row.primer_rev_sequence
                fasta_information_obj['fasta_file_name'] = row.fasta_file_name
            # add this sample_instance ###########
            fasta_info_instance_list.append(fasta_information_obj)
            sample_information_df = pandas.DataFrame.from_records(data=fasta_info_instance_list).drop_duplicates()
        return sample_information_df

    def get_fasta_information_record_list(self, tag_primer_fasta_information=False):
        """Based on the Fasta information TSV, returns a list of dictionnaries with run_id, marker_id, biosample_id
        and replicate entries (See return)

        :param tag_primer_fasta_information: Boolean. Default=False. If True, will also return tag and primer sequences and fasta file name
        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate': 1}, {'run_id': 1, ...
        """
        fasta_info_instance_list = []
        for row in self.df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            # replicate_name = row.replicate_name
            replicate = row.replicate
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
                # get replicate ###########
                # stmt_select_replicate = sqlalchemy.select([self.__replicate_model.__table__.c.id]).where(self.__replicate_model.__table__.c.name == replicate_name)
                # replicate = conn.execute(stmt_select_replicate).first()[0]
            fasta_information_obj = {'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id, 'replicate': replicate}
            if tag_primer_fasta_information:
                fasta_information_obj['tag_fwd_sequence'] = row.tag_fwd_sequence
                fasta_information_obj['primer_fwd_sequence'] = row.primer_fwd_sequence
                fasta_information_obj['tag_rev_sequence'] = row.tag_rev_sequence
                fasta_information_obj['primer_rev_sequence'] = row.primer_rev_sequence
                fasta_information_obj['fasta_file_name'] = row.fasta_file_name
            # add this sample_instance ###########
            fasta_info_instance_list.append(fasta_information_obj)
        return fasta_info_instance_list

    def get_variant_read_count_df(self, variant_read_count_like_model, filter_id=None):
        """Based on the Fasta samples and the variant_read_count_model, returns the variant_read_count_df

        :param variant_read_count_like_model: SQLalchemy models with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count
        :param filter_id:
        :return: DataFrame with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count
        """

        variant_read_count_like_table = variant_read_count_like_model.__table__

        variant_read_count_list = []
        for sample_instance in self.get_fasta_information_record_list():
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate = sample_instance['replicate']
            stmt_select = sqlalchemy.select([variant_read_count_like_table.c.run_id,
                                  variant_read_count_like_table.c.marker_id,
                                  variant_read_count_like_table.c.biosample_id,
                                  variant_read_count_like_table.c.replicate,
                                  variant_read_count_like_table.c.variant_id,
                                  variant_read_count_like_table.c.read_count]).distinct()\
                                    .where(variant_read_count_like_table.c.run_id == run_id)\
                                    .where(variant_read_count_like_table.c.marker_id == marker_id)\
                                    .where(variant_read_count_like_table.c.biosample_id == biosample_id)\
                                    .where(variant_read_count_like_table.c.replicate == replicate)
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
            columns=['run_id', 'marker_id', 'biosample_id', 'replicate', 'variant_id', 'read_count'])

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            Logger.instance().warning(VTAMexception("No variants available for this Filter: {}. "
                                                    "The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)
        return variant_read_count_df

    def get_full_table_df(self, model):
        """Based on the Fasta samples returns all the columns of the table as df

        :param model: SQLalchemy models
        :return: DataFrame with all columns and rows of fasta info
        """

        table = model.__table__

        variant_read_count_list = []
        for sample_instance in self.get_fasta_information_record_list():
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate = sample_instance['replicate']
            stmt_select = table.select().distinct()\
                                    .where(table.c.run_id == run_id)\
                                    .where(table.c.marker_id == marker_id)\
                                    .where(table.c.biosample_id == biosample_id)\
                                    .where(table.c.replicate == replicate)

            with self.__engine.connect() as conn:
                for row in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row)
        #
        table_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=[column.key for column in table.c])

        # Exit if no variants for analysis
        try:
            assert table_df.shape[0] > 0
        except AssertionError:
            Logger.instance().warning(VTAMexception("No variants available for this Filter: {}. "
                                                    "The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)
        return table_df

    def get_variant_df(self, variant_read_count_like_model, variant_model, filter_id=None):
        """Based on the Fasta information TSV and variant_model, returns the variant_df

        :return: DataFrame with columns: index, sequence
        """

        variant_read_count_like_df = self.get_variant_read_count_df(variant_read_count_like_model, filter_id)

        record_list = []
        variant_model_table = variant_model.__table__
        for variant_id in variant_read_count_like_df.variant_id.unique().tolist():
            stmt_select = sqlalchemy.select([variant_model_table.c.sequence]).where(variant_model_table.c.id == variant_id)
            with self.__engine.connect() as conn:
                for sequence in conn.execute(stmt_select).first():
                    record_list.append({'id': variant_id, 'sequence': sequence})
        variant_df = pandas.DataFrame.from_records(record_list, index='id')

        return variant_df

    def get_biosample_df(self, variant_read_count_like_model, filter_id=None):
        """Based on the Fasta information TSV and biosample_model, returns the biosample_df

        :return: DataFrame with columns: index, name
        """

        variant_read_count_like_df = self.get_variant_read_count_df(variant_read_count_like_model, filter_id)

        record_list = []
        biosample_model_table = self.__biosample_model.__table__
        for biosample_id in variant_read_count_like_df.biosample_id.unique().tolist():
            stmt_select = sqlalchemy.select([biosample_model_table.c.name]).where(biosample_model_table.c.id == biosample_id)
            with self.__engine.connect() as conn:
                for biosample_name in conn.execute(stmt_select).first():
                    record_list.append({'id': biosample_id, 'name': biosample_name})

        biosample_df = pandas.DataFrame.from_records(record_list, index='id')

        # Sort biosample names according to fastainfo biosample order
        biosample_name_list = self.df.biosample_name.drop_duplicates(keep='first').tolist()
        biosample_df['name'] = pandas.Categorical(biosample_df['name'], biosample_name_list)
        biosample_df.sort_values(by='name', inplace=True)
        biosample_df.name = biosample_df.name.astype(str)

        return biosample_df

    def get_marker_df(self, variant_read_count_like_model, filter_id=None):
        """Based on the Fasta information TSV and marker_model, returns the marker_df

        :return: DataFrame with columns: index, name
        """

        variant_read_count_like_df = self.get_variant_read_count_df(variant_read_count_like_model, filter_id)

        record_list = []
        marker_model_table = self.__marker_model.__table__
        for marker_id in variant_read_count_like_df.marker_id.unique().tolist():
            stmt_select = sqlalchemy.select([marker_model_table.c.name]).where(marker_model_table.c.id == marker_id)
            with self.__engine.connect() as conn:
                for marker_name in conn.execute(stmt_select).first():
                    record_list.append({'id': marker_id, 'name': marker_name})

        marker_df = pandas.DataFrame.from_records(record_list, index='id')

        return marker_df

    def get_run_df(self, variant_read_count_like_model, filter_id=None):
        """Based on the Fasta information TSV and run_model, returns the run_df

        :return: DataFrame with columns: index, name
        """

        variant_read_count_like_df = self.get_variant_read_count_df(variant_read_count_like_model, filter_id)

        record_list = []
        run_model_table = self.__run_model.__table__
        for run_id in variant_read_count_like_df.run_id.unique().tolist():
            stmt_select = sqlalchemy.select([run_model_table.c.name]).where(run_model_table.c.id == run_id)
            with self.__engine.connect() as conn:
                for run_name in conn.execute(stmt_select).first():
                    record_list.append({'id': run_id, 'name': run_name})

        run_df = pandas.DataFrame.from_records(record_list, index='id')

        return run_df

    def get_variant_to_chimera_borderline_df(self, filter_chimera_borderline_model):
        """Based on the Fasta information TSV and run_model, returns the run_df

        :return: DataFrame with columns: index, name
        """

        # fasta_info = FastaInformation(self.fasta_info_tsv, self.engine, self.run_model, self.marker_model, self.biosample_model)

        filter_chimera_borderline_df = self.get_full_table_df(model=filter_chimera_borderline_model)

        variant_to_chimera_borderline_df = filter_chimera_borderline_df[['run_id', 'marker_id', 'variant_id', 'filter_delete']]
        variant_to_chimera_borderline_df = variant_to_chimera_borderline_df.sort_values(by='filter_delete', ascending=False)
        variant_to_chimera_borderline_df = variant_to_chimera_borderline_df\
            .drop_duplicates(subset=['run_id', 'marker_id', 'variant_id'], keep='first', inplace=False)
        variant_to_chimera_borderline_df.rename({'filter_delete': 'chimera_borderline'}, axis=1, inplace=True)

        return variant_to_chimera_borderline_df

