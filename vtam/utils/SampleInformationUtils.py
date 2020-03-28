import pandas
import sqlalchemy
import sys

from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception

from vtam.models.Run import Run as run_model
from vtam.models.Marker import Marker as marker_model
from vtam.models.Biosample import Biosample as biosample_model

class SampleInformationUtils(object):
    """Takes an object as in the SampleInformation model and carries out several computations on and the DB"""

    def __init__(self, engine, sample_information_df):
        """Takes the engine and the sample_information_df to carry out several computations"""
        self.__engine = engine
        self.sample_information_df = sample_information_df
        self.sample_record_list = list(self.sample_information_df.T.to_dict().values())

    def get_variant_read_count_df(self, variant_read_count_like_model, filter_id=None):
        """Based on the SortedReadFile samples and the variant_read_count_model, returns the variant_read_count_df

        :param variant_read_count_like_model: SQLalchemy models with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count
        :param filter_id:
        :return: DataFrame with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count
        """

        variant_read_count_like_table = variant_read_count_like_model.__table__

        variant_read_count_list = []
        # for sample_instance in self.get_fasta_information_record_list():
        for sample_instance_row in self.sample_information_df.itertuples():
            run_id = sample_instance_row.run_id
            marker_id = sample_instance_row.marker_id
            biosample_id = sample_instance_row.biosample_id
            replicate = sample_instance_row.replicate
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
            Logger.instance().warning(VTAMexception("No variants available after the Codon Stop filter. "
                                                    "The pipeline will stop here.".format(self.__class__.__name__)))
            sys.exit(0)
        return variant_read_count_df

    def get_full_table_df(self, model):
        """Based on the SortedReadFile samples returns all the columns of the table as variant_read_count_df

        :param model: SQLalchemy models
        :return: DataFrame with all columns and rows of fasta info
        """

        table = model.__table__

        variant_read_count_list = []
        # for sample_instance in self.get_fasta_information_record_list():
        for sample_instance_row in self.sample_information_df.itertuples():
            run_id = sample_instance_row.run_id
            marker_id = sample_instance_row.marker_id
            biosample_id = sample_instance_row.biosample_id
            replicate = sample_instance_row.replicate
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
        """Based on the SortedReadFile information TSV and variant_model, returns the variant_df

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

    def get_biosample_df(self, biosample_model):
        """Returns biosample_df for this sample_information_id

        :return: DataFrame with columns: index, name
        """

        record_list = []
        for biosample_id in self.sample_information_df.biosample_id.unique().tolist():
            stmt_select = sqlalchemy.select([biosample_model.__table__.c.name]).where(biosample_model.__table__.c.id == biosample_id)
            with self.__engine.connect() as conn:
                for biosample_name in conn.execute(stmt_select).first():
                    record_list.append({'id': biosample_id, 'name': biosample_name})

        biosample_df = pandas.DataFrame.from_records(record_list, index='id')

        # Sort biosample names according to fastainfo biosample order
        biosample_id_list = self.sample_information_df.biosample_id.drop_duplicates(keep='first').tolist()
        biosample_df.index = pandas.Categorical(biosample_df.index, biosample_id_list)
        biosample_df.sort_index(inplace=True)
        biosample_df.index = biosample_df.index.astype(int)

        return biosample_df

    def get_marker_df(self, marker_model):
        """Returns marker_df for this sample_information_id

        :return: DataFrame with columns: index, name
        """

        record_list = []
        for marker_id in self.sample_information_df.marker_id.unique().tolist():
            stmt_select = sqlalchemy.select([marker_model.__table__.c.name]).where(marker_model.__table__.c.id == marker_id)
            with self.__engine.connect() as conn:
                for marker_name in conn.execute(stmt_select).first():
                    record_list.append({'id': marker_id, 'name': marker_name})

        marker_df = pandas.DataFrame.from_records(record_list, index='id')

        return marker_df

    def get_run_df(self, run_model):
        """Based on the SortedReadFile information TSV and run_model, returns the run_df

        :return: DataFrame with columns: index, name
        """

        record_list = []
        for run_id in self.sample_information_df.run_id.unique().tolist():
            stmt_select = sqlalchemy.select([run_model.__table__.c.name]).where(run_model.__table__.c.id == run_id)
            with self.__engine.connect() as conn:
                for run_name in conn.execute(stmt_select).first():
                    record_list.append({'id': run_id, 'name': run_name})

        run_df = pandas.DataFrame.from_records(record_list, index='id')

        return run_df

    def get_variant_to_chimera_borderline_df(self, filter_chimera_borderline_model):
        """Based on the SortedReadFile information TSV and run_model, returns the run_df

        :return: DataFrame with columns: index, name
        """

        filter_chimera_borderline_df = self.get_full_table_df(model=filter_chimera_borderline_model)

        variant_to_chimera_borderline_df = filter_chimera_borderline_df[['run_id', 'marker_id', 'variant_id', 'filter_delete']]
        variant_to_chimera_borderline_df = variant_to_chimera_borderline_df.sort_values(by='filter_delete', ascending=False)
        variant_to_chimera_borderline_df = variant_to_chimera_borderline_df\
            .drop_duplicates(subset=['run_id', 'marker_id', 'variant_id'], keep='first', inplace=False)
        variant_to_chimera_borderline_df.rename({'filter_delete': 'chimera_borderline'}, axis=1, inplace=True)

        return variant_to_chimera_borderline_df


class FastaInformationTSV(SampleInformationUtils):
    """Reads fasta information TSV file to a sample_information_id object"""

    def __init__(self, fasta_info_tsv, engine, include_tag_primer_fasta=False):
        """A new instance needs the path to the fasta information path information TSV file as wel as DB information to interact with the DB"""
        #
        self.fasta_information_df = pandas.read_csv(fasta_info_tsv, sep="\t", header=0)
        self.__engine = engine
        sample_information_id_df = self.__get_sample_information_df(include_tag_primer_fasta=include_tag_primer_fasta)
        super().__init__(engine, sample_information_id_df)


    def __get_sample_information_df(self, include_tag_primer_fasta=False):
        """Based on the SortedReadFile information TSV, returns a list of dictionnaries with run_id, marker_id, biosample_id
        and replicate entries (See return)

        :param tag_primer_fasta_information: Boolean. Default=False. If True, will also return tag and primer sequences and fasta file name
        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate': 1}, {'run_id': 1, ...
        """
        fasta_info_instance_list = []

        for row in self.fasta_information_df.itertuples():
            marker_name = row.Marker
            run_name = row.Run
            biosample_name = row.Biosample
            # replicate_name = row.replicate_name
            replicate = row.Replicate
            with self.__engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select([run_model.__table__.c.id]).where(run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = sqlalchemy.select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
            fasta_information_obj = {'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id, 'replicate': replicate}
            if include_tag_primer_fasta:
                fasta_information_obj['tag_fwd_sequence'] = row.tag_fwd_sequence
                fasta_information_obj['primer_fwd_sequence'] = row.primer_fwd_sequence
                fasta_information_obj['tag_rev_sequence'] = row.tag_rev_sequence
                fasta_information_obj['primer_rev_sequence'] = row.primer_rev_sequence
                fasta_information_obj['fasta_file_name'] = row.fasta_file_name
            # add this sample_instance ###########
            fasta_info_instance_list.append(fasta_information_obj)

        sample_information_id_df = pandas.DataFrame.from_records(data=fasta_info_instance_list).drop_duplicates(inplace=False)

        return sample_information_id_df


    def get_ids_df(self, include_tag_primer_fasta=False):

        """Based on the SortedReadFile information TSV, returns a list of dictionnaries with run_id, marker_id, biosample_id
        and replicate entries (See return)

        :param tag_primer_fasta_information: Boolean. Default=False. If True, will also return tag and primer sequences and fasta file name
        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate': 1}, {'run_id': 1, ...
        """

        fasta_info_ids_df = pandas.DataFrame()
        fasta_info_ids_df['run_id'] = None
        fasta_info_ids_df['marker_id'] = None
        fasta_info_ids_df['biosample_id'] = None

        # fasta_info_instance_list = []

        with self.__engine.connect() as conn:
            for k, row in self.fasta_information_df.iterrows():
                marker_name = row.Marker
                run_name = row.Run
                biosample_name = row.Biosample
                replicate = int(row.Replicate)
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select([run_model.__table__.c.id]).where(run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = sqlalchemy.select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                new_row_dic = dict(row)
                new_row_dic['run_id'] = run_id
                new_row_dic['marker_id'] = marker_id
                new_row_dic['biosample_id'] = biosample_id
                new_row_dic['replicate'] = replicate
                del new_row_dic["Run"]
                del new_row_dic["Marker"]
                del new_row_dic["Biosample"]
                del new_row_dic["Replicate"]

                fasta_info_ids_df = fasta_info_ids_df.append(pandas.DataFrame(new_row_dic, index=[0]), sort=False)
                # fasta_information_obj = {'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id, 'replicate': replicate}
                # add this sample_instance ###########
                # fasta_info_instance_list.append(fasta_information_obj)

        # sample_information_id_df = pandas.DataFrame.from_records(data=fasta_info_instance_list).drop_duplicates(inplace=False)
        return fasta_info_ids_df
