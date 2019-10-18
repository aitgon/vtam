from sqlalchemy import select
from vtam import Logger, VTAMexception

import pandas
import sys


class FilterCommon(object):


    def __init__(self, engine, run_model, marker_model, biosample_model, replicate_model, input_filter_model, output_filter_models):
        self.engine = engine
        self.run_model = run_model
        self.marker_model = marker_model
        self.biosample_model = biosample_model
        self.replicate_model = replicate_model
        self.input_filter_model = input_filter_model
        self.output_filter_models = output_filter_models

    def get_fastainfo_instance_list_with_ids(self, input_file_fastainfo):
        """
        Returns a list of dictionnaries with the the ids of the TSV file with the fasta informations. Example
    [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate_id': 1}, {'run_id': 1, 'marker_id': 1, 'biosample_id': 2, 'replicate_id': 1}, {'run_id': 1, 'marker_id': 1, 'biosample_id': 3, 'replicate_id': 1}, {'run_id': 1, 'marker_id': 1, 'biosample_id': 4, 'replicate_id': 1}, {'run_id': 1, 'marker_id': 1, 'biosample_id': 5, 'replicate_id': 1}, {'run_id': 1, 'marker_id': 1, 'biosample_id': 6, 'replicate_id': 1}, {'run_id': 1, 'marker_id': 1, 'biosample_id': 7, 'replicate_id': 1}, {'run_id': 1, 'marker_id': 1, 'biosample_id': 8, 'replicate_id': 1}

        :param input_file_fastainfo: TSV file with the fasta information
        :return: list of dictionnaries
        """
        fastainfo_df = pandas.read_csv(input_file_fastainfo, sep="\t", header=0,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
            'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        fastainfo_instance_list = []
        for row in fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with self.engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([self.run_model.__table__.c.id]).where(self.run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([self.marker_model.__table__.c.id]).where(self.marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([self.biosample_model.__table__.c.id]).where(self.biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([self.replicate_model.__table__.c.id]).where(self.replicate_model.__table__.c.name==replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                fastainfo_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id,
                                             'replicate_id':replicate_id})
        return fastainfo_instance_list



    def delete_output_filter_model(self, fastainfo_instance_list):
        """Deletes the entries in the output filter model based on a list of instances (dicts) defined by, run_id, marker_id, biosample_id, replicate_id"""
        with self.engine.connect() as conn:
            if isinstance(self.output_filter_models, list):
                for output_filter_model in self.output_filter_models:
                    conn.execute(output_filter_model.__table__.delete(), fastainfo_instance_list)
            else:
                conn.execute(self.output_filter_models.__table__.delete(), fastainfo_instance_list)



    def get_variant_read_count_model(self, fastainfo_instance_list, filter_id=None):
        """Get variant_read_count df from input filter model

        :param input_file_fastainfo: TSV file with the fasta information
        :return: DataFrame with columns: run_id, marker_id, biosample_id, replicate_id, variant_id, read_count
        """
        filter_min_replicate_number_table = self.input_filter_model.__table__

        variant_read_count_list = []
        for sample_instance in fastainfo_instance_list:
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate_id = sample_instance['replicate_id']
            if filter_id is None:
                stmt_select = select([filter_min_replicate_number_table.c.run_id,
                                      filter_min_replicate_number_table.c.marker_id,
                                      filter_min_replicate_number_table.c.biosample_id,
                                      filter_min_replicate_number_table.c.replicate_id,
                                      filter_min_replicate_number_table.c.variant_id,
                                      filter_min_replicate_number_table.c.read_count]).distinct()\
                                        .where(filter_min_replicate_number_table.c.run_id == run_id)\
                                        .where(filter_min_replicate_number_table.c.marker_id == marker_id)\
                                        .where(filter_min_replicate_number_table.c.biosample_id == biosample_id)\
                                        .where(filter_min_replicate_number_table.c.replicate_id == replicate_id)\
                                        .where(filter_min_replicate_number_table.c.filter_delete == 0)
            else:
                stmt_select = select([filter_min_replicate_number_table.c.run_id,
                                      filter_min_replicate_number_table.c.marker_id,
                                      filter_min_replicate_number_table.c.biosample_id,
                                      filter_min_replicate_number_table.c.replicate_id,
                                      filter_min_replicate_number_table.c.variant_id,
                                      filter_min_replicate_number_table.c.read_count]).distinct()\
                                        .where(filter_min_replicate_number_table.c.run_id == run_id)\
                                        .where(filter_min_replicate_number_table.c.marker_id == marker_id)\
                                        .where(filter_min_replicate_number_table.c.biosample_id == biosample_id)\
                                        .where(filter_min_replicate_number_table.c.replicate_id == replicate_id)\
                                        .where(filter_min_replicate_number_table.c.filter_id == filter_id)\
                                        .where(filter_min_replicate_number_table.c.filter_delete == 0)
            with self.engine.connect() as conn:
                for row2 in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row2)
        #
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count'])

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            Logger.instance().warning(VTAMexception("No variants available for this Filter: {}. The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)
        return variant_read_count_df



    ##########################################################
    #
    # Convert DF to list of dictionaries to use in an sqlalchemy core insert
    #
    ##########################################################
    @staticmethod
    def filter_delete_df_to_dict(filter_df):
        """Convert DF to list of dictionaries to use in an sqlalchemy core insert"""
        records = []
        # import pdb; pdb.set_trace()
        for row in filter_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            variant_id = row.variant_id
            read_count = row.read_count
            instance = {'run_id': run_id, 'marker_id': marker_id,
                                                         'variant_id': variant_id, 'biosample_id': biosample_id,
                                                         'read_count': read_count}
            if 'filter_delete' in dir(row):
                instance['filter_delete'] = row.filter_delete
            if 'filter_id' in dir(row):
                instance['filter_id'] = row.filter_id
            if 'replicate_id' in dir(row):
                instance['replicate_id'] = row.replicate_id
            if 'replicate_count' in dir(row):
                instance['replicate_count'] = row.replicate_count
            if 'read_count_average' in dir(row):
                instance['read_count_average'] = row.read_count_average
            records.append(instance)
        return records


