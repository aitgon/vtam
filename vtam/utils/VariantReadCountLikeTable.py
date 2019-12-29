from sqlalchemy import select
from vtam import Logger, VTAMexception

import pandas
import sys


class VariantReadCountLikeTable(object):
    """Takes a any type of VariantReadCount models/table with at least run_id, marker_id, biosample_id, replicate, variant_id
    attributes/columns and performs various operations on it"""

    def __init__(self, engine, variant_read_count_like_model):
        # self.filter_name = filter_name
        self.engine = engine
        self.variant_read_count_like_model = variant_read_count_like_model

    def delete_records(self, record_list):
        """Deletes the entries in the output filter models based on a list of instances (dicts) defined by, run_id, marker_id, biosample_id, replicate"""

        with self.engine.connect() as conn:
            conn.execute(self.variant_read_count_like_model.__table__.delete(), record_list)

    def get_variant_read_count_df(self, fastainfo_instance_list, filter_id=None):
        """Get variant_read_count df from input filter models

        :param fasta_info_tsv: TSV file with the fasta_path information
        :return: DataFrame with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count
        """
        filter_min_replicate_number_table = self.input_variant_read_count_like_model.__table__

        variant_read_count_list = []
        for sample_instance in fastainfo_instance_list:
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate = sample_instance['replicate']
            if filter_id is None:
                stmt_select = select([filter_min_replicate_number_table.c.run_id,
                                      filter_min_replicate_number_table.c.marker_id,
                                      filter_min_replicate_number_table.c.biosample_id,
                                      filter_min_replicate_number_table.c.replicate,
                                      filter_min_replicate_number_table.c.variant_id,
                                      filter_min_replicate_number_table.c.read_count]).distinct()\
                                        .where(filter_min_replicate_number_table.c.run_id == run_id)\
                                        .where(filter_min_replicate_number_table.c.marker_id == marker_id)\
                                        .where(filter_min_replicate_number_table.c.biosample_id == biosample_id)\
                                        .where(filter_min_replicate_number_table.c.replicate == replicate)\
                                        .where(filter_min_replicate_number_table.c.filter_delete == 0)
            else:
                stmt_select = select([filter_min_replicate_number_table.c.run_id,
                                      filter_min_replicate_number_table.c.marker_id,
                                      filter_min_replicate_number_table.c.biosample_id,
                                      filter_min_replicate_number_table.c.replicate,
                                      filter_min_replicate_number_table.c.variant_id,
                                      filter_min_replicate_number_table.c.read_count]).distinct()\
                                        .where(filter_min_replicate_number_table.c.run_id == run_id)\
                                        .where(filter_min_replicate_number_table.c.marker_id == marker_id)\
                                        .where(filter_min_replicate_number_table.c.biosample_id == biosample_id)\
                                        .where(filter_min_replicate_number_table.c.replicate == replicate)\
                                        .where(filter_min_replicate_number_table.c.filter_id == filter_id)\
                                        .where(filter_min_replicate_number_table.c.filter_delete == 0)
            with self.engine.connect() as conn:
                for row2 in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row2)
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
            if 'replicate' in dir(row):
                instance['replicate'] = row.replicate
            if 'replicate_count' in dir(row):
                instance['replicate_count'] = row.replicate_count
            if 'read_count_average' in dir(row):
                instance['read_count_average'] = row.read_count_average
            records.append(instance)
        return records


