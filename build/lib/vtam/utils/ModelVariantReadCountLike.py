import pandas
import sqlalchemy


class ModelVariantReadCountLike(object):
    """Takes a any type of VariantReadCount models/table with at least run_id, marker_id, sample_id, replicate, variant_id
    attributes/columns and performs various operations on it"""

    def __init__(self, engine, variant_read_count_like_model):
        # self.filter_name = filter_name
        self.engine = engine
        self.variant_read_count_like_model = variant_read_count_like_model

    def delete_from_db(self, sample_record_list):
        """Deletes the entries in the output filter models based on a list of instances (dicts) defined by, run_id,
         marker_id, sample_id, replicate"""

        sample_column_list = pandas.DataFrame(
            sample_record_list).columns.tolist()

        with self.engine.connect() as conn:
            stmt = self.variant_read_count_like_model.__table__.delete()
            stmt = stmt.where(
                self.variant_read_count_like_model.__table__.c.run_id == sqlalchemy.bindparam('run_id'))
            stmt = stmt.where(
                self.variant_read_count_like_model.__table__.c.marker_id == sqlalchemy.bindparam('marker_id'))
            if 'sample_id' in sample_column_list:
                stmt = stmt.where(
                    self.variant_read_count_like_model.__table__.c.sample_id == sqlalchemy.bindparam('sample_id'))
            if 'replicate' in sample_column_list and 'replicate' in [
                    col.key for col in self.variant_read_count_like_model.__table__.columns]:
                stmt = stmt.where(
                    self.variant_read_count_like_model.__table__.c.replicate == sqlalchemy.bindparam('replicate'))
            conn.execute(stmt, sample_record_list)

    ##########################################################
    #
    # Convert DF to list of dictionaries to use in an sqlalchemy core insert
    #
    ##########################################################

    @staticmethod
    def filter_delete_df_to_dict(filter_df):
        """Convert DF to list of dictionaries to use in an sqlalchemy core insert"""
        records = []
        for row in filter_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            sample_id = row.sample_id
            variant_id = row.variant_id
            read_count = row.read_count
            instance = {'run_id': run_id, 'marker_id': marker_id,
                                                         'variant_id': variant_id, 'sample_id': sample_id,
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
