import pandas
import sqlalchemy
import sys

from vtam import Logger, VTAMexception


class FastaInfo(object):

    def __init__(self, fasta_info_tsv, engine, run_model, marker_model, biosample_model, replicate_model):
        self.engine = engine
        self.run_model = run_model
        self.marker_model = marker_model
        self.biosample_model = biosample_model
        self.replicate_model = replicate_model
        #
        self.fastainfo_df = pandas.read_csv(fasta_info_tsv, sep="\t", header=0,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
                                                    'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        #
        self.fastainfo_instance_list = self.get_ids_of_run_marker_biosample_replicate()


    def get_ids_of_run_marker_biosample_replicate(self):
        """Returns a list of dictionnaries with run_id, marker_id, biosample_id and replicate_id entries (See return)

        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate_id': 1}, {'run_id': 1, ...
        """
        fastainfo_instance_list = []
        for row in self.fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with self.engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select([self.run_model.__table__.c.id]).where(self.run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([self.marker_model.__table__.c.id]).where(self.marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = sqlalchemy.select([self.biosample_model.__table__.c.id]).where(self.biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = sqlalchemy.select([self.replicate_model.__table__.c.id]).where(self.replicate_model.__table__.c.name==replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                fastainfo_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id,
                                             'replicate_id':replicate_id})
        return fastainfo_instance_list


    def get_variant_read_count_df(self, variant_read_count_model):
        """Get variant_read_count df from variant_read_count_model

        :param variant_read_count_model: SQLalchemy model with columns: run_id, marker_id, biosample_id, replicate_id, variant_id, read_count
        :return: DataFrame with columns: run_id, marker_id, biosample_id, replicate_id, variant_id, read_count
        """

        variant_read_count_table = variant_read_count_model.__table__

        variant_read_count_list = []
        for sample_instance in self.fastainfo_instance_list:
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate_id = sample_instance['replicate_id']
            stmt_select = sqlalchemy.select([variant_read_count_table.c.run_id,
                                  variant_read_count_table.c.marker_id,
                                  variant_read_count_table.c.biosample_id,
                                  variant_read_count_table.c.replicate_id,
                                  variant_read_count_table.c.variant_id,
                                  variant_read_count_table.c.read_count]).distinct()\
                                    .where(variant_read_count_table.c.run_id == run_id)\
                                    .where(variant_read_count_table.c.marker_id == marker_id)\
                                    .where(variant_read_count_table.c.biosample_id == biosample_id)\
                                    .where(variant_read_count_table.c.replicate_id == replicate_id)
            with self.engine.connect() as conn:
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
