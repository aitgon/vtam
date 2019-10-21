import pandas
import sqlalchemy


class FastaInfo(object):

    def __init__(self, fasta_info_tsv, engine):
        self.fastainfo_df = pandas.read_csv(fasta_info_tsv, sep="\t", header=0,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
                                                    'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        ################################################################################################################
        #
        # Compute IDs and return
        #
        ################################################################################################################



    def get_ids_of_run_marker_biosample_replicate(self, engine, run_model, marker_model, biosample_model, replicate_model):
        """Returns a list of dictionnaries with run_id, marker_id, biosample_id and replicate_id entries (See return)

        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate_id': 1}, {'run_id': 1, ...
        """
        fastainfo_instance_list = []
        for row in self.fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select([run_model.__table__.c.id]).where(run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = sqlalchemy.select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = sqlalchemy.select([replicate_model.__table__.c.id]).where(replicate_model.__table__.c.name==replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                fastainfo_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id,
                                             'replicate_id':replicate_id})
        return fastainfo_instance_list
