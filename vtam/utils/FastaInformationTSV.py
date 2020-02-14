import pandas

from vtam.utils.SampleInformationUtils import SampleInformationUtils


class FastaInformationTSV(SampleInformationUtils):
    """Reads fasta information TSV file to a sample_information_id object"""

    def __init__(self, fasta_info_tsv, engine, include_tag_primer_fasta=False):
        """A new instance needs the path to the fasta information path information TSV file as wel as DB information to interact with the DB"""
        #
        self.fasta_information_df = pandas.read_csv(fasta_info_tsv, sep="\t", header=0,
                                                    names=['tag_fwd_sequence', 'primer_fwd_sequence',
                                                           'tag_rev_sequence', 'primer_rev_sequence', 'marker_name',
                                                           'biosample_name', 'replicate', 'run_name', 'fastq_fwd',
                                                           'fastq_rev', 'fasta_file_name'])
        self.__engine = engine
        sample_information_id_df = self.__get_sample_information_df(include_tag_primer_fasta=include_tag_primer_fasta)
        super().__init__(engine, sample_information_id_df)


    def __get_sample_information_df(self, run_model, marker_model, biosample_model, include_tag_primer_fasta=False):
        """Based on the Fasta information TSV, returns a list of dictionnaries with run_id, marker_id, biosample_id
        and replicate entries (See return)

        :param tag_primer_fasta_information: Boolean. Default=False. If True, will also return tag and primer sequences and fasta file name
        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate': 1}, {'run_id': 1, ...
        """
        fasta_info_instance_list = []

        for row in self.fasta_information_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            # replicate_name = row.replicate_name
            replicate = row.replicate
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


