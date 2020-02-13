import pandas


class UserSampleSelection(object):
    """Carries out various operations in the DB based on the user sample selection"""

    def __init__(self, sample_selection_tsv):
        """The entry point is a TSV file given by the user with any combinations of columns among: Run, Marker,
        Biosample, Replicate"""

        self.sample_selection_df = pandas.read_csv(sample_selection_tsv, sep="\t", header=0)

        # self.__engine = engine
        # sample_information_id_df = self.__get_sample_information_df(run_model, marker_model, biosample_model, include_tag_primer_fasta=include_tag_primer_fasta)
        # super().__init__(engine, sample_information_id_df)

