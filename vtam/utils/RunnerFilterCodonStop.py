import Bio


class RunnerFilterCodonStop(object):

    def __init__(self, variant_read_count_df):
        """Carries out a chimera analysis"""
        self.variant_read_count_df = variant_read_count_df
        self.genetic_code = None

    def get_variant_read_count_delete_df(
            self,
            variant_df,
            genetic_code,
            skip_filter_codon_stop):

        variant_read_count_delete_df = self.variant_read_count_df.copy()
        variant_read_count_delete_df['filter_delete'] = False

        if not skip_filter_codon_stop:

            variant_has_stop_codon_df = self.annotate_stop_codon_count(
                variant_df, genetic_code)
            variants_with_stop_codons_list = variant_has_stop_codon_df.index[variant_has_stop_codon_df['has_stop_codon'] == 1].tolist(
            )

            variant_read_count_delete_df = self.variant_read_count_df.copy()
            variant_read_count_delete_df['filter_delete'] = False
            variant_read_count_delete_df.loc[variant_read_count_delete_df.variant_id.isin(
                variants_with_stop_codons_list), 'filter_delete'] = True

        return variant_read_count_delete_df

    def annotate_stop_codon_count(self, variant_df, genetic_code):
        """Takes a stats_df of variants and add the number of stop codons for each open reading frame

        Returns
        -------
        pandas DF
            VariantReadCountLikeModel are id, sequence, codon_stop_nb_frame1, codon_stop_nb_frame2, codon_stop_nb_frame3"""

        # For safety, convert to upper
        variant_df['sequence'] = variant_df['sequence'].str.upper()

        variant_has_stop_codon_df = variant_df.copy()
        variant_has_stop_codon_df['has_stop_codon'] = 0

        for row in variant_df.iterrows():
            id = row[0]
            sequence = row[1].sequence
            #
            variant_has_stop_codon_frame1 = self.seq_has_codon_stop(sequence=sequence, frame=1, genetic_code=genetic_code)
            variant_has_stop_codon_frame2 = self.seq_has_codon_stop(sequence=sequence, frame=2, genetic_code=genetic_code)
            variant_has_stop_codon_frame3 = self.seq_has_codon_stop(sequence=sequence, frame=3, genetic_code=genetic_code)

            # Check if all frames have stop codons
            if variant_has_stop_codon_frame1 and variant_has_stop_codon_frame2 and variant_has_stop_codon_frame3:
                variant_has_stop_codon_df.loc[variant_has_stop_codon_df.index ==
                                              id, 'has_stop_codon'] = 1
        return variant_has_stop_codon_df

    def seq_has_codon_stop(self, sequence, frame, genetic_code):
        """Takes one sequence and returns whether it has a stop codon or not

        Parameters
        ----------
        sequence: str
            DNA sequence in upper case
        frame : int
            Open reading frame index, 1,2,3

        Returns
        -------
        bool
            Has stop codon

        """
        if self.count_sequence_codon_stops(sequence, frame, genetic_code) > 0:
            return True
        else:
            return False

    def count_sequence_codon_stops(self, sequence, frame, genetic_code):
        """Takes one sequence and counts the number of stop codons for a given open reading frame and genetic genetic_code

        Parameters
        ----------
        sequence: str
            DNA sequence in upper case
        frame : int
            Open reading frame index, 1,2,3
        genetic_code : int
            NCBI genetic_codes: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

        Returns
        -------
        int
            Number of stop codons in the given open reading frame

        """
        # genetic table number 5: 'stop_codons': ['TAA', 'UAA', 'TAG', 'UAG']
        stop_codon_list = Bio.Data.CodonTable.generic_by_id[genetic_code].__dict__[
            'stop_codons']
        codon_list = [sequence[i:i + 3]
                      for i in range(frame - 1, len(sequence), 3)]
        codon_stop_count = sum([codon_list.count(stop_codon)
                                for stop_codon in stop_codon_list])
        return codon_stop_count
