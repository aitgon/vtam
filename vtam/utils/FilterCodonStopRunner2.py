import Bio


class FilterCodonStopRunner2(object):

    def __init__(self, code):
        """Carries out a chimera analysis"""
        self.code = code

    def annotate_stop_codon_count(self, variant_df):
        """Takes a df of variants and add the number of stop codons for each open reading frame

        Returns
        -------
        pandas DF
            Columns are id, sequence, codon_stop_nb_frame1, codon_stop_nb_frame2, codon_stop_nb_frame3"""


        variant_df['sequence'] = variant_df['sequence'].str.upper()  # For safety, convert to upper

        variant_has_stop_codon_df = variant_df.copy()
        variant_has_stop_codon_df['has_stop_codon'] = 0

        for row in variant_df.iterrows():
            id = row[1].id
            seq = row[1].sequence
            #
            if self.seq_has_codon_stop(seq=seq, frame=1):
                variant_has_stop_codon_df.loc[variant_has_stop_codon_df.id == id, 'has_stop_codon'] = 1
                continue
            elif self.seq_has_codon_stop(seq=seq, frame=2):
                variant_has_stop_codon_df.loc[variant_has_stop_codon_df.id == id, 'has_stop_codon'] = 1
                continue
            else:
                variant_has_stop_codon_df.loc[variant_has_stop_codon_df.id == id, 'has_stop_codon'] = 1
            #
            # stop_codon_count = self.count_codon_stop_nb_one_seq(seq=seq, frame=2)
            # variant_stop_codon_count_df.loc[variant_stop_codon_count_df.id == id, 'codon_stop_nb_frame2'] = stop_codon_count
            # #
            # stop_codon_count = self.count_codon_stop_nb_one_seq(seq=seq, frame=3)
            # variant_stop_codon_count_df.loc[variant_stop_codon_count_df.id == id, 'codon_stop_nb_frame3'] = stop_codon_count

        return variant_has_stop_codon_df


    def seq_has_codon_stop(self, seq, frame):
        """Takes one sequence and returns whether it has a stop codon or not

        Parameters
        ----------
        seq: str
            DNA seq in upper case
        frame : int
            Open reading frame index, 1,2,3

        Returns
        -------
        bool
            Has stop codon

        """
        if self.count_codon_stop_nb_one_seq(seq, frame) > 0:
            return True
        else:
            return False

    def count_codon_stop_nb_one_seq(self, seq, frame):
        """Takes one sequence and count the number of stop codons for a given open reading frame and genetic code

        Parameters
        ----------
        seq: str
            DNA seq in upper case
        frame : int
            Open reading frame index, 1,2,3

        Returns
        -------
        int
            Number of stop codons in the given open reading frame

        """
        # genetic table number 5: 'stop_codons': ['TAA', 'UAA', 'TAG', 'UAG']
        stop_codon_list = Bio.Data.CodonTable.generic_by_id[self.code].__dict__['stop_codons']
        codon_list = [seq[i:i + 3] for i in range(frame - 1, len(seq), 3)]
        codon_stop_count = sum([codon_list.count(stop_codon) for stop_codon in stop_codon_list])
        return codon_stop_count
