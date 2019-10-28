class FastaIO:
    """
    Defines functions for the input and output to a fasta file
    """

    @staticmethod
    def variant_df_to_fasta_file(variant_df, fasta_path, add_column=None):
        """
        Takes variant DF with two columns (id, sequence) or optionally a third column (add_column)
        and returns a path to the fasta file

        Args:
            variant_df (pandas.DataFrame): DF with two columns (id, sequence)
            add_column (str): The add_column with be added to the fasta header as: >id;add_column_label=add_column_value
            fasta_path (str): Path to FASTA file

        Returns:
            None

        """
        with open(fasta_path, "w") as fout:
            for row in variant_df.itertuples():
                fasta_line = ">{}".format(row.id)
                if not add_column is None:
                    fasta_line += ";{}={}".format(add_column, getattr(row, add_column))
                fasta_line += "\n".format(row.id)
                fasta_line += "{}\n".format(row.sequence)
                fout.write(fasta_line)
