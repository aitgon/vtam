class DataframeVariant:
    """
    Functions to manipulate a variant_df
    """

    def __init__(self, variant_df):
        """Initiates a new object based on a variant_df (id, sequence)"""
        self.__variant_df = variant_df

    def to_fasta(self, fasta_path, add_column=None):
        """
        Takes variant DF with two columns (id, sequence) or optionally a third column (add_column)
        and returns a tsv_path to the fasta_path file

        Args:
            variant_df (pandas.DataFrame): DF with two columns (id, sequence)
            add_column (str): The add_column with be added to the fasta_path header as: >id;add_column_label=add_column_value
            fasta_path (str): Path to FASTA file

        Returns:
            None

        """
        with open(fasta_path, "w") as fout:
            for row in self.__variant_df.itertuples():
                fasta_line = ">{}".format(row.Index)
                if add_column is not None:
                    fasta_line += ";{}={}".format(add_column,
                                                  getattr(row, add_column))
                fasta_line += "\n"
                fasta_line += "{}\n".format(row.sequence)
                fout.write(fasta_line)
