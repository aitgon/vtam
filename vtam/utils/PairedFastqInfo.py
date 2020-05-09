import argparse
import os
import pandas


class PairedFastqInfo(object):
    """Class for TSV file about the paired FASTQ files"""

    @staticmethod
    def read_tsv_into_df(fastqinfo_tsv_path):
        """
        Updated: Mai 8, 2020

        Parameters
        ----------
        fastqinfo_tsv_path: str
            Path to the FASTQ information files in TSV format

        Returns
        -------
        pandas.DataFrame

        """

        fastqinfo_df = pandas.read_csv(fastqinfo_tsv_path, sep='\t', header=0)
        fastqinfo_df.columns = fastqinfo_df.columns.str.lower()

        return fastqinfo_df

    @classmethod
    def check_fastqinfo_tsv(cls, fastqinfo_tsv_path):
        """
        Checks that the this input file is in the correct format. It is also used by the argparser.
        Updated: Mai 8, 2020

        Parameters
        ----------
        fastqinfo_tsv_path: str
            Path to the known_occurrences files in TSV format

        Returns
        -------

        """

        """Check fastqinfo_tsv_path format

        :param fastqinfo_tsv_path: Valid non-empty file fastqinfo_tsv_path
        :return: void

        """

        if not os.path.isfile(fastqinfo_tsv_path):
            raise argparse.ArgumentTypeError("The file '{}' does not exist. Please fix it.".format(
                fastqinfo_tsv_path))
        elif not os.stat(fastqinfo_tsv_path).st_size > 0:
            raise argparse.ArgumentTypeError("The file '{}' is empty!".format(fastqinfo_tsv_path))
        header_lower = {'marker', 'run', 'biosample', 'replicate', 'tagfwd', 'primerfwd', 'tagrev', 'primerrev', 'fastqfwd', 'fastqrev'}
        known_occurrences_df = cls.read_tsv_into_df(fastqinfo_tsv_path)
        if set(known_occurrences_df.columns) >= header_lower:  # contains at least the 'header_lower' columns
            return fastqinfo_tsv_path  # return the fastqinfo_tsv_path
        else:
            raise argparse.ArgumentTypeError("The format of file '{}' is wrong. Please fix it.".format(
                fastqinfo_tsv_path))
