import argparse
import os
import sys

import pandas
import sqlalchemy

from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception

from vtam.models.Run import Run
from vtam.models.Marker import Marker
from vtam.models.Biosample import Biosample


class SampleInformationFile:
    """Sample information file (paired fastq, merged fasta and sorted read fasta)"""

    def __init__(self, tsv_path):
        self.tsv_path = tsv_path

    def read_tsv_into_df(self):
        """
        Updated: Mai 10, 2020

        Parameters
        ----------
        tsv_path: str
            Path to the information files in TSV format

        Returns
        -------
        pandas.DataFrame

        """

        fastqinfo_df = pandas.read_csv(self.tsv_path, sep='\t', header=0)
        fastqinfo_df.columns = fastqinfo_df.columns.str.lower()

        return fastqinfo_df

    def check_args(self, header):
        """
        Checks that the this input file is in the correct format. It is also used by the argparser.
        Updated: Mai 8, 2020

        Parameters
        ----------
        header: list
            List of header fields in lower case of SampleInformationFile

        Returns
        -------

        """

        """Check tsv_path format

        :param tsv_path: Valid non-empty file tsv_path
        :return: void

        """

        if not os.path.isfile(self.tsv_path):
            raise argparse.ArgumentTypeError("The file '{}' does not exist. Please fix it.".format(
                self.tsv_path))
        elif not os.stat(self.tsv_path).st_size > 0:
            raise argparse.ArgumentTypeError("The file '{}' is empty!".format(self.tsv_path))
        sample_info_df = self.read_tsv_into_df()
        if set(sample_info_df.columns) >= header:  # contains at least the 'header_lower' columns
            return self.tsv_path  # return the tsv_path
        else:
            raise argparse.ArgumentTypeError("The format of file '{}' is wrong. Please fix it.".format(
                self.tsv_path))

    def to_identifier_df(self, engine):
        """Takes the sample information df and replaces names with ids

        Returns
        -------
        pandas.DataFrame
        DF with ids instead of names and columns run_id, marker_id, biosample_id

        """

        sample_info_df = pandas.DataFrame()
        sample_info_df['run_id'] = None
        sample_info_df['marker_id'] = None
        sample_info_df['biosample_id'] = None

        with engine.connect() as conn:
            for k, row in self.read_tsv_into_df().iterrows():
                marker_name = row.marker
                run_name = row.run
                biosample_name = row.biosample
                replicate = int(row.replicate)
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select([Run.__table__.c.id]).where(Run.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([Marker.__table__.c.id]).where(Marker.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = sqlalchemy.select([Biosample.__table__.c.id]).where(Biosample.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                new_row_dic = dict(row)
                new_row_dic['run_id'] = run_id
                new_row_dic['marker_id'] = marker_id
                new_row_dic['biosample_id'] = biosample_id
                new_row_dic['replicate'] = replicate
                del new_row_dic["run"]
                del new_row_dic["marker"]
                del new_row_dic["biosample"]

                sample_info_df = sample_info_df.append(pandas.DataFrame(new_row_dic, index=[0]), sort=False)

        sample_info_df.columns = sample_info_df.columns.str.lower()
        return sample_info_df

    def get_variant_read_count_df(self, engine, variant_read_count_like_model, filter_id=None):
        """
        Based on the SortedReadFile samples and the variant_read_count_model, returns the variant_read_count_input_df

        Parameters
        ----------
        engine SQLalchemy engine
        variant_read_count_like_model: SQLaclchemy declarative model
        SQLalchemy models with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count

        filter_id

        Returns
        -------
        DataFrame with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count

        """

        variant_read_count_like_table = variant_read_count_like_model.__table__

        variant_read_count_list = []
        for sample_instance_row in self.to_identifier_df().itertuples():
            run_id = sample_instance_row.run_id
            marker_id = sample_instance_row.marker_id
            biosample_id = sample_instance_row.biosample_id
            replicate = sample_instance_row.replicate
            stmt_select = sqlalchemy.select([variant_read_count_like_table.c.run_id,
                                                variant_read_count_like_table.c.marker_id,
                                                variant_read_count_like_table.c.biosample_id,
                                                variant_read_count_like_table.c.replicate,
                                                variant_read_count_like_table.c.variant_id,
                                                variant_read_count_like_table.c.read_count]).distinct()\
                                                .where(variant_read_count_like_table.c.run_id == run_id)\
                                                .where(variant_read_count_like_table.c.marker_id == marker_id)\
                                                .where(variant_read_count_like_table.c.biosample_id == biosample_id)\
                                                .where(variant_read_count_like_table.c.replicate == replicate)
            # Used for filters tables where filter_delete attribute exists
            if 'filter_delete' in [column.key for column in variant_read_count_like_table.columns]:
                stmt_select = stmt_select.where(variant_read_count_like_table.c.filter_delete == 0)
            # used for filter lfn where filter_id = 8 is necessary (do not pass all filters)
            if not filter_id is None:
                stmt_select = stmt_select.where(variant_read_count_like_table.c.filter_id == filter_id)

            with engine.connect() as conn:
                for row in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row)
        #
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['run_id', 'marker_id', 'biosample_id', 'replicate', 'variant_id', 'read_count'])

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            Logger.instance().warning(VTAMexception("No variants available after the Codon Stop filter. "
                                                    "The pipeline will stop here.".format(self.__class__.__name__)))
            sys.exit(0)
        return variant_read_count_df

