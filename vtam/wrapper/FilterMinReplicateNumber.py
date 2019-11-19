from vtam.utils.FastaInformation import FastaInformation
from vtam.utils.VariantReadCountDF import VariantReadCountDF
from vtam.utils.VariantReadCountLikeTable import VariantReadCountLikeTable
from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

import pandas
import sys


class FilterVarious(object):
    pass


class FilterMinReplicateNumber(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterMinReplicateNumber"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_variant_filter_lfn = "FilterLFN"
    # Output table
    __output_table_filter_min_replicate_number = "FilterMinReplicateNumber"


    def specify_input_file(self):
        return[
            FilterMinReplicateNumber.__input_file_fastainfo,

        ]

    def specify_input_table(self):
        return [
            FilterMinReplicateNumber.__input_table_run,
            FilterMinReplicateNumber.__input_table_marker,
            FilterMinReplicateNumber.__input_table_biosample,
            FilterMinReplicateNumber.__input_table_replicate,
            FilterMinReplicateNumber.__input_table_variant_filter_lfn,
        ]


    def specify_output_table(self):
        return [
            FilterMinReplicateNumber.__output_table_filter_min_replicate_number,
        ]

    def specify_params(self):
        return {
            "min_replicate_number": "int",
            # "log_verbosity": "int",
            # "log_file": "str",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input files
        input_file_fastainfo = self.input_file(FilterMinReplicateNumber.__input_file_fastainfo)
        #
        # Input tables
        run_model = self.input_table(FilterMinReplicateNumber.__input_table_run)
        marker_model = self.input_table(FilterMinReplicateNumber.__input_table_marker)
        biosample_model = self.input_table(FilterMinReplicateNumber.__input_table_biosample)
        replicate_model = self.input_table(FilterMinReplicateNumber.__input_table_replicate)
        input_filter_lfn_model = self.input_table(FilterMinReplicateNumber.__input_table_variant_filter_lfn)
        #
        # Options
        min_replicate_number = self.option("min_replicate_number")
        #
        # Output tables
        output_filter_min_replicate_model = self.output_table(FilterMinReplicateNumber.__output_table_filter_min_replicate_number)

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################

        fasta_info = FastaInformation(input_file_fastainfo, engine, run_model, marker_model, biosample_model, replicate_model)
        fasta_info_record_list = fasta_info.get_fasta_information_record_list()

        ##########################################################
        #
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        variant_read_count_like_utils = VariantReadCountLikeTable(variant_read_count_like_model=output_filter_min_replicate_model, engine=engine)
        variant_read_count_like_utils.delete_output_filter_model(fasta_info_record_list=fasta_info_record_list)

        ##########################################################
        #
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        filter_id = 8 # Is filter_id from all FilterLFN all
        variant_read_count_df = fasta_info.get_variant_read_count_df(variant_read_count_like_model=input_filter_lfn_model, filter_id=filter_id)

        ##########################################################
        #
        # 4. Run Filter
        #
        ##########################################################
        filter_output_df = f9_delete_min_replicate_number(variant_read_count_df, min_replicate_number)

        ##########################################################
        #
        # Write to DB
        #
        ##########################################################

        record_list = VariantReadCountLikeTable.filter_delete_df_to_dict(filter_output_df)
        with engine.connect() as conn:
            conn.execute(output_filter_min_replicate_model.__table__.insert(), record_list)

        ##########################################################
        #
        # Exit vtam if all variants delete
        #
        ##########################################################

        try:
            assert not filter_output_df.filter_delete.sum() == filter_output_df.shape[0]
        except AssertionError:
            Logger.instance().warning(VTAMexception("This filter has deleted all the variants: {}. The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)


def f9_delete_min_replicate_number(variant_read_count_df, min_replicate_number=2):
    """
    This filter deletes variants if present in less than min_replicate_number replicates

    :param min_replicate_number: minimal number of replicates in which the variant must be present
    :return: None
    Non Low frequency noise filter minimum remplicant  (min_replicate_number) with a single threshold or several.
    Function IDs: 9 (min_replicate_number is 2)

    This filters deletes the variant if the count of the combinaison variant i and biosample j
    is low then the min_replicate_number.
    The deletion condition is: count(comb (N_ij) < min_replicate_number.


    Pseudo-algorithm of this function:

    1. Compute count(comb (N_ij)
    2. Set variant/biosample/replicate for deletion if count  column is low the min_replicate_number


    Updated:
    February 28, 2019

    Args:
       min_repln(float): Default deletion threshold


    Returns:
        None: The output of this filter is added to the 'self.variant_read_count_filter_delete_df'
        with filter_id='9' and 'filter_delete'=1 or 0


    """
    #
    df_filter_output=variant_read_count_df.copy()
    # replicate count
    df_grouped = variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).count().reset_index()
    df_grouped = df_grouped[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id']] # keep columns
    df_grouped = df_grouped.rename(columns={'replicate_id': 'replicate_count'})
    #
    df_filter_output['filter_delete'] = False
    df_filter_output = pandas.merge(df_filter_output, df_grouped, on=['run_id', 'marker_id', 'variant_id', 'biosample_id'], how='inner')
    df_filter_output.loc[df_filter_output.replicate_count < min_replicate_number, 'filter_delete'] = True
    #
    df_filter_output = df_filter_output[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count', 'filter_delete']]
    return df_filter_output

