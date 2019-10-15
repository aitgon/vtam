import inspect
import os
import sys

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper


from sqlalchemy import select
import pandas

from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils.utilities import filter_delete_df_to_dict


class FilterMinReplicateNumber(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterMinReplicateNumber"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
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
            FilterMinReplicateNumber.__input_table_marker,
            FilterMinReplicateNumber.__input_table_run,
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
            "log_verbosity": "int",
            "log_file": "str",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file output
        input_file_fastainfo = self.input_file(FilterMinReplicateNumber.__input_file_fastainfo)
        #
        # Input table models
        marker_model = self.input_table(FilterMinReplicateNumber.__input_table_marker)
        run_model = self.input_table(FilterMinReplicateNumber.__input_table_run)
        biosample_model = self.input_table(FilterMinReplicateNumber.__input_table_biosample)
        replicate_model = self.input_table(FilterMinReplicateNumber.__input_table_replicate)
        variant_filter_lfn_model = self.input_table(FilterMinReplicateNumber.__input_table_variant_filter_lfn)
        #
        # Options
        min_replicate_number = self.option("min_replicate_number")
        #
        # Output table models
        filter_min_replicate_number_model = self.output_table(FilterMinReplicateNumber.__output_table_filter_min_replicate_number)

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        fastainfo_df = pandas.read_csv(input_file_fastainfo, sep="\t", header=0,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
            'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        sample_instance_list = []
        for row in fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id])\
                    .where(run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id])\
                    .where(marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id])\
                    .where(biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([replicate_model.__table__.c.id])\
                    .where(replicate_model.__table__.c.name==replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id,
                                             'replicate_id': replicate_id})

        ##########################################################
        #
        # 2. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################
        with engine.connect() as conn:
            conn.execute(filter_min_replicate_number_model.__table__.delete(), sample_instance_list)

        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################
        variant_filter_lfn_model_table = variant_filter_lfn_model.__table__
        variant_read_count_list = []
        for sample_instance in sample_instance_list:
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate_id = sample_instance['replicate_id']
            stmt_select = select([variant_filter_lfn_model_table.c.run_id,
                                  variant_filter_lfn_model_table.c.marker_id,
                                  variant_filter_lfn_model_table.c.biosample_id,
                                  variant_filter_lfn_model_table.c.replicate_id,
                                  variant_filter_lfn_model_table.c.variant_id,
                                  variant_filter_lfn_model_table.c.read_count]).distinct()\
                                    .where(variant_filter_lfn_model_table.c.run_id == run_id)\
                                    .where(variant_filter_lfn_model_table.c.marker_id == marker_id)\
                                    .where(variant_filter_lfn_model_table.c.biosample_id == biosample_id)\
                                    .where(variant_filter_lfn_model_table.c.replicate_id == replicate_id)\
                                    .where(variant_filter_lfn_model_table.c.filter_id == 8)\
                                    .where(variant_filter_lfn_model_table.c.filter_delete == 0)
            with engine.connect() as conn:
                for row2 in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row2)
        #
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count'])

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            sys.stderr.write("Error: No variants available for this filter: {}".format(os.path.basename(__file__)))
            sys.exit(1)

        ##########################################################
        #
        # 4. Run Filter
        #
        ##########################################################
        filter_output_df = f9_delete_min_replicate_number(variant_read_count_df, min_replicate_number)

        # ##########################################################
        # #
        # # 5. Insert Filter data
        # #
        # ##########################################################
        # records = df_filter_output.to_dict('records')
        # with engine.connect() as conn:
        #     conn.execute(filter_min_replicate_number_model.__table__.insert(), records)

        ############################################
        # Write to DB
        ############################################
        records = filter_delete_df_to_dict(filter_output_df)
        with engine.connect() as conn:
            conn.execute(filter_min_replicate_number_model.__table__.insert(), records)

        ##########################################################
        #
        # 6. Exit vtam if all variants delete
        #
        ##########################################################
        # Exit if no variants for analysis
        try:
            assert not filter_output_df.filter_delete.sum() == filter_output_df.shape[0]
        except AssertionError:
            Logger.instance().info(VTAMexception("Error: This filter has deleted all the variants"))
            sys.exit(1)


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

