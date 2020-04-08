from sqlalchemy import bindparam

from vtam.utils.SampleInformationUtils import FastaInformationTSV
from vtam.utils.VariantReadCountLikeTable import VariantReadCountLikeTable
from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception
from wopmars.models.ToolWrapper import ToolWrapper

import pandas
import sys


class FilterVarious(object):
    pass


class FilterMinReplicateNumber(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterMinReplicateNumber"
    }

    # Input file
    __input_file_readinfo = "readinfo"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_variant_filter_lfn = "FilterLFN"
    # Output table
    __output_table_filter_min_replicate_number = "FilterMinReplicateNumber"


    def specify_input_file(self):
        return[
            FilterMinReplicateNumber.__input_file_readinfo,

        ]

    def specify_input_table(self):
        return [
            FilterMinReplicateNumber.__input_table_run,
            FilterMinReplicateNumber.__input_table_marker,
            FilterMinReplicateNumber.__input_table_biosample,
            FilterMinReplicateNumber.__input_table_variant_filter_lfn,
        ]


    def specify_output_table(self):
        return [
            FilterMinReplicateNumber.__output_table_filter_min_replicate_number,
        ]

    def specify_params(self):
        return {
            "min_replicate_number": "int",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        ################################################################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ################################################################################################################
        #
        # Input files
        fasta_info_tsv = self.input_file(FilterMinReplicateNumber.__input_file_readinfo)
        #
        # Input tables
        run_model = self.input_table(FilterMinReplicateNumber.__input_table_run)
        marker_model = self.input_table(FilterMinReplicateNumber.__input_table_marker)
        biosample_model = self.input_table(FilterMinReplicateNumber.__input_table_biosample)
        input_filter_lfn_model = self.input_table(FilterMinReplicateNumber.__input_table_variant_filter_lfn)
        #
        # Options
        min_replicate_number = self.option("min_replicate_number")
        # input_filter_lfn = self.option("input_filter_lfn")
        #
        # Output tables
        output_filter_min_replicate_model = self.output_table(FilterMinReplicateNumber.__output_table_filter_min_replicate_number)

        ################################################################################################################
        #
        # 1. Read readinfo to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ################################################################################################################

        fasta_info_tsv = FastaInformationTSV(engine=engine, fasta_info_tsv=fasta_info_tsv)

        ################################################################################################################
        #
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        variant_read_count_like_utils = VariantReadCountLikeTable(variant_read_count_like_model
                                                                  =output_filter_min_replicate_model, engine=engine)
        variant_read_count_like_utils.delete_from_db(sample_record_list=fasta_info_tsv.sample_record_list)

        ################################################################################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        # If filter_id exists, then previous Filter is FilterLFN and we select filter_id=8
        if 'filter_id' in [c.key for c in input_filter_lfn_model.__table__.columns]:
            filter_id = 8
        # Previous filter is not FilterLFN
        else:
            filter_id = None
        variant_read_count_df = fasta_info_tsv.get_variant_read_count_df(
            variant_read_count_like_model=input_filter_lfn_model, filter_id=filter_id)

        ################################################################################################################
        #
        # 4. Run Filter
        #
        ################################################################################################################

        filter_output_df = f9_delete_min_replicate_number(variant_read_count_df, min_replicate_number)

        ################################################################################################################
        #
        # Write to DB
        #
        ################################################################################################################

        record_list = VariantReadCountLikeTable.filter_delete_df_to_dict(filter_output_df)

        with engine.connect() as conn:

            # Delete instances that will be inserted
            del_stmt = output_filter_min_replicate_model.__table__.delete() \
                .where(output_filter_min_replicate_model.run_id == bindparam('run_id')) \
                .where(output_filter_min_replicate_model.marker_id == bindparam('marker_id')) \
                .where(output_filter_min_replicate_model.biosample_id == bindparam('biosample_id')) \
                .where(output_filter_min_replicate_model.replicate == bindparam('replicate'))
            conn.execute(del_stmt, record_list)

            # Insert new instances
            conn.execute(output_filter_min_replicate_model.__table__.insert(), record_list)

        ################################################################################################################
        #
        # Touch output tables, to update modification date
        #
        ################################################################################################################

        for output_table_i in self.specify_output_table():
            declarative_meta_i = self.output_table(output_table_i)
            obj = session.query(declarative_meta_i).order_by(declarative_meta_i.id.desc()).first()
            session.query(declarative_meta_i).filter_by(id=obj.id).update({'id': obj.id})
            session.commit()

        ################################################################################################################
        #
        # Exit vtam if all variants deleted
        #
        ################################################################################################################

        if filter_output_df.filter_delete.sum() == filter_output_df.shape[0]:
            Logger.instance().warning(VTAMexception("This filter has deleted all the variants: {}. "
                                                    "The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)


def f9_delete_min_replicate_number(variant_read_count_df, min_replicate_number=2):
    """
    This filter deletes variants if present in less than min_replicate_number replicates

    This filters deletes the variant if the count of the combinaison variant i and biosample j
    is low then the min_replicate_number.
    The deletion condition is: count(comb (N_ij) < min_replicate_number.

    Pseudo-algorithm of this function:

    1. Compute count(comb (N_ij)
    2. Set variant/biosample/replicate for deletion if count  column is low the min_replicate_number

    Updated:
    Jan 5, 2020

    :param variant_read_count_df: Variant read count dataframe
    :type variant_read_count_df: pandas.DataFrame

    :param min_replicate_number: Minimal number of replicates
    :type variant_read_count_input_df: int

    :return: The output of this filter is added to the 'self.variant_read_count_filter_delete_df' with 'filter_delete'=1 or 0
    :rtype: None
    """
    #
    df_filter_output=variant_read_count_df.copy()
    # replicate count
    df_grouped = variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).count().reset_index()
    df_grouped = df_grouped[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate']] # keep columns
    df_grouped = df_grouped.rename(columns={'replicate': 'replicate_count'})
    #
    df_filter_output['filter_delete'] = False
    df_filter_output = pandas.merge(df_filter_output, df_grouped, on=['run_id', 'marker_id', 'variant_id', 'biosample_id'], how='inner')
    df_filter_output.loc[df_filter_output.replicate_count < min_replicate_number, 'filter_delete'] = True
    #
    df_filter_output = df_filter_output[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate', 'read_count', 'filter_delete']]
    return df_filter_output

