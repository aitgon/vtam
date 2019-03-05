import sqlalchemy
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger


from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner
from sqlalchemy import select
import pandas


class FilterLFN(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FilterLFN"
    }

    # Input file
    __input_file_sample2fasta = "sample2fasta"
    # Input table
    __input_table_variant_read_count = "VariantReadCount"
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    # Output table
    __output_table_variant_filter_lfn = "FilterLFN"


    def specify_input_file(self):
        return[
            FilterLFN.__input_file_sample2fasta,
        ]

    def specify_input_table(self):
        return [
            FilterLFN.__input_table_variant_read_count,
            FilterLFN.__input_table_marker,
            FilterLFN.__input_table_run,
            FilterLFN.__input_table_biosample,
            FilterLFN.__input_table_replicate,
        ]


    def specify_output_table(self):
        return [
            FilterLFN.__output_table_variant_filter_lfn,
        ]

    def specify_params(self):
        return {
            "lfn_per_variant_threshold": "float",
            "lfn_per_replicate_threshold": "float",
            "lfn_per_biosample_per_replicate_threshold": "float",
            "lfn_read_count_threshold": "float",
        }

    # "min_replicate_number": "int",
    # "pcr_error_var_prop": "float",
    # "renkonen_number_of_replicate": "int",
    # "renkonen_renkonen_tail": "float"

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        #
        # Input file path
        __input_file_sample2fasta = self.input_file(FilterLFN.__input_file_sample2fasta)
        #
        # Input table models
        variant_read_count_model = self.input_table(FilterLFN.__input_table_variant_read_count)
        marker_model = self.input_table(FilterLFN.__input_table_marker)
        run_model = self.input_table(FilterLFN.__input_table_run)
        biosample_model = self.input_table(FilterLFN.__input_table_biosample)
        replicate_model = self.input_table(FilterLFN.__input_table_replicate)
        #
        # Output table models
        variant_filter_lfn_model = self.output_table(FilterLFN.__output_table_variant_filter_lfn)
        #
        # Create variant_biosample_replicate_df to run the filters with:
        # variant_id, biosample_id, replicate_id, read_count, variant_sequence
        # variant_model_table = variant_model.__table__
        variant_read_count_model_table = variant_read_count_model.__table__
        stmt_marker_id = select([variant_read_count_model_table.c.marker_id]).distinct()
        ##########################################################
        #
        # Delete marker/run/biosample/replicate from variant_read_count_model
        # Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################
        # Get ids for larter deletion
        sample2fasta_df = pandas.read_csv(__input_file_sample2fasta, sep="\t", header=None,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
            'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        variant_read_count_list = []
        for row in sample2fasta_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([replicate_model.__table__.c.id]).where(replicate_model.__table__.c.name==replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # Delete marker/run/biosample/replicate from variant_read_count_model
                conn.execute(variant_filter_lfn_model.__table__.delete()
                        .where(variant_filter_lfn_model.__table__.c.run_id == run_id)
                        .where(variant_filter_lfn_model.__table__.c.biosample_id == biosample_id)
                        .where(variant_filter_lfn_model.__table__.c.replicate_id == replicate_id)
                )
                # Select marker/run/biosample/replicate from variant_read_count_model
                stmt_select = select([col for col in variant_read_count_model.__table__.columns]).distinct()\
                        .where(variant_read_count_model.__table__.c.run_id == run_id)\
                        .where(variant_read_count_model.__table__.c.biosample_id == biosample_id)\
                        .where(variant_read_count_model.__table__.c.replicate_id == replicate_id)
                for row2 in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row2)
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['id', 'run_id', 'variant_id', 'marker_id', 'biosample_id', 'replicate_id', 'read_count'])
            #
        lfn_filter_runner = FilterLFNRunner(variant_read_count_df)
        #
        # FilterMinReplicateNumber parameters
        lfn_per_variant_threshold = self.option("lfn_per_variant_threshold")
        lfn_per_replicate_threshold = self.option("lfn_per_replicate_threshold")
        lfn_per_replicate_series_threshold = self.option("lfn_per_replicate_series_threshold")
        lfn_read_count_threshold = self.option("lfn_read_count_threshold")
        lfn_per_biosample_per_replicate_threshold = self.option("lfn_per_biosample_per_replicate_threshold")
        #
        Logger.instance().info("Launching LFN filter:")
        #
        ############################################
        # FilterMinReplicateNumber 2: f2_f4_lfn_delete_per_sum_variant
        ############################################
        lfn_filter_runner.f2_f4_lfn_delete_per_sum_variant(lfn_per_variant_threshold)
        #
        ############################################
        # FilterMinReplicateNumber  3: f3_f5_lfn_delete_per_sum_variant_replicate
        ############################################
        lfn_filter_runner.f3_f5_lfn_delete_per_sum_variant_replicate(lfn_per_replicate_threshold)

        ############################################
        # FilterMinReplicateNumber 6:  f6_lfn_delete_per_sum_biosample_replicate_delete
        ############################################

        lfn_filter_runner.f6_lfn_delete_per_sum_biosample_replicate(
            lfn_per_biosample_per_replicate_threshold)

        ############################################
        # FilterMinReplicateNumber  7:f7_lfn_delete_absolute_read_count
        ############################################
        lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)

        ############################################
        # FilterMinReplicateNumber 8:f8_lfn_delete_do_not_pass_all_filters
        ############################################
        lfn_filter_runner.f8_lfn_delete_do_not_pass_all_filters()

        ############################################
        # Write all LFN Filters
        ############################################
        records = lfn_filter_runner.delete_variant_df.to_dict('records')
        with engine.connect() as conn:
            conn.execute(variant_filter_lfn_model.__table__.insert(), records)
