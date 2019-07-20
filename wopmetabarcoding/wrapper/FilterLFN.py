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
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_variant_read_count = "VariantReadCount"
    # Output table
    __output_table_filter_lfn = "FilterLFN"


    def specify_input_file(self):
        return[
            FilterLFN.__input_file_sample2fasta,
        ]

    def specify_input_table(self):
        return [
            FilterLFN.__input_table_marker,
            FilterLFN.__input_table_run,
            FilterLFN.__input_table_biosample,
            FilterLFN.__input_table_replicate,
            FilterLFN.__input_table_variant_read_count,
        ]


    def specify_output_table(self):
        return [
            FilterLFN.__output_table_filter_lfn,
        ]

    def specify_params(self):
        return {
            "filter_lfn_variant": "required|int",
            "lfn_variant_threshold": "float",
            "lfn_variant_replicate_threshold": "float",
            "lfn_biosample_replicate_threshold": "required|float",
            "lfn_read_count_threshold": "required|float",
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
        # Input file path
        input_file_sample2fasta = self.input_file(FilterLFN.__input_file_sample2fasta)
        #
        # Input table models
        run_model = self.input_table(FilterLFN.__input_table_run)
        marker_model = self.input_table(FilterLFN.__input_table_marker)
        biosample_model = self.input_table(FilterLFN.__input_table_biosample)
        replicate_model = self.input_table(FilterLFN.__input_table_replicate)
        variant_read_count_model = self.input_table(FilterLFN.__input_table_variant_read_count)
        #
        # Output table models
        variant_filter_lfn_model = self.output_table(FilterLFN.__output_table_filter_lfn)
        #
        # Options
        filter_lfn_variant = int(self.option("filter_lfn_variant"))
        lfn_variant_threshold = self.option("lfn_variant_threshold")
        lfn_variant_replicate_threshold = self.option("lfn_variant_replicate_threshold")
        lfn_biosample_replicate_threshold = self.option("lfn_biosample_replicate_threshold")
        lfn_read_count_threshold = self.option("lfn_read_count_threshold")
        #
        ##########################################################
        #
        # 1. Read sample2fasta to get run_id, marker_id, biosample_id, replicate_id for current analysis
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        # 4. Apply filters
        # 5. Write filters to variant_filter_lfn_model
        #
        ##########################################################

        ##########################################################
        #
        # 1. Read sample2fasta to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        sample2fasta_df = pandas.read_csv(input_file_sample2fasta, sep="\t", header=None,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
            'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        sample_instance_list = []
        for row in sample2fasta_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([replicate_model.__table__.c.id]).where(replicate_model.__table__.c.name==replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id, 'replicate_id':replicate_id})

        ##########################################################
        #
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################
        with engine.connect() as conn:
            conn.execute(variant_filter_lfn_model.__table__.delete(), sample_instance_list)

        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################
        variant_read_count_list = []
        for sample_instance in sample_instance_list:
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate_id = sample_instance['replicate_id']
            stmt_select = select([col for col in variant_read_count_model.__table__.columns]).distinct()\
                    .where(variant_read_count_model.__table__.c.run_id == run_id)\
                    .where(variant_read_count_model.__table__.c.marker_id == marker_id)\
                    .where(variant_read_count_model.__table__.c.biosample_id == biosample_id)\
                    .where(variant_read_count_model.__table__.c.replicate_id == replicate_id)
            with engine.connect() as conn:
                for row2 in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row2)
        #
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['id', 'run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count'])
        #
        lfn_filter_runner = FilterLFNRunner(variant_read_count_df)
        #
        Logger.instance().info("Launching LFN filter:")
        #
        ############################################
        # TaxAssign 2: f2_f4_lfn_delete_variant
        # Or
        # TaxAssign  3: f3_f5_lfn_delete_variant_replicate
        ############################################
        if bool(filter_lfn_variant):
            lfn_filter_runner.f2_f4_lfn_delete_variant(lfn_variant_threshold)
        else:
            lfn_filter_runner.f3_f5_lfn_delete_variant_replicate(lfn_variant_replicate_threshold)
        #

        ############################################
        # TaxAssign 6:  f6_lfn_delete_biosample_replicate_delete
        ############################################

        lfn_filter_runner.f6_lfn_delete_biosample_replicate(lfn_biosample_replicate_threshold)

        ############################################
        # TaxAssign  7:f7_lfn_delete_absolute_read_count
        ############################################
        lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)

        ############################################
        # TaxAssign 8:f8_lfn_delete_do_not_pass_all_filters
        ############################################
        lfn_filter_runner.f8_lfn_delete_do_not_pass_all_filters()

        ############################################
        # Delete current sample in TaxAssign
        ############################################

        ############################################
        # Write all LFN Filters
        ############################################
        records = lfn_filter_runner.delete_variant_df.to_dict('records')
        with engine.connect() as conn:
            conn.execute(variant_filter_lfn_model.__table__.insert(), records)
