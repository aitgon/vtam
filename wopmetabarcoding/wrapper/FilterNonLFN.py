import sqlalchemy
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger


from wopmetabarcoding.wrapper.FilterNonLFNutilities import FilterNonLFNRunner
from sqlalchemy import select
import pandas


class FilterNonLFN(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FilterLFN"
    }

    # Input file
    __input_file_sample2fasta = "sample2fasta"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_variant_read_count = "VariantReadCount"
    __input_table_variant = "Variant"
    # __input_table_min_replicate_number = "FilterMinReplicateNumber"

    # Output table
    __output_table_filter_non_lfn = "FilterNonLFN"


    def specify_input_file(self):
        return[
            FilterNonLFN.__input_file_sample2fasta,
        ]

    def specify_input_table(self):
        return [
            FilterNonLFN.__input_table_marker,
            FilterNonLFN.__input_table_run,
            FilterNonLFN.__input_table_biosample,
            FilterNonLFN.__input_table_replicate,
            FilterNonLFN.__input_table_variant_read_count,
            FilterNonLFN.__input_table_variant,
            # FilterNonLFN.__input_table_min_replicate_number,
        ]


    def specify_output_table(self):
        return [
            FilterNonLFN.__output_table_filter_non_lfn,
        ]

    def specify_params(self):
        return {
            "pcr_error_var_prop": "float",

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
        input_file_sample2fasta = self.input_file(FilterNonLFN.__input_file_sample2fasta)
        #
        # Input table models
        marker_model = self.input_table(FilterNonLFN.__input_table_marker)
        run_model = self.input_table(FilterNonLFN.__input_table_run)
        biosample_model = self.input_table(FilterNonLFN.__input_table_biosample)
        replicate_model = self.input_table(FilterNonLFN.__input_table_replicate)
        variant_read_count_model = self.input_table(FilterNonLFN.__input_table_variant_read_count)
        variant_model = self.input_table(FilterNonLFN.__input_table_variant)
        # filterminreplicatenumber_model =self.input_table(FilterNonLFN.__input_table_min_replicate_number)

        #
        # Output table models
        filter_non_lfn_model = self.output_table(FilterNonLFN.__output_table_filter_non_lfn)
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
                #get
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id, 'replicate_id':replicate_id})

        ##########################################################
        #
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################
        with engine.connect() as conn:
            conn.execute(filter_non_lfn_model.__table__.delete(), sample_instance_list)

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
        #select id/sequence from variant_model
        variant_list = []
        for sample_instance in sample_instance_list:
            id = sample_instance['id']
            sequence = sample_instance['sequence']
            stmt_select = select([col for col in variant_read_count_model.__table__.columns]).distinct() \
                .where(variant_model.__table__.c.id == id) \
                .where(variant_model.__table__.c.sequence == sequence) \

            with engine.connect() as conn:
                for row2 in conn.execute(stmt_select).fetchall():
                    variant_list.append(row2)
        #
        variant_df = pandas.DataFrame.from_records(variant_read_count_list,
                                                              columns=['id', 'sequence'])
        #

        non_lfn_filter_runner = FilterNonLFNRunner(variant_read_count_df,marker_id,variant_df)

        #
        # FilterLFN parameters
        pcr_error_var_prop = self.option("pcr_error_var_prop")




        #
        Logger.instance().info("Launching Non LFN filter:")
        #
        ############################################
        # FilterNonLFN 10: f10_pcr_error
        ############################################
        non_lfn_filter_runner.f10_pcr_error(pcr_error_var_prop)
        #

        ############################################
        # Write all LFN Filters
        ############################################
        records = non_lfn_filter_runner.delete_variant_df.to_dict('records') # modify delete_variant_df
        with engine.connect() as conn:
            conn.execute(filter_non_lfn_model.__table__.insert(), records)
