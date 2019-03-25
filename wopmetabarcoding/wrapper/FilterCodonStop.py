import inspect

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from sqlalchemy import select
import pandas





class FilterCodonStop(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FilterCodonStop"
    }

    # Input file
    __input_file_sample2fasta = "sample2fasta"
    __input_file_genetic_codes = "genetic_codes"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_filter_indel = "FilterIndel"
    __input_table_Variant = "Variant"
    # Output table
    __output_table_filter_codon_stop = "FilterCodonStop"



    def specify_input_file(self):
        return[
            FilterCodonStop.__input_file_sample2fasta,
            FilterCodonStop.__input_file_genetic_codes,
        ]

    def specify_input_table(self):
        return [
            FilterCodonStop.__input_table_marker,
            FilterCodonStop.__input_table_run,
            FilterCodonStop.__input_table_biosample,
            FilterCodonStop.__input_table_replicate,
            FilterCodonStop.__input_table_filter_indel,
            FilterCodonStop.__input_table_Variant,
        ]


    def specify_output_table(self):
        return [
            FilterCodonStop.__output_table_filter_codon_stop,

        ]



    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        #
        # Input file path
        input_file_sample2fasta = self.input_file(FilterCodonStop.__input_file_sample2fasta)
        input_file_genetic_codes = self.input_file(FilterCodonStop.__input_file_genetic_codes)
        #
        # Input table models
        marker_model = self.input_table(FilterCodonStop.__input_table_marker)
        run_model = self.input_table(FilterCodonStop.__input_table_run)
        indel_model = self.input_table(FilterCodonStop.__input_table_filter_indel)
        biosample_model = self.input_table(FilterCodonStop.__input_table_biosample)
        replicate_model = self.input_table(FilterCodonStop.__input_table_replicate)
        variant_model = self.input_table(FilterCodonStop.__input_table_Variant)
        #

        #
        # Output table models
        filter_codon_stop_model = self.output_table(FilterCodonStop.__output_table_filter_codon_stop)


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
        # 2. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################
        with engine.connect() as conn:
            conn.execute(filter_codon_stop_model.__table__.delete(), sample_instance_list)


        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        indel_model_table = indel_model.__table__
        stmt_variant_filter_lfn = select([indel_model_table.c.marker_id,
                                          indel_model_table.c.run_id,
                                          indel_model_table.c.variant_id,
                                          indel_model_table.c.biosample_id,
                                          indel_model_table.c.replicate_id,
                                          indel_model_table.c.filter_id,
                                          indel_model_table.c.filter_delete,
                                          indel_model_table.c.read_count])\
            .where(indel_model_table.c.filter_id == 13)\
            .where(indel_model_table.c.filter_delete == 0)
        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_read_count_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                    columns=['marker_id','run_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count'])
        # run_id, marker_id, variant_id, biosample_id, replicate_id, read_count, filter_id, filter_delete
        variant_model_table = variant_model.__table__
        stmt_variant = select([variant_model_table.c.id,
                               variant_model_table.c.sequence])

        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                                                              columns=['id', 'sequence'])
        ##########################################################
        #
        # 4. Run Filter
        #
        ##########################################################
        df_out = f14_filter_chimera(variant_read_count_df, variant_df)
        # df_filter_output = df_filter.copy()
        ##########################################################
        #
        # 5. Insert Filter data
        #
        ##########################################################
        records = df_out.to_dict('records')
        with engine.connect() as conn:
                conn.execute(filter_codon_stop_model.__table__.insert(), records)
        #


def f14_filter_chimera(variant_read_count_df, variant_df):
    """
    filter chimera
    """



