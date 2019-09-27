import inspect

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.utils.PathManager import PathFinder
from wopmetabarcoding.utils.VSearch import VSearch1, Vsearch2, Vsearch3
from Bio import SeqIO
import os

from sqlalchemy import select
import pandas

from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.utils.utilities import create_step_tmp_dir, tempdir


class FilterChimera(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FilterChimera"
    }

    # Input file
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_file_sample2fasta = "sample2fasta"
    __input_table_filter_pcr_error = "FilterPCRError"
    __input_table_Variant = "Variant"
    # Output table
    __output_table_filter_chimera = "FilterChimera"
    __output_table_filter_chimera_borderline = "FilterChimeraBorderline"


    def specify_input_file(self):
        return[
            FilterChimera.__input_file_sample2fasta,

        ]

    def specify_input_table(self):
        return [
            FilterChimera.__input_table_marker,
            FilterChimera.__input_table_run,
            FilterChimera.__input_table_biosample,
            FilterChimera.__input_table_replicate,
            FilterChimera.__input_table_filter_pcr_error,
            FilterChimera.__input_table_Variant,
        ]


    def specify_output_table(self):
        return [
            FilterChimera.__output_table_filter_chimera,
            FilterChimera.__output_table_filter_chimera_borderline,
        ]



    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        this_step_tmp_dir = create_step_tmp_dir(__file__)

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file path
        input_file_sample2fasta = self.input_file(FilterChimera.__input_file_sample2fasta)
        #
        # Input table models
        marker_model = self.input_table(FilterChimera.__input_table_marker)
        run_model = self.input_table(FilterChimera.__input_table_run)
        pcr_error_model = self.input_table(FilterChimera.__input_table_filter_pcr_error)
        biosample_model = self.input_table(FilterChimera.__input_table_biosample)
        replicate_model = self.input_table(FilterChimera.__input_table_replicate)
        variant_model = self.input_table(FilterChimera.__input_table_Variant)
        #
        # Output table models
        filter_chimera_model = self.output_table(FilterChimera.__output_table_filter_chimera)
        filter_chimera_borderline_model = self.output_table(FilterChimera.__output_table_filter_chimera_borderline)

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
            conn.execute(filter_chimera_model.__table__.delete(), sample_instance_list)
        #
        with engine.connect() as conn:
            conn.execute(filter_chimera_borderline_model.__table__.delete(), sample_instance_list)

        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        pcr_error_model_table = pcr_error_model.__table__
        stmt_variant_filter_lfn = select([pcr_error_model_table.c.marker_id,
                                          pcr_error_model_table.c.run_id,
                                          pcr_error_model_table.c.variant_id,
                                          pcr_error_model_table.c.biosample_id,
                                          pcr_error_model_table.c.replicate_id,
                                          pcr_error_model_table.c.read_count])\
            .where(pcr_error_model_table.c.filter_delete == 0)
        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_read_count_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                    columns=['marker_id','run_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count'])
        # run_id, marker_id, variant_id, biosample_id, replicate_id, read_count, filter_delete
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
        df_chimera, df_chimera_borderline = f11_filter_chimera(variant_read_count_df, variant_df, this_step_tmp_dir)
        # df_filter_output = df_filter.copy()
        ##########################################################
        #
        # 5. Insert Filter data
        #
        ##########################################################
        with engine.connect() as conn:
                conn.execute(filter_chimera_model.__table__.insert(), df_chimera.to_dict('records'))
        #
        with engine.connect() as conn:
                conn.execute(filter_chimera_borderline_model.__table__.insert(), df_chimera_borderline.to_dict('records'))


def f11_filter_chimera(variant_read_count_df, variant_df, this_step_tmp_dir):
    """
    filter chimera
    """

    logger.debug(
        "file: {}; line: {}".format(__file__, inspect.currentframe().f_lineno))
    #
    filter_output_df = variant_read_count_df.copy()
    filter_output_df['filter_delete'] = False
    #
    filter_borderline_output_df = variant_read_count_df.copy()
    filter_borderline_output_df['filter_delete'] = False
    #
    df_grouped_variant_read_count = variant_read_count_df.groupby(
        by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).sum().reset_index()
    df_grouped_variant_read_count = df_grouped_variant_read_count[
        ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'read_count']]  # keep columns
    #
    df_grouped_run_marker_biosample = variant_read_count_df.groupby(
        by=['run_id', 'marker_id', 'biosample_id']).sum().reset_index()
    df_grouped_run_marker_biosample = df_grouped_run_marker_biosample[
        ['run_id', 'marker_id', 'biosample_id']]  # keep columns
    ###################################################################
    #
    # 1. Make a fasta file with all variants for each run, marker, biosample
    #
    ###################################################################
    for row_run_marker_biosample in df_grouped_run_marker_biosample.iterrows():
        run_id = row_run_marker_biosample[1].run_id
        marker_id = row_run_marker_biosample[1].marker_id
        biosample_id = row_run_marker_biosample[1].biosample_id
        #
        df_grouped_biosample_id = df_grouped_variant_read_count.loc[
            df_grouped_variant_read_count.biosample_id == biosample_id, ['variant_id', 'read_count']]

        ###################################################################
        #
        # 2. Sort variants by abundance and write to fasta
        #
        ###################################################################
        df_grouped_biosample_id.sort_values(by='read_count', ascending=False, inplace=True)
        #
        chimera1_fasta = os.path.join(this_step_tmp_dir,
                                      'chimera1_biosample_id_{}.fasta'.format(biosample_id))
        #
        #  Prepare 1 fasta file
        with open(chimera1_fasta, 'w') as fasta_fout:
            for variant_id in df_grouped_biosample_id.variant_id.unique():
                variant_sequence = variant_df.loc[variant_df.id == variant_id, 'sequence'].values[0]
                read_count = df_grouped_biosample_id.loc[
                    df_grouped_variant_read_count.variant_id == variant_id, 'read_count'].values[0]
                fasta_fout.write('>{};size={}\n{}\n'.format(variant_id, read_count, variant_sequence))

        ###################################################################
        #
        # 3. Run uchime_denovo
        #
        ###################################################################
        chimera2_borderline_fasta = os.path.join(this_step_tmp_dir,
                                                 'chimera2_biosample_id_{}_borderline.fasta'.format(biosample_id))
        chimera2_nonchimeras_fasta = os.path.join(this_step_tmp_dir,
                                                  'chimera2_biosample_id_{}_nonchimeras.fasta'.format(biosample_id))
        chimera2_chimeras_fasta = os.path.join(this_step_tmp_dir,
                                               'chimera2_biosample_id_{}_chimeras.fasta'.format(biosample_id))
        #
        vsearch_chimera_args = {
            "uchime_denovo": chimera1_fasta,
            "borderline": chimera2_borderline_fasta,
            "nonchimeras": chimera2_nonchimeras_fasta,
            "chimeras": chimera2_chimeras_fasta
        }
        vsearch_chimera = Vsearch3(**vsearch_chimera_args)
        vsearch_chimera.run()

        ###################################################################
        #
        # 4. Delete variant from replicate/sample if chimeras
        #
        ###################################################################
        with open(chimera2_chimeras_fasta, "r") as handle:
            for chimera_seqrecord in SeqIO.parse(handle, "fasta"):
                variant_id = int(chimera_seqrecord.id.split(';')[0])
                filter_output_df.loc[(filter_output_df['run_id'] == run_id)
                                     & (filter_output_df['marker_id'] == marker_id)
                                     & (filter_output_df['biosample_id'] == biosample_id)
                                     & (filter_output_df['variant_id'] == variant_id), 'filter_delete'] = True

        with open(chimera2_borderline_fasta, "r") as handle:
            for chimera_seqrecord in SeqIO.parse(handle, "fasta"):
                variant_id = int(chimera_seqrecord.id.split(';')[0])
                filter_borderline_output_df.loc[(filter_borderline_output_df['run_id'] == run_id)
                                     & (filter_borderline_output_df['marker_id'] == marker_id)
                                     & (filter_borderline_output_df['biosample_id'] == biosample_id)
                                     & (filter_borderline_output_df['variant_id'] == variant_id), 'filter_delete'] = True
        #

    return filter_output_df, filter_borderline_output_df

