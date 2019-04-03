import inspect

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
import os
from math import floor
from sqlalchemy import select
import pandas
from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.VSearch import VSearch1


class FilterPCRError(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FilterPCRError"
    }

    # Input file
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_variant = "Variant"
    __input_table_filter_min_replicate_number = "FilterMinReplicateNumber"
    __input_file_sample2fasta = "sample2fasta"
    # Output table
    __output_table_filter_pcr_error = "FilterPCRError"


    def specify_input_file(self):
        return[
            FilterPCRError.__input_file_sample2fasta,

        ]

    def specify_input_table(self):
        return [
            FilterPCRError.__input_table_marker,
            FilterPCRError.__input_table_run,
            FilterPCRError.__input_table_biosample,
            FilterPCRError.__input_table_replicate,
            FilterPCRError.__input_table_variant,
            FilterPCRError.__input_table_filter_min_replicate_number,
        ]


    def specify_output_table(self):
        return [
            FilterPCRError.__output_table_filter_pcr_error,
        ]

    def specify_params(self):
        return {
            "pcr_error_var_prop": "float",

        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        #
        # Input file path
        input_file_sample2fasta = self.input_file(FilterPCRError.__input_file_sample2fasta)
        #
        # Input table models
        marker_model = self.input_table(FilterPCRError.__input_table_marker)
        run_model = self.input_table(FilterPCRError.__input_table_run)
        biosample_model = self.input_table(FilterPCRError.__input_table_biosample)
        replicate_model = self.input_table(FilterPCRError.__input_table_replicate)
        variant_model = self.input_table(FilterPCRError.__input_table_variant)
        filter_min_replicate_number_model = self.input_table(FilterPCRError.__input_table_filter_min_replicate_number)
        #
        # Options
        pcr_error_var_prop = self.option("pcr_error_var_prop")
        #
        # Output table models
        filter_pcr_error_model = self.output_table(FilterPCRError.__output_table_filter_pcr_error)

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
            conn.execute(filter_pcr_error_model.__table__.delete(), sample_instance_list)

        ##########################################################
        #
        # 3. Select variant_read_count_model
        #
        ##########################################################

        filter_min_replicate_number_table = filter_min_replicate_number_model.__table__
        stmt_variant_filter_lfn = select([filter_min_replicate_number_table.c.marker_id,
                                          filter_min_replicate_number_table.c.run_id,
                                          filter_min_replicate_number_table.c.variant_id,
                                          filter_min_replicate_number_table.c.biosample_id,
                                          filter_min_replicate_number_table.c.replicate_id,
                                          filter_min_replicate_number_table.c.read_count])\
            .where(filter_min_replicate_number_table.c.filter_id == 9)\
            .where(filter_min_replicate_number_table.c.filter_delete == 0)
        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_read_count_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                    columns=['marker_id','run_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count'])

        ##########################################################
        #
        # 4. Select variant sequence
        #
        ##########################################################
        variant_model_table = variant_model.__table__
        stmt_variant = select([variant_model_table.c.id,
                               variant_model_table.c.sequence])
        # Select to DataFrame
        variant_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_list.append(row)
        variant_df = pandas.DataFrame.from_records(variant_list,
                    columns=['id', 'sequence'])

        ##########################################################
        #
        # 5. Run Vsearch
        #
        ##########################################################
        vsearch_output_df = f10_pcr_error_run_vsearch(variant_df)

        ##########################################################
        #
        # 6. Analyze vsearch output
        #
        ##########################################################
        df_filter_output = f10_pcr_error_analyze_vsearch_output_df(variant_read_count_df, vsearch_output_df, pcr_error_var_prop)

        ##########################################################
        #
        # 7. Insert Filter data
        #
        ##########################################################
        df_filter_output
        records = df_filter_output.to_dict('records')
        with engine.connect() as conn:
            conn.execute(filter_pcr_error_model.__table__.insert(), records)
def f10_pcr_error_run_vsearch(variant_df):
    # length of smallest sequence
    length_min = min(variant_df.sequence.apply(len).tolist())
    # calcul identity
    identity = floor((length_min - 1) / length_min * 100) / 100

    #
    ###################################################################
    # 5-1. Make a fasta file with all variants of the sample or replicate
    ###################################################################

    PathFinder.mkdir_p(os.path.join(tempdir, "f10_pcr_error"))
    variant_fasta = os.path.join(tempdir, "f10_pcr_error", '{}.fasta'.format("pcr_error"))
    with open(variant_fasta, 'w') as fout:
        for row in variant_df.itertuples():
            id = row.id
            sequence = row.sequence
            fout.write(">{}\n{}\n".format(id, sequence))
        ###################################################################
        # 5-2 Detect all pairs of variants with only 1 difference in the sequences and strong difference in abundance (readcounts)
        # 5-2.1 vsearch
        ###################################################################
    # import pdb; pdb.set_trace()
    sample_tsv = os.path.join(tempdir, '{}.tsv'.format("pcr_error"))
    vsearch_usearch_global_args = {'db': variant_fasta,
                                   'usearch_global': variant_fasta,
                                   'id': str(identity),
                                   'maxrejects': 0,
                                   'maxaccepts': 0,
                                   'userout': sample_tsv,
                                   'userfields': "query+target+alnlen+ids+mism+gaps",
                                   }

    vsearch_usearch_global = VSearch1(**vsearch_usearch_global_args)
    vsearch_usearch_global.run()
    column_names = ['query', 'target', 'alnlen', 'ids', 'mism', 'gaps']
    vsearch_output_df = pandas.read_csv(sample_tsv, sep='\t', names=column_names)
    return vsearch_output_df

def f10_pcr_error_analyze_vsearch_output_df(variant_read_count_df, vsearch_output_df, pcr_error_var_prop):
    # Output of filter pcr_error
    filter_output_df = variant_read_count_df.copy()
    filter_output_df['filter_id'] = 10
    filter_output_df['filter_delete'] = False
    #
    # Aggregate by biosample
    variant_read_count_grouped_df = variant_read_count_df.groupby(
        by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).sum().reset_index()
    # import pdb; pdb.set_trace()
    variant_read_count_grouped_df.drop(columns=['replicate_id'])
    # Compute sum mismatch and gap
    vsearch_output_df[
        'sum_mism_gaps'] = vsearch_output_df.mism + vsearch_output_df.gaps
    # Keep when sum mismatch and gap equals 1
    check_read_count_df = vsearch_output_df.loc[
        vsearch_output_df.sum_mism_gaps == 1, ['query', 'target']]
    # Add two colum the first for the variant id sequence query and the second for the target sequance variant id
    check_read_count_df = variant_read_count_grouped_df.merge(check_read_count_df, left_on=['variant_id'],
                                                              right_on=['query'])
    check_read_count_df = check_read_count_df.merge(variant_read_count_grouped_df,
                                                    left_on=['run_id', 'marker_id', 'biosample_id', 'target'],
                                                    right_on=['run_id', 'marker_id', 'biosample_id', 'variant_id'])

    check_read_count_df.drop(columns=['query'])
    check_read_count_df.drop(columns=['target'])
    check_read_count_df = check_read_count_df.rename(columns={'variant_id_x': 'variant_id'})
    check_read_count_df = check_read_count_df.rename(columns={'variant_id_y': 'variant_id_target'})
    check_read_count_df = check_read_count_df.rename(columns={'read_count_x': 'read_count'})
    check_read_count_df = check_read_count_df.rename(columns={'read_count_y': 'read_count_target'})

    # Add two column for the two expected ratio cases ratio 1 and ratio 2
    check_read_count_df['read_count_ratio'] = check_read_count_df.read_count / check_read_count_df.read_count_target
    check_read_count_df = (check_read_count_df.loc[check_read_count_df.read_count_ratio < pcr_error_var_prop])
    for row in check_read_count_df.itertuples():
        filter_output_df.loc[(filter_output_df['run_id'] == row.run_id)
                             & (filter_output_df['marker_id'] == row.marker_id)
                             & (filter_output_df['biosample_id'] == row.biosample_id)
                             & (filter_output_df['variant_id'] == row.variant_id), 'filter_delete'] = True
    return filter_output_df
