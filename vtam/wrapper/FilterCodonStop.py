import os
import sys

import Bio
import inspect

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from sqlalchemy import select
import pandas

from vtam import OptionManager, VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils.utilities import filter_delete_df_to_dict


class FilterCodonStop(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterCodonStop"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
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
            FilterCodonStop.__input_file_fastainfo,
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

    def specify_params(self):
        return {
            "genetic_table_number": "int",
            "log_verbosity": "int",
            "log_file": "str"

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
        input_file_fastainfo = self.input_file(FilterCodonStop.__input_file_fastainfo)
        #
        # Input table models
        marker_model = self.input_table(FilterCodonStop.__input_table_marker)
        run_model = self.input_table(FilterCodonStop.__input_table_run)
        indel_model = self.input_table(FilterCodonStop.__input_table_filter_indel)
        biosample_model = self.input_table(FilterCodonStop.__input_table_biosample)
        replicate_model = self.input_table(FilterCodonStop.__input_table_replicate)
        variant_model = self.input_table(FilterCodonStop.__input_table_Variant)
        #options
        genetic_table_number = int(self.option("genetic_table_number"))
        #
        # Output table models
        filter_codon_stop_model = self.output_table(FilterCodonStop.__output_table_filter_codon_stop)


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

        variant_read_count_list = []
        for sample_instance in sample_instance_list:
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate_id = sample_instance['replicate_id']
            stmt_select = select([indel_model_table.c.run_id,
                                  indel_model_table.c.marker_id,
                                  indel_model_table.c.biosample_id,
                                  indel_model_table.c.replicate_id,
                                  indel_model_table.c.variant_id,
                                  indel_model_table.c.read_count]).distinct()\
                                    .where(indel_model_table.c.run_id == run_id)\
                                    .where(indel_model_table.c.marker_id == marker_id)\
                                    .where(indel_model_table.c.biosample_id == biosample_id)\
                                    .where(indel_model_table.c.replicate_id == replicate_id)\
                                    .where(indel_model_table.c.filter_delete == 0)
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
        #
        ##########################################################

        # else:
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
        df_out = f14_filter_codon_stop(variant_read_count_df, variant_df, genetic_table_number)

        # ##########################################################
        # #
        # # 5. Insert Filter data
        # #
        # ##########################################################
        # records = df_out.to_dict('records')
        # with engine.connect() as conn:
        #         conn.execute(filter_codon_stop_model.__table__.insert(), df_out.to_dict('records'))
        #
        # # df_out.to_sql(name='FilterCodonStop', con=engine.connect(), if_exists='replace')
        #


        ############################################
        # Write to DB
        ############################################
        records = filter_delete_df_to_dict(df_out)
        with engine.connect() as conn:
            conn.execute(filter_codon_stop_model.__table__.insert(), records)

        ##########################################################
        #
        # 6. Exit vtam if all variants delete
        #
        ##########################################################
        # Exit if no variants for analysis
        try:
            assert df_out.shape[0] == 0
        except AssertionError:
            Logger.instance().info(VTAMexception("Error: This filter has deleted all the variants"))
            sys.exit(1)


def f14_filter_codon_stop(variant_read_count_df, variant_df, genetic_table_number=5):
    """
    filter chimera
    """
    df = variant_df.copy()
    df2 = variant_read_count_df.copy()
    df2['filter_delete'] = False
    #
    df['orf1_codon_stop_nb'] = 0
    df['orf2_codon_stop_nb'] = 0
    df['orf3_codon_stop_nb'] = 0
    df['min_nb_codon_stop'] = 1
    #

    #
    orf_frame_index = 1  #  1-based
    #

    for row in df.iterrows():
        id = row[1].id
        sequence = row[1].sequence
        #
        sequence_orf1 = sequence[orf_frame_index - 1:] # get 1st orf sequence
        sequence_orf1 = sequence_orf1[0:len(sequence_orf1) - (len(sequence_orf1) % 3)] # trimming for module 3
        orf1_nb_codon_stop = str(Seq(sequence_orf1, IUPAC.unambiguous_dna).translate(
            Bio.Data.CodonTable.generic_by_id[genetic_table_number])).count('*')
        df.loc[df.id == id, 'orf1_codon_stop_nb'] = orf1_nb_codon_stop
        #
        sequence_orf2 = sequence[orf_frame_index:] # get 2nd orf sequence
        sequence_orf2 = sequence_orf2[0:len(sequence_orf2) - (len(sequence_orf2) % 3)] # trimming for module 3
        orf2_nb_codon_stop = str(Seq(sequence_orf2, IUPAC.unambiguous_dna).translate(
            Bio.Data.CodonTable.generic_by_id[genetic_table_number])).count('*')
        df.loc[df.id == id, 'orf2_codon_stop_nb'] = orf2_nb_codon_stop
        #
        #
        sequence_orf3 = sequence[orf_frame_index + 1:] # get 2nd orf sequence
        sequence_orf3 = sequence_orf3[0:len(sequence_orf3) - (len(sequence_orf3) % 3)] # trimming for module 3
        orf3_nb_codon_stop = str(Seq(sequence_orf3, IUPAC.unambiguous_dna).translate(
            Bio.Data.CodonTable.generic_by_id[genetic_table_number])).count('*')
        df.loc[df.id == id, 'orf3_codon_stop_nb'] = orf3_nb_codon_stop
        # if min_nb_codon_stop =0 so the variant is OK
        minimum = min(orf1_nb_codon_stop, orf2_nb_codon_stop, orf3_nb_codon_stop)
        if (minimum == 0) :
           df.loc[df.id == id,'min_nb_codon_stop'] = 0

    #
    #list of variant id that are Not OK
    list_variant_not_ok = df.id[df['min_nb_codon_stop'] == 1].tolist()
    df2.loc[df2.variant_id.isin(list_variant_not_ok),'filter_delete'] = 1

    return df2