from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from sqlalchemy import select
from vtam.utils.FilterCommon import FilterCommon
from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

import Bio
import pandas
import sys


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
        if not self.option("log_verbosity") is None:
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
        biosample_model = self.input_table(FilterCodonStop.__input_table_biosample)
        replicate_model = self.input_table(FilterCodonStop.__input_table_replicate)
        variant_model = self.input_table(FilterCodonStop.__input_table_Variant)
        input_filter_model = self.input_table(FilterCodonStop.__input_table_filter_indel)
        #options
        genetic_table_number = int(self.option("genetic_table_number"))
        #
        # Output table models
        output_filter_models = self.output_table(FilterCodonStop.__output_table_filter_codon_stop)


        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################

        filter_various = FilterCommon(self.__class__.__name__, engine, run_model, marker_model, biosample_model, replicate_model,
                                      input_filter_model,
                                      output_filter_models=output_filter_models)
        fastainfo_instance_list = filter_various.get_fastainfo_instance_list_with_ids(input_file_fastainfo)

        ##########################################################
        #
        # 2. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################

        filter_various.delete_output_filter_model(fastainfo_instance_list)


        ##########################################################
        #
        # 3. Select variant_read_count_model
        #
        ##########################################################

        variant_read_count_df = filter_various.get_variant_read_count_model(fastainfo_instance_list)

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
        filter_output_df = f14_filter_codon_stop(variant_read_count_df, variant_df, genetic_table_number)


        ##########################################################
        #
        # Write to DB
        #
        ##########################################################

        records = FilterCommon.filter_delete_df_to_dict(filter_output_df)
        with engine.connect() as conn:
            conn.execute(output_filter_models.__table__.insert(), records)

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