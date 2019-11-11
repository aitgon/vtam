import inspect
import pandas
import sqlalchemy

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from vtam.utils.FastaInformation import FastaInformation
from vtam.utils.Logger import Logger
from vtam.utils.VariantReadCountDF import VariantReadCountDF
from vtam.utils.OptionManager import OptionManager
from vtam.utils.FilterLFNrunner import FilterLFNrunner
from vtam.utils.VariantKnown import VariantKnown
from vtam.wrapper.FilterMinReplicateNumber import f9_delete_min_replicate_number


class OptimizeLFNreadCountAndLFNvariant(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.OptimizeLFNreadCountAndLFNvariant"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
    __input_file_variant_known = "variant_known"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn_read_count_and_lfn_variant = "optimize_lfn_read_count_and_lfn_variant"
    __output_file_optimize_lfn_variant_specific = "optimize_lfn_variant_specific"

    def specify_input_file(self):
        return [
            OptimizeLFNreadCountAndLFNvariant.__input_file_fastainfo,
            OptimizeLFNreadCountAndLFNvariant.__input_file_variant_known,
        ]

    def specify_input_table(self):
        return [
            OptimizeLFNreadCountAndLFNvariant.__input_table_marker,
            OptimizeLFNreadCountAndLFNvariant.__input_table_run,
            OptimizeLFNreadCountAndLFNvariant.__input_table_biosample,
            OptimizeLFNreadCountAndLFNvariant.__input_table_replicate,
            OptimizeLFNreadCountAndLFNvariant.__input_table_variant,
            OptimizeLFNreadCountAndLFNvariant.__input_table_variant_read_count,
        ]

    def specify_output_file(self):
        return [
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_read_count_and_lfn_variant,
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_variant_specific,
        ]

    def specify_params(self):
        return {
            "is_optimize_lfn_variant_replicate": "int",
            "lfn_variant_or_variant_replicate_threshold": "float",
            "lfn_biosample_replicate_threshold": "float",
            "lfn_read_count_threshold": "float",
            "min_replicate_number": "int",
            # "log_verbosity": "int",
            # "log_file": "str",
        }

    def run(self, f7_lfn_delete_absolute_read_count=None):
        """
        Algorithm (Updated Oct 13, 2019)

        1. Read file with known variants (Mock/tolerate, delete and real)
        2. Control if user variants and sequence are consistent in the database
        3. Get variant_read_count of this run-marker-biosample-replicate experiment
        5. Compute maximal lfn_read_count_threshold that keeps all 'keep' variants with the 'lfn_read_count_and_lfn_variant' algorithm
        6. Compute maximal lfn_variant_threshold that keeps all 'keep' variants with the 'lfn_read_count_and_lfn_variant' algorithm (See below)
        7. Loop between default and lfn_read_count_threshold and lfn_read_count_and_lfn_variant parameters
            7.1 Compute number of keep variants. Should be always maximal.
            7.2 Compute number of delete variants Should decrease.
        8. Compute variant(-replicate) specific threshold for delete variants
            8.1 For each variant i (Or variant-replicate i-k ),
                get N_ijk_max and use it to computer variant specific threshold

        Description of the 'lfn_read_count_and_lfn_variant' algorithm

        1. Remove if does not pass these filter
            1.1 Filter lfn_variant (Or lfn_variant_replicate)
            1.2 Filter lfn_biosample_replicate
            1.3 Filter absolute read count
        2. Filter if not min replicate number

        """
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        # if not self.option("log_verbosity") is None:
        #     OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file output
        variant_known_tsv = self.input_file(OptimizeLFNreadCountAndLFNvariant.__input_file_variant_known)
        fasta_info_tsv = self.input_file(OptimizeLFNreadCountAndLFNvariant.__input_file_fastainfo)
        #
        # Input table models
        run_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_run)
        marker_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_marker)
        biosample_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_biosample)
        replicate_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_replicate)
        variant_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_variant)

        variant_read_count_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_variant_read_count)
        #
        # Output file output
        output_file_optimize_lfn_tsv = self.output_file(
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_read_count_and_lfn_variant)
        output_file_lfn_variant_specific_threshold_tsv = self.output_file(
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_variant_specific)
        #
        # Options
        is_optimize_lfn_variant_replicate = bool(int(self.option("is_optimize_lfn_variant_replicate")))
        min_replicate_number = self.option("min_replicate_number")
        lfn_biosample_replicate_threshold = self.option("lfn_biosample_replicate_threshold")
        lfn_read_count_threshold = self.option("lfn_read_count_threshold")
        lfn_variant_or_variant_replicate_threshold = self.option("lfn_variant_or_variant_replicate_threshold")

        ################################################################################################################
        #
        # Read user known variant information and verify consistency with DB
        #
        ################################################################################################################

        variant_known = VariantKnown(variant_known_tsv, fasta_info_tsv, engine, variant_model, run_model, marker_model,
                                     biosample_model, replicate_model)
        run_marker_biosample_variant_keep_df = variant_known.get_run_marker_biosample_variant_keep_df()
        # variant_known_ids_df = variant_known.variant_known_ids_df

        ################################################################################################################
        #
        # Get variant read count df
        #
        ################################################################################################################

        fasta_info = FastaInformation(fasta_info_tsv, engine, run_model, marker_model, biosample_model, replicate_model)
        variant_read_count_df = fasta_info.get_variant_read_count_df(variant_read_count_model)

        ################################################################################################
        #
        # Get delete variants, that are not keep in mock samples
        #
        ################################################################################################

        # These columns: run_id  marker_id  biosample_id  variant_id
        variant_delete_df, variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df \
            = variant_known.get_run_marker_biosample_variant_delete_df(variant_read_count_df)

        ################################################################################################
        #
        # Search maximal value of lfn_read_count_threshold
        #
        ################################################################################################
        #
        #
        Logger.instance().debug("Search maximal value of lfn_variant_or_variant_replicate")

        count_keep = 0
        count_keep_max = 0
        #
        variant_read_count_df = variant_read_count_df[
            ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count']]
        lfn_read_count_threshold_previous = 10
        #  loop over lfn_read_count_threshold
        for lfn_read_count_threshold in list(range(lfn_read_count_threshold_previous, 1001, 10)):
            Logger.instance().debug(
                "file: {}; line: {}; lfn_read_count_threshold: {}; count_keep_max: {}; count_keep: {} ----------------------"
                    .format(__file__, inspect.currentframe().f_lineno, lfn_read_count_threshold, count_keep_max, count_keep))

            variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(
                is_optimize_lfn_variant_replicate, variant_read_count_df, run_marker_biosample_variant_keep_df,
                lfn_biosample_replicate_threshold,
                lfn_read_count_threshold, min_replicate_number, lfn_variant_or_variant_replicate_threshold)
            if count_keep > count_keep_max:
                count_keep_max = count_keep
            elif count_keep < count_keep_max:
                break  # stop when count_keep starts to decrease
            lfn_read_count_threshold_previous = lfn_read_count_threshold

        lfn_read_count_threshold_max = lfn_read_count_threshold_previous  # upper border of lfn_read_count_threshold

        ################################################################################################
        #
        # Search maximal value of lfn_variant_or_variant_replicate
        #
        ################################################################################################
        #

        Logger.instance().debug("Search maximal value of lfn_variant_or_variant_replicate")
        count_keep = 0
        count_keep_max = 0
        #
        variant_read_count_df = variant_read_count_df[
            ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count']].drop_duplicates(
            inplace=False)
        lfn_variant_or_variant_replicate_previous = 0.001  # is divided by 1000
        # loop over lfn_variant_or_variant_replicate: 0.001, 0.002, ...
        for lfn_variant_or_variant_replicate in [i / 1000 for i in range(int(lfn_variant_or_variant_replicate_previous * 1000), 101, 1)]:
            Logger.instance().debug(
                "file: {}; line: {}; lfn_variant_or_variant_replicate_threshold: {}; count_keep_max: {}; count_keep: {} ----------------------"
                    .format(__file__, inspect.currentframe().f_lineno, lfn_variant_or_variant_replicate, count_keep_max, count_keep))

            variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(
                is_optimize_lfn_variant_replicate,
                variant_read_count_df, run_marker_biosample_variant_keep_df, lfn_biosample_replicate_threshold,
                lfn_read_count_threshold, min_replicate_number, lfn_variant_or_variant_replicate_threshold)

            if count_keep > count_keep_max:
                count_keep_max = count_keep
            elif count_keep < count_keep_max:
                break  # stop when count_keep starts to decrease
            lfn_variant_or_variant_replicate_previous = lfn_variant_or_variant_replicate

        lfn_variant_or_variant_replicate_max = lfn_variant_or_variant_replicate_previous  # upper border of lfn_read_count_threshold

        ################################################################################################
        #
        # Optimize together two parameters to previous define borders
        #
        ################################################################################################
        #
        out_lfn_variant_list = []

        variant_read_count_df = variant_read_count_df[
            ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count']].drop_duplicates(
            inplace=False)
        # loop over lfn_read_count_threshold
        for lfn_read_count_threshold in list(range(10, lfn_read_count_threshold_max + 1, 5)):
            # loop over lfn_variant_or_variant_replicate_threshold: 0.001, 0.002, ...
            for lfn_variant_or_variant_replicate_threshold in [i / 1000 for i in range(1, int(lfn_variant_or_variant_replicate_max * 1000) + 1, 1)]:
                Logger.instance().debug(
                    "file: {}; line: {}; lfn_read_count_threshold: {}; lfn_variant_or_variant_replicate_threshold: {}; "
                    "count_keep_max: {}; count_keep: {} ============="
                        .format(__file__, inspect.currentframe().f_lineno, lfn_read_count_threshold,
                                lfn_variant_or_variant_replicate_threshold, count_keep_max, count_keep))

                variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(
                    is_optimize_lfn_variant_replicate, variant_read_count_df, run_marker_biosample_variant_keep_df,
                    lfn_biosample_replicate_threshold,
                    lfn_read_count_threshold, min_replicate_number, lfn_variant_or_variant_replicate_threshold)

                ##########################################################
                #
                # Count delete
                #
                ##########################################################
                variant_read_count_remained_delete_negative_df = variant_read_count_remained_df.merge(variant_delete_df,
                                                                                                      on=['run_id',
                                                                                                          'marker_id',
                                                                                                          'biosample_id',
                                                                                                          'variant_id']).drop_duplicates(
                    inplace=False)
                count_delete = variant_read_count_remained_delete_negative_df.shape[0]

                ##########################################################
                #
                # Store results
                #
                ##########################################################
                if count_keep >= count_keep_max:  # Store results if count_keep maximal
                    out_lfn_variant_row_dic = {"lfn_variant_or_variant_replicate_threshold": lfn_variant_or_variant_replicate_threshold,
                                               "lfn_read_count_threshold": lfn_read_count_threshold,
                                               "occurence_nb_keep": count_keep, "occurence_nb_delete": count_delete}
                    out_lfn_variant_list.append(out_lfn_variant_row_dic)

                if count_keep > count_keep_max:
                    count_keep_max = count_keep



        ##########################################################
        #
        # Write TSV file
        #
        ##########################################################
        # From list of dics to df
        out_lfn_variant_or_variant_replicate_df = pandas.DataFrame(out_lfn_variant_list)
        # List of columns in order
        column_names = ['occurence_nb_keep', 'occurence_nb_delete', 'lfn_read_count_threshold', 'lfn_variant_or_variant_replicate_threshold']
        # Reorder columns
        out_lfn_variant_or_variant_replicate_df = out_lfn_variant_or_variant_replicate_df[column_names]
        # Sort columns
        out_lfn_variant_or_variant_replicate_df.sort_values(by=column_names,
                                                            ascending=[False, True, True, True], inplace=True)
        # Rename columns depending on whether this is optimize_lfn_variant or is_optimize_lfn_variant_replicate
        if not is_optimize_lfn_variant_replicate:  # optimize lfn variant
            out_lfn_variant_or_variant_replicate_df = out_lfn_variant_or_variant_replicate_df\
                .rename(columns={'lfn_variant_or_variant_replicate_threshold': 'lfn_variant_threshold'})
        else:  # optimize lfn variant replicate
            out_lfn_variant_or_variant_replicate_df = out_lfn_variant_or_variant_replicate_df\
                .rename(columns={'lfn_variant_or_variant_replicate_threshold': 'lfn_variant_replicate_threshold'})

        out_lfn_variant_or_variant_replicate_df.to_csv(output_file_optimize_lfn_tsv, header=True, sep='\t', index=False,
                                                       float_format='%.6f')


        ##########################################################
        #
        # Variant delete: variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df
        #
        ##########################################################
        variant_read_count_delete_df = variant_read_count_df.merge(variant_delete_df, on=['run_id',
                                                                                                    'marker_id',
                                                                                                    'biosample_id',
                                                                                                    'variant_id'])
        variant_read_count_instance = VariantReadCountDF(variant_read_count_delete_df)

        if not is_optimize_lfn_variant_replicate:  # optimize lfn variant
            N_i_df = variant_read_count_instance.get_N_i_df()

            lfn_variant_or_variant_replicate_specific_threshold_df = variant_read_count_delete_df.merge(N_i_df,
                                                                                                        on=['run_id',
                                                                                                            'marker_id',
                                                                                                            'variant_id'])
            del (N_i_df)
            lfn_variant_or_variant_replicate_specific_threshold_df[
                'lfn_variant_threshold'] = lfn_variant_or_variant_replicate_specific_threshold_df.read_count \
                                           / lfn_variant_or_variant_replicate_specific_threshold_df.N_i
            lfn_variant_or_variant_replicate_specific_threshold_df.sort_values(by='lfn_variant_threshold',
                                                                               ascending=False, inplace=True)
            lfn_variant_or_variant_replicate_specific_threshold_df.drop_duplicates('variant_id', keep='first',
                                                                                   inplace=True)
            lfn_variant_or_variant_replicate_specific_threshold_df = (
            lfn_variant_or_variant_replicate_specific_threshold_df[
                ['run_id', 'marker_id', 'variant_id', 'read_count', 'N_i', 'lfn_variant_threshold']]).drop_duplicates(
                inplace=False)
        else:  # optimize lfn variant replicate
            N_ik_df = variant_read_count_instance.get_N_ik_df()
            lfn_variant_or_variant_replicate_specific_threshold_df = variant_read_count_delete_df.merge(N_ik_df,
                                                                                                        on=['run_id',
                                                                                                            'marker_id',
                                                                                                            'variant_id',
                                                                                                            'replicate_id'])
            lfn_variant_or_variant_replicate_specific_threshold_df[
                'lfn_variant_replicate_threshold'] = lfn_variant_or_variant_replicate_specific_threshold_df.read_count \
                                                     / lfn_variant_or_variant_replicate_specific_threshold_df.N_ik
            lfn_variant_or_variant_replicate_specific_threshold_df.sort_values(by='lfn_variant_replicate_threshold',
                                                                               ascending=False, inplace=True)
            lfn_variant_or_variant_replicate_specific_threshold_df.drop_duplicates(['variant_id', 'replicate_id'],
                                                                                   keep='first', inplace=True)
            lfn_variant_or_variant_replicate_specific_threshold_df = (
            lfn_variant_or_variant_replicate_specific_threshold_df[
                ['run_id', 'marker_id', 'variant_id', 'replicate_id', 'read_count', 'N_ik',
                 'lfn_variant_replicate_threshold']]).drop_duplicates(inplace=False)

        ##########################################################
        #
        # Annotate with delete not in mock
        #
        ##########################################################
        lfn_variant_delete_mock_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df.merge(
            variant_delete_mock_df,
            on=['run_id',
                'marker_id',
                'variant_id'])
        lfn_variant_delete_mock_specific_threshold_df['biosample_type'] = 'mock'
        lfn_variant_delete_mock_specific_threshold_df['action'] = ''

        ##########################################################
        #
        # Annotate with delete not in mock
        #
        ##########################################################
        lfn_variant_delete_negative_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df.merge(
            variant_delete_negative_df,
            on=['run_id',
                'marker_id',
                'variant_id'])
        lfn_variant_delete_negative_specific_threshold_df['biosample_type'] = 'negative'
        lfn_variant_delete_negative_specific_threshold_df['action'] = ''

        ##########################################################
        #
        # Annotate with delete in real
        #
        ##########################################################
        lfn_variant_delete_real_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df.merge(
            variant_delete_real_df,
            on=['run_id',
                'marker_id',
                'variant_id'])
        lfn_variant_delete_real_specific_threshold_df['biosample_type'] = 'real'
        lfn_variant_delete_real_specific_threshold_df['action'] = 'delete'

        lfn_variant_or_variant_replicate_specific_threshold_df = pandas.concat(
            [lfn_variant_delete_mock_specific_threshold_df, lfn_variant_delete_negative_specific_threshold_df,
             lfn_variant_delete_real_specific_threshold_df])

        lfn_variant_or_variant_replicate_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df.rename(
            columns={'read_count': 'N_ijk_max'})
        if not is_optimize_lfn_variant_replicate:  # optimize lfn variant
            lfn_variant_or_variant_replicate_specific_threshold_df.groupby(
                by=['run_id', 'marker_id', 'variant_id', 'N_ijk_max', 'N_i', 'lfn_variant_threshold']).agg(
                {'biosample_type': lambda x: ','.join(sorted(set(x)))}).reset_index()
        else:  # optimize lfn variant replicate
            lfn_variant_or_variant_replicate_specific_threshold_df.groupby(
                by=['run_id', 'marker_id', 'variant_id', 'N_ijk_max', 'N_ik', 'lfn_variant_replicate_threshold']).agg(
                {'biosample_type': lambda x: ','.join(sorted(set(x)))}).reset_index()

        ##########################################################
        #
        # Convert run_id, marker_id and biosample_id to their names
        #
        ##########################################################

        with engine.connect() as conn:
            run_id_to_name = conn.execute(
                sqlalchemy.select([run_model.__table__.c.id, run_model.__table__.c.name])).fetchall()
            marker_id_to_name = conn.execute(
                sqlalchemy.select([marker_model.__table__.c.id, marker_model.__table__.c.name])).fetchall()
            if is_optimize_lfn_variant_replicate:  # optimize lfn variant
                replicate_id_to_name = conn.execute(
                    sqlalchemy.select([replicate_model.__table__.c.id, replicate_model.__table__.c.name])).fetchall()

        run_id_to_name_df = pandas.DataFrame.from_records(data=run_id_to_name, columns=['run_id', 'run_name'])
        lfn_variant_or_variant_replicate_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df.merge(run_id_to_name_df, on='run_id')

        marker_id_to_name_df = pandas.DataFrame.from_records(data=marker_id_to_name, columns=['marker_id', 'marker_name'])
        lfn_variant_or_variant_replicate_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df.merge(marker_id_to_name_df, on='marker_id')


        if not is_optimize_lfn_variant_replicate:  # optimize lfn variant
            lfn_variant_or_variant_replicate_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df[
                ['run_name', 'marker_name', 'variant_id', 'N_ijk_max', 'N_i', 'lfn_variant_threshold', 'biosample_type']]
            column_names = ['run_name', 'marker_name',
                               'variant_id', 'N_ijk_max', 'N_i', 'lfn_variant_threshold', 'biosample_type']
            lfn_variant_or_variant_replicate_specific_threshold_df.sort_values(by=column_names, inplace=True)
        else:  # optimize lfn variant replicate
            replicate_id_to_name_df = pandas.DataFrame.from_records(data=replicate_id_to_name, columns=['replicate_id', 'replicate_name'])
            lfn_variant_or_variant_replicate_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df.merge(replicate_id_to_name_df, on='replicate_id')
            column_names = ['run_name', 'marker_name', 'variant_id', 'replicate_name',
                               'N_ijk_max', 'N_ik', 'lfn_variant_replicate_threshold', 'biosample_type']
            lfn_variant_or_variant_replicate_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df[column_names]
            lfn_variant_or_variant_replicate_specific_threshold_df.sort_values(by=column_names, inplace=True)

        #
        # Print only when N_ijk_max larger than default lfn_read_count_threshold
        lfn_variant_or_variant_replicate_specific_threshold_df \
            = lfn_variant_or_variant_replicate_specific_threshold_df.loc[
            lfn_variant_or_variant_replicate_specific_threshold_df.N_ijk_max > lfn_read_count_threshold]

        # Aggregate biosample types
        lfn_variant_or_variant_replicate_specific_threshold_df \
            = lfn_variant_or_variant_replicate_specific_threshold_df.groupby(
            column_names[:-1]).agg({'biosample_type': lambda x: ','.join(sorted(set(x)))}).reset_index()

        lfn_variant_or_variant_replicate_specific_threshold_df.to_csv(
            output_file_lfn_variant_specific_threshold_tsv, header=True, sep='\t', index=False, float_format='%.6f')


def lfn_read_count_and_lfn_variant(is_optimize_lfn_variant_replicate, variant_read_count_df, variant_keep_df,
                                   lfn_biosample_replicate_threshold,
                                   lfn_read_count_threshold, min_replicate_number,
                                   lfn_variant_or_variant_replicate_threshold):
    lfn_filter_runner = FilterLFNrunner(variant_read_count_df)

    ###################
    #
    # Filter lfn_variant
    #
    ####################

    if not is_optimize_lfn_variant_replicate:  # optimize lfn variant replicate
        lfn_filter_runner.f2_f4_lfn_delete_variant(lfn_variant_or_variant_replicate_threshold)
    else:  # optimize lfn variant replicate
        lfn_filter_runner.f3_f5_lfn_delete_variant_replicate(lfn_variant_or_variant_replicate_threshold)

    ###################
    #
    # Filter lfn_biosample_replicate
    #
    ####################

    lfn_filter_runner.f6_lfn_delete_biosample_replicate(lfn_biosample_replicate_threshold)

    ###################
    #
    # Filter absolute read count
    #
    ####################

    lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)

    ###################
    #
    # f8_lfn_delete_do_not_pass_all_filters
    #
    ####################

    lfn_filter_runner.f8_lfn_delete_do_not_pass_all_filters()

    variant_read_count_remained_df = lfn_filter_runner.variant_read_count_filter_delete_df

    variant_read_count_remained_df = variant_read_count_remained_df.loc[
        (variant_read_count_remained_df.filter_id == 8) &
        (variant_read_count_remained_df.filter_delete == 0)]

    ##########################################################
    #
    # f9_delete_min_replicate_number
    #
    ##########################################################

    variant_read_count_remained_df = f9_delete_min_replicate_number(variant_read_count_remained_df,
                                                                    min_replicate_number)
    variant_read_count_remained_df = variant_read_count_remained_df.loc[
        (variant_read_count_remained_df.filter_delete == 0)]
    variant_read_count_remained_df.drop('filter_delete', axis=1, inplace=True)

    ##########################################################
    #
    # Count keep
    #
    ##########################################################

    variant_read_count_remained_df = variant_read_count_remained_df[
        ['run_id', 'marker_id', 'biosample_id', 'variant_id']]
    variant_read_count_remained_df.drop_duplicates(inplace=True)

    variant_read_count_remained_keep_df = variant_read_count_remained_df.merge(variant_keep_df,
                                                                               on=['run_id', 'marker_id',
                                                                                   'biosample_id', 'variant_id'])
    variant_read_count_remained_keep_df = variant_read_count_remained_keep_df[
        ['run_id', 'marker_id', 'variant_id', 'biosample_id']].drop_duplicates()
    count_keep = variant_read_count_remained_keep_df.shape[0]

    # Delete object
    del (lfn_filter_runner)
    return variant_read_count_remained_df, count_keep
