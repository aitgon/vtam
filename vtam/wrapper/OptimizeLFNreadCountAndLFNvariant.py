import inspect
import pandas

from sqlalchemy import select, literal

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.FilterLFNrunner import FilterLFNrunner
from vtam.wrapper.FilterMinReplicateNumber import f9_delete_min_replicate_number
from vtam.wrapper.OptimizeUtillities import KnownVariantAnalyzer


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
            "log_verbosity": "int",
            "log_file": "str",
        }

    def run(self, f7_lfn_delete_absolute_read_count=None):
        """
        Algorithm

        1. Read file with known variants (Mock, delete and real)
        2. Control if user variants and sequence are consistent in the database
        3. Get variant_read_count of the known variants
        4. Get variant read count of keep, delete-negative and delete-real
        5. Compute maximal lfn_read_count_threshold that keeps all 'keep' variants with the 'lfn_read_count_and_lfn_variant' algorithm
        6. Compute maximal lfn_variant_threshold that keeps all 'keep' variants with the 'lfn_read_count_and_lfn_variant' algorithm (See below)
        7. Compute nb of keep and delete variants until for different values of lfn_read_count_threshold and lfn_variant_threshold with the 'lfn_read_count_and_lfn_variant' algorithm

        Description of the 'lfn_read_count_and_lfn_variant' algorithm

        """
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
        # Input file path
        input_file_variant_known = self.input_file(OptimizeLFNreadCountAndLFNvariant.__input_file_variant_known)
        input_file_fastainfo = self.input_file(OptimizeLFNreadCountAndLFNvariant.__input_file_fastainfo)
        #
        # Input table models
        run_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_run)
        marker_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_marker)
        biosample_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_biosample)
        replicate_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_replicate)
        variant_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_variant)

        variant_read_count_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_variant_read_count)
        #
        # Output file path
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
        OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # Read "variants_optimize" to get run_id, marker_id, biosample_id, variant_id for current analysis
        #
        ##########################################################

        variant_known_df = pandas.read_csv(input_file_variant_known, sep="\t", header=0, \
                                           names=['marker_name', 'run_name', 'biosample_name', 'biosample_type',
                                                  'variant_id', 'action', 'variant_sequence', 'note'], index_col=False)

        ########################
        #
        # Control if known variants and sequence are consistent in the database
        #
        ########################

        variant_control_df = variant_known_df[['variant_id', 'variant_sequence']].drop_duplicates()
        variant_control_df = variant_control_df.loc[~variant_control_df.variant_id.isnull()]
        with engine.connect() as conn:
            for row in variant_control_df.itertuples():
                variant_id = row.variant_id
                variant_sequence = row.variant_sequence
                stmt_select = select([variant_model.__table__.c.id, variant_model.__table__.c.sequence]) \
                    .where(variant_model.__table__.c.id == variant_id) \
                    .where(variant_model.__table__.c.sequence == variant_sequence)
                if conn.execute(stmt_select).first() is None:
                    Logger.instance().error("Variant {} and its sequence are not coherent with the VTAM database"
                                            .format(variant_id))
                    exit(1)

        ##########################################################
        #
        # Read variant known df and get IDs from database
        #
        ##########################################################
        sample_instance_list = []
        for row in variant_known_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            variant_id = row.variant_id
            biosample_type = row.biosample_type
            action = row.action
            with engine.connect() as conn:
                stmt_select_variant_known = select([
                    run_model.__table__.c.id,
                    marker_model.__table__.c.id,
                    biosample_model.__table__.c.id,
                    literal(variant_id),
                    literal(biosample_type),
                    literal(action)]) \
                    .distinct() \
                    .where(run_name == run_model.__table__.c.name) \
                    .where(marker_name == marker_model.__table__.c.name) \
                    .where(biosample_name == biosample_model.__table__.c.name)

                sample_instance_list += conn.execute(stmt_select_variant_known).fetchall()
        variant_known_df = pandas.DataFrame(sample_instance_list,
                                            columns=['run_id', 'marker_id', 'biosample_id', 'variant_id',
                                                     'biosample_type', 'action'])
        # Change types to int
        variant_known_df.variant_id = pandas.to_numeric(variant_known_df.variant_id)

        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        variant_read_count_model = variant_read_count_model.__table__
        stmt_variant_filter_lfn = select([variant_read_count_model.c.marker_id,
                                          variant_read_count_model.c.run_id,
                                          variant_read_count_model.c.variant_id,
                                          variant_read_count_model.c.biosample_id,
                                          variant_read_count_model.c.replicate_id,
                                          variant_read_count_model.c.read_count])
        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_read_count_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                                                              columns=['marker_id', 'run_id', 'variant_id',
                                                                       'biosample_id', 'replicate_id', 'read_count'])

        ##########################################################
        #
        # Get keep variants, that is variants marked as keep in either mock or real biosamples
        #
        ##########################################################
        known_variant_analyzer = KnownVariantAnalyzer(variant_known_df=variant_known_df,
                                                      variant_read_count_df=variant_read_count_df)

        # These columns: run_id  marker_id  biosample_id  variant_id
        variant_keep_df = known_variant_analyzer.get_variant_keep_df()

        ##########################################################
        #
        # Get delete variants, that are not keep in mock samples
        #
        ##########################################################

        # These columns: run_id  marker_id  biosample_id  variant_id
        variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df, variant_delete_df \
            = known_variant_analyzer.get_variant_delete_df()

        ##############
        #
        # Search maximal value of lfn_read_count_threshold
        #
        ##############
        #
        #
        count_keep = 0
        count_keep_max = 0
        #
        variant_read_count_df = variant_read_count_df[
            ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count']]
        lfn_read_count_threshold_previous = 10
        #  loop over lfn_read_count_threshold
        for lfn_read_count_threshold in list(range(lfn_read_count_threshold_previous, 1001, 10)):
            Logger.instance().debug(
                "file: {}; line: {}; lfn_read_count_threshold: {} ----------------------"
                    .format(__file__, inspect.currentframe().f_lineno, lfn_read_count_threshold))

            variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(
                is_optimize_lfn_variant_replicate, variant_read_count_df, variant_keep_df,
                lfn_biosample_replicate_threshold,
                lfn_read_count_threshold, min_replicate_number, lfn_variant_or_variant_replicate_threshold)
            if count_keep > count_keep_max:
                count_keep_max = count_keep
            elif count_keep < count_keep_max:
                break  # stop when count_keep starts to decrease
            lfn_read_count_threshold_previous = lfn_read_count_threshold

        lfn_read_count_threshold_max = lfn_read_count_threshold_previous  # upper border of lfn_read_count_threshold
        #
        ##############
        #
        # Set maximal value of lfn_variant_or_variant_replicate
        #
        ##############
        #

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
                "file: {}; line: {}; lfn_variant_or_variant_replicate: {} ----------------------"
                    .format(__file__, inspect.currentframe().f_lineno, lfn_variant_or_variant_replicate))

            variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(
                is_optimize_lfn_variant_replicate,
                variant_read_count_df, variant_keep_df, lfn_biosample_replicate_threshold,
                lfn_read_count_threshold, min_replicate_number, lfn_variant_or_variant_replicate_threshold)

            if count_keep > count_keep_max:
                count_keep_max = count_keep
            elif count_keep < count_keep_max:
                break  # stop when count_keep starts to decrease
            lfn_variant_or_variant_replicate_previous = lfn_variant_or_variant_replicate

        lfn_variant_or_variant_replicate_max = lfn_variant_or_variant_replicate_previous  # upper border of lfn_read_count_threshold

        ##############
        #
        # Optimize together two parameters to previous define borders
        #
        ##############
        #
        out_lfn_variant_list = []

        variant_read_count_df = variant_read_count_df[
            ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count']].drop_duplicates(
            inplace=False)
        # lfn_read_count_threshold_default = lfn_read_count_threshold
        # lfn_variant_or_variant_replicate_default = lfn_variant_or_variant_replicate
        # while count_keep >= count_keep_max:
        # loop over lfn_read_count_threshold
        for lfn_read_count_threshold in list(range(10, lfn_read_count_threshold_max + 1, 5)):
            # loop over lfn_variant_or_variant_replicate_threshold: 0.001, 0.002, ...
            for lfn_variant_or_variant_replicate_threshold in [i / 1000 for i in range(1, int(lfn_variant_or_variant_replicate_max * 1000) + 1, 1)]:
                Logger.instance().debug(
                    "file: {}; line: {}; lfn_read_count_threshold: {}; lfn_variant_or_variant_replicate_threshold: {} ============="
                        .format(__file__, inspect.currentframe().f_lineno, lfn_read_count_threshold,
                                lfn_variant_or_variant_replicate_threshold))

                variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(
                    is_optimize_lfn_variant_replicate, variant_read_count_df, variant_keep_df,
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
                                               "variant_nb_keep": count_keep, "variant_nb_delete": count_delete}
                    out_lfn_variant_list.append(out_lfn_variant_row_dic)

                if count_keep > count_keep_max:
                    count_keep_max = count_keep

        ##########################################################
        #
        # Write TSV file
        #
        ##########################################################
        if not is_optimize_lfn_variant_replicate:  # optimize lfn variant
            out_lfn_variant_or_variant_replicate_df = pandas.DataFrame(out_lfn_variant_list, columns=['variant_nb_keep', 'variant_nb_delete',
                                                                                 'lfn_read_count_threshold',
                                                                                 'lfn_variant_or_variant_replicate_threshold'])
            out_lfn_variant_or_variant_replicate_df.columns = ['variant_nb_keep', 'variant_nb_delete',
                                                                                 'lfn_read_count_threshold',
                                                                                 'lfn_variant_threshold']
        else:  # optimize lfn variant replicate
            out_lfn_variant_or_variant_replicate_df = pandas.DataFrame(out_lfn_variant_list, columns=['variant_nb_keep', 'variant_nb_delete',
                                                                                 'lfn_read_count_threshold',
                                                                                 'lfn_variant_or_variant_replicate_threshold'])
            out_lfn_variant_or_variant_replicate_df.columns = ['variant_nb_keep', 'variant_nb_delete',
                                                                                 'lfn_read_count_threshold',
                                                                                 'lfn_variant_replicate_threshold']

        out_lfn_variant_or_variant_replicate_df.sort_values(by=["variant_nb_keep", "variant_nb_delete"],
                                                            ascending=[False, True], inplace=True)

        out_lfn_variant_or_variant_replicate_df.to_csv(output_file_optimize_lfn_tsv, header=True, sep='\t', index=False)

        ##########################################################
        #
        # LFN variant specific threshold
        #
        ##########################################################
        # variant_delete_df = pandas.concat([variant_read_count_delete_negative_df, variant_delete_real_df], axis=0)
        # lfn_variant_specific_threshold_df = variant_delete_df.copy()

        ##########################################################
        #
        # Variant delete: variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df
        #
        ##########################################################
        variant_read_count_delete_mock_df = variant_read_count_df.merge(variant_delete_mock_df, on=['run_id',
                                                                                                    'marker_id',
                                                                                                    'biosample_id',
                                                                                                    'variant_id'])
        # variant_read_count_delete_mock_df['biosample_type'] = 'mock'
        variant_read_count_negative_df = variant_read_count_df.merge(variant_delete_negative_df, on=['run_id',
                                                                                                     'marker_id',
                                                                                                     'biosample_id',
                                                                                                     'variant_id'])
        # variant_read_count_delete_mock_df['biosample_type'] = 'negative'
        variant_read_count_delete_real_df = variant_read_count_df.merge(variant_delete_real_df, on=['run_id',
                                                                                                    'marker_id',
                                                                                                    'biosample_id',
                                                                                                    'variant_id'])
        # variant_read_count_delete_mock_df['biosample_type'] = 'real'
        # variant_read_count_delete_mock_df['action'] = 'delete'
        variant_read_count_delete_df \
            = pandas.concat(
            [variant_read_count_delete_mock_df, variant_read_count_negative_df, variant_read_count_delete_real_df])

        if not is_optimize_lfn_variant_replicate:  # optimize lfn variant
            N_i_df = variant_read_count_delete_df[['run_id', 'marker_id', 'variant_id', 'read_count']] \
                .groupby(by=['run_id', 'marker_id', 'variant_id']) \
                .sum().reset_index()
            N_i_df = N_i_df.rename(columns={'read_count': 'N_i'})
            N_i_df.drop_duplicates(inplace=True)
            lfn_variant_or_variant_replicate_specific_threshold_df = variant_read_count_delete_df.merge(N_i_df,
                                                                                                        on=['run_id',
                                                                                                            'marker_id',
                                                                                                            'variant_id'])
            del (N_i_df)
            lfn_variant_or_variant_replicate_specific_threshold_df[
                'lfn_variant_threshold'] = lfn_variant_or_variant_replicate_specific_threshold_df.read_count \
                                           / lfn_variant_or_variant_replicate_specific_threshold_df.N_i
            # TODO 20191006 Discuss with Emese with ascending=True or ascending=False
            lfn_variant_or_variant_replicate_specific_threshold_df.sort_values(by='lfn_variant_threshold',
                                                                               ascending=True, inplace=True)
            lfn_variant_or_variant_replicate_specific_threshold_df.drop_duplicates('variant_id', keep='first',
                                                                                   inplace=True)
            lfn_variant_or_variant_replicate_specific_threshold_df = (
            lfn_variant_or_variant_replicate_specific_threshold_df[
                ['run_id', 'marker_id', 'variant_id', 'read_count', 'N_i', 'lfn_variant_threshold']]).drop_duplicates(
                inplace=False)
        else:  # optimize lfn variant replicate
            N_ik_df = variant_read_count_delete_df[['run_id', 'marker_id', 'variant_id', 'replicate_id', 'read_count']] \
                .groupby(by=['run_id', 'marker_id', 'variant_id', 'replicate_id']) \
                .sum().reset_index()
            N_ik_df = N_ik_df.rename(columns={'read_count': 'N_ik'})
            N_ik_df.drop_duplicates(inplace=True)
            lfn_variant_or_variant_replicate_specific_threshold_df = variant_read_count_delete_df.merge(N_ik_df,
                                                                                                        on=['run_id',
                                                                                                            'marker_id',
                                                                                                            'variant_id',
                                                                                                            'replicate_id'])
            del (N_ik_df)
            lfn_variant_or_variant_replicate_specific_threshold_df[
                'lfn_variant_replicate_threshold'] = lfn_variant_or_variant_replicate_specific_threshold_df.read_count \
                                                     / lfn_variant_or_variant_replicate_specific_threshold_df.N_ik
            # TODO Discuss with Emese with ascending=True or ascending=False
            lfn_variant_or_variant_replicate_specific_threshold_df.sort_values(by='lfn_variant_replicate_threshold',
                                                                               ascending=True, inplace=True)
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

        if not is_optimize_lfn_variant_replicate:  # optimize lfn variant
            lfn_variant_or_variant_replicate_specific_threshold_df = \
            lfn_variant_or_variant_replicate_specific_threshold_df[
                ['run_id', 'marker_id', 'variant_id', 'read_count', 'N_i', 'lfn_variant_threshold',
                 'biosample_type']].groupby(
                by=['run_id', 'marker_id', 'variant_id', 'read_count', 'N_i', 'lfn_variant_threshold'])[
                'biosample_type'].apply(lambda x: ','.join(set(x))).reset_index()
        else:  # optimize lfn variant replicate
            lfn_variant_or_variant_replicate_specific_threshold_df = \
            lfn_variant_or_variant_replicate_specific_threshold_df[
                ['run_id', 'marker_id', 'variant_id', 'replicate_id', 'read_count', 'N_ik', 'lfn_variant_replicate_threshold',
                 'biosample_type']].groupby(
                by=['run_id', 'marker_id', 'variant_id', 'replicate_id', 'read_count', 'N_ik', 'lfn_variant_replicate_threshold'])[
                'biosample_type'].apply(lambda x: ','.join(set(x))).reset_index()
            lfn_variant_or_variant_replicate_specific_threshold_df = lfn_variant_or_variant_replicate_specific_threshold_df.rename(
                columns={'read_count': 'N_ijk_max'})
            lfn_variant_or_variant_replicate_specific_threshold_df.to_csv(
                output_file_lfn_variant_specific_threshold_tsv, header=True, sep='\t', index=False)


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
