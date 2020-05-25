import pandas
from vtam.utils.FilterLFNreplicateRemainRunner import FilterLFNreplicateRemainRunner

from vtam.utils.FilterLFNRunner import FilterLFNrunner
from vtam.utils.FilterMinReplicateNumberRunner import FilterMinReplicateNumberRunner
from vtam.utils.Logger import Logger
from vtam.utils.VariantReadCountLikeDF import VariantReadCountLikeDF


class OptimizeLFNreadCountAndLFNvariantRunner:

    def __init__(self, nijk_df, known_occurrences_df):

        self.nijk_df = nijk_df
        self.known_occurrences_df = known_occurrences_df

    def get_optimize_df(self, lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff,
                        min_replicate_number):

        filter_kwargs = {"lfn_ni_cutoff": lfn_ni_cutoff,
                         "lfn_nik_cutoff": lfn_nik_cutoff,
                         "lfn_njk_cutoff": lfn_njk_cutoff,
                         "lfn_nijk_cutoff": lfn_nijk_cutoff,
                         'min_replicate_number': min_replicate_number,
                      }

        lfn_ni_nik_cutoff = lfn_ni_cutoff  # current max
        which_filter = 'lfn_variant_cutoff'
        if not (lfn_nik_cutoff is None):
            lfn_ni_nik_cutoff = lfn_nik_cutoff  # current max
            which_filter = 'lfn_variant_replicate_cutoff'

        ############################################################################################
        #
        # Group and run_name this genetic_code by run_name/marker_name combination
        #  Loop by run_name/marker_name
        #
        ############################################################################################

        out_optimize_df = pandas.DataFrame()
        out_optimize2_df = pandas.DataFrame()

        for row in self.known_occurrences_df[['run_id', 'marker_id']].drop_duplicates().itertuples():

            run_id = row.run_id
            marker_id = row.marker_id
            known_occurrs_run_marker_df = self.known_occurrences_df.loc[
                (self.known_occurrences_df.run_id == run_id) & (
                            self.known_occurrences_df.marker_id == marker_id),]
            nijk_run_marker_df = self.nijk_df.loc[
                (self.nijk_df.run_id == run_id) & (self.nijk_df.marker_id == marker_id),]
            nijk_run_marker_df = nijk_run_marker_df[
                ['run_id', 'marker_id', 'biosample_id', 'replicate', 'variant_id',
                 'read_count']].drop_duplicates(inplace=False)

            ########################################################################################
            #
            # Keep and delete variants
            #
            ########################################################################################

            keep_run_marker_biosample_variant_df = known_occurrs_run_marker_df.loc[
                known_occurrs_run_marker_df.action == 'keep', ['run_id', 'marker_id',
                                                               'biosample_id', 'variant_id']]
            count_keep_max = keep_run_marker_biosample_variant_df.shape[
                0]  # count_keep_max taken from known_occurrences.tsv

            delete_run_marker_biosample_variant_df = known_occurrs_run_marker_df.loc[
                self.known_occurrences_df.action == 'delete', ['run_id', 'marker_id', 'biosample_id',
                                                          'variant_id']]

            ########################################################################################
            #
            # Search maximal value of lfn_nijk_cutoff
            #
            ########################################################################################

            Logger.instance().debug("Searching upper limit of lfn_nijk_cutoff")

            lfn_nijk_cutoff_local_max = lfn_nijk_cutoff  # current max
            lfn_nijk_cutoff_global_max = 100  # max value
            filter_kwargs_local = filter_kwargs.copy()

            #  loop over lfn_nijk_cutoff
            for lfn_nijk_cutoff_item in list(
                    range(lfn_nijk_cutoff, lfn_nijk_cutoff_global_max + 1, 10)):
                filter_kwargs_local['lfn_nijk_cutoff'] = lfn_nijk_cutoff_item

                nijk_remain_df = FilterLFNreplicateRemainRunner(
                    nijk_df=nijk_run_marker_df).get_nijk_remain_df(**filter_kwargs_local)

                count_keep = nijk_remain_df.merge(
                    keep_run_marker_biosample_variant_df,
                    on=['run_id', 'marker_id', 'biosample_id', 'variant_id']).drop_duplicates(
                    inplace=False).shape[0]

                Logger.instance().debug(
                    "lfn_nijk_cutoff: {} (Max {}); count_keep: {} (Max {})"
                        .format(lfn_nijk_cutoff_item, lfn_nijk_cutoff_global_max, count_keep,
                                count_keep_max))

                if count_keep < count_keep_max:
                    break  # stops when count_keep decreases below count_keep_max

                lfn_nijk_cutoff_local_max = lfn_nijk_cutoff_item  # new max lfn_nijk_cutoff

            Logger.instance().debug(
                "Upper limit of lfn_nijk_cutoff: {}".format(lfn_nijk_cutoff_local_max))

            ########################################################################################
            #
            # Search maximal value of lfn_ni_nik_cutoff
            #
            ########################################################################################

            Logger.instance().debug(
                "Searching upper limit of lfn_variant_threshold or lfn_variant_replicate_threshold")

            lfn_ni_njk_cutoff_global_max = 0.05  # max value
            filter_kwargs_local = filter_kwargs.copy()

            # loop over lfn_variant_or_variant_replicate: 0.001, 0.002, ...
            lfn_ni_nik_cutoff_list = [
                i / 1000 for i in range(int(lfn_ni_nik_cutoff * 1000), int(lfn_ni_njk_cutoff_global_max * 1000) + 1, 1)]
            for lfn_ni_nik_cutoff_item in lfn_ni_nik_cutoff_list:

                if lfn_nik_cutoff is None:
                    filter_kwargs_local['lfn_ni_cutoff'] = lfn_ni_nik_cutoff_item
                else:
                    filter_kwargs_local['lfn_nik_cutoff'] = lfn_ni_nik_cutoff_item

                nijk_remain_df = FilterLFNreplicateRemainRunner(
                    nijk_df=nijk_run_marker_df).get_nijk_remain_df(**filter_kwargs_local)

                count_keep = nijk_remain_df.merge(
                    keep_run_marker_biosample_variant_df,
                    on=['run_id', 'marker_id', 'biosample_id', 'variant_id']).drop_duplicates(
                    inplace=False).shape[0]

                Logger.instance().debug(
                    "{}: {} (Max {}); count_keep: {} (Max {})"
                        .format(which_filter, lfn_ni_nik_cutoff_item,
                                lfn_ni_njk_cutoff_global_max,
                                count_keep, count_keep_max))

                if count_keep < count_keep_max:
                    break  # stops when count_keep decreases below count_keep_max

                lfn_ni_njk_cutoff_local_max = lfn_ni_nik_cutoff_item


            ########################################################################################
            #
            # Optimize together two parameters to previous define borders
            #
            ########################################################################################

            out_lfn_variant_list = []

            nb_points = 10
            lfn_nijk_cutoff_list = [*range(lfn_nijk_cutoff, lfn_nijk_cutoff_local_max + 1, int(
                (lfn_nijk_cutoff_local_max - lfn_nijk_cutoff) / nb_points))]
            lfn_ni_nik_cutoff_list = [i / 1000 for i in range(int(
                lfn_ni_nik_cutoff * 1000), int(
                lfn_ni_njk_cutoff_local_max * 1000) + 1, int(
                (lfn_ni_njk_cutoff_local_max * 1000 - lfn_ni_nik_cutoff * 1000) / nb_points))]

            filter_kwargs_local = filter_kwargs.copy()

            # loop over lfn_nijk_cutoff
            for lfn_nijk_cutoff_item in lfn_nijk_cutoff_list:
                # loop over lfn_ni_nik_cutoff: 0.001, 0.002, ...
                for lfn_ni_nik_cutoff_item in lfn_ni_nik_cutoff_list:

                    ################################################################################
                    #
                    # run filters, then count keep and delete
                    #
                    ################################################################################

                    filter_kwargs_local['lfn_nijk_cutoff'] = lfn_nijk_cutoff_item
                    # filter_kwargs_local['lfn_ni_nik_cutoff'] = lfn_ni_nik_cutoff_item

                    if lfn_nik_cutoff is None:
                        filter_kwargs_local['lfn_ni_cutoff'] = lfn_ni_nik_cutoff_item
                    else:
                        filter_kwargs_local['lfn_nik_cutoff'] = lfn_ni_nik_cutoff_item

                    nijk_remain_df = FilterLFNreplicateRemainRunner(
                        nijk_df=nijk_run_marker_df).get_nijk_remain_df(**filter_kwargs_local)

                    count_keep = nijk_remain_df.merge(
                        keep_run_marker_biosample_variant_df,
                        on=['run_id', 'marker_id', 'biosample_id', 'variant_id']).drop_duplicates(
                        inplace=False).shape[0]

                    count_delete = nijk_remain_df.merge(
                        delete_run_marker_biosample_variant_df,
                        on=['run_id', 'marker_id', 'biosample_id', 'variant_id']).drop_duplicates(
                        inplace=False).shape[0]

                    ################################################################################
                    #
                    # Store results
                    #
                    ################################################################################

                    Logger.instance().debug(
                        "lfn_nijk_cutoff: {} (Max {}); {}: {} (Max {}); "
                        "count_keep_max: {}; count_keep: {}".format(
                            lfn_nijk_cutoff_item, lfn_nijk_cutoff_local_max, which_filter,
                            lfn_ni_nik_cutoff_item,
                            lfn_ni_njk_cutoff_local_max, count_keep_max, count_keep))

                    if count_keep >= count_keep_max:  # Store results if count_keep maximal
                        out_lfn_variant_row_dic = {
                            "lfn_ni_nik_cutoff": lfn_ni_nik_cutoff_item,
                            "lfn_nijk_cutoff": lfn_nijk_cutoff_item,
                            "occurrence_nb_keep": count_keep, "occurrence_nb_delete": count_delete}
                        out_lfn_variant_list.append(out_lfn_variant_row_dic)
                    else:
                        break  # stops when count_keep decreases below count_keep_max

            ########################################################################################
            #
            # Prepare output of this run-marker
            #
            ########################################################################################

            # From list of dics to variant_read_count_input_df
            out_optimize_run_marker_df = pandas.DataFrame(out_lfn_variant_list)
            # List of columns in order
            column_names = ['occurrence_nb_keep', 'occurrence_nb_delete', 'lfn_nijk_cutoff',
                            'lfn_ni_nik_cutoff']
            # Reorder columns
            out_optimize_run_marker_df = out_optimize_run_marker_df[column_names]
            # Sort columns
            out_optimize_run_marker_df.sort_values(by=column_names,
                                                   ascending=[False, True, True, True],
                                                   inplace=True)
            # Rename columns depending on whether this is optimize_lfn_variant or is_optimize_lfn_variant_replicate
            if lfn_nik_cutoff is None:  # optimize lfn variant
                out_optimize_run_marker_df = out_optimize_run_marker_df \
                    .rename(columns={'lfn_ni_nik_cutoff': 'lfn_variant_threshold'})
            else:  # optimize lfn variant replicate
                out_optimize_run_marker_df = out_optimize_run_marker_df \
                    .rename(columns={'lfn_ni_nik_cutoff': 'lfn_variant_replicate_threshold'})

            ########################################################################################
            #
            # Concat
            #
            ########################################################################################

            out_optimize_run_marker_df['run_id'] = run_id
            out_optimize_run_marker_df['marker_id'] = marker_id
            out_optimize_df = pandas.concat([out_optimize_df, out_optimize_run_marker_df], axis=0)

            ########################################################################################
            #
            # Variant delete-specific thresholds
            #
            ########################################################################################

            nijk_run_marker_delete_df = nijk_run_marker_df.merge(
                delete_run_marker_biosample_variant_df, on=['run_id', 'marker_id', 'biosample_id',
                                                            'variant_id'])
            nijk_df_i_obj = VariantReadCountLikeDF(nijk_run_marker_df)

            if lfn_nik_cutoff is None:  # optimize lfn variant
                N_i_df = nijk_df_i_obj.get_N_i_df()

                lfn_ni_or_nik_specific_threshold_df = nijk_run_marker_delete_df.merge(
                    N_i_df, on=['run_id', 'marker_id', 'variant_id'])
                lfn_ni_or_nik_specific_threshold_df[
                    'lfn_variant_threshold'] = lfn_ni_or_nik_specific_threshold_df.read_count \
                                               / lfn_ni_or_nik_specific_threshold_df.N_i
                lfn_ni_or_nik_specific_threshold_df.sort_values(by='lfn_variant_threshold',
                                                                ascending=False, inplace=True)
                lfn_ni_or_nik_specific_threshold_df.drop_duplicates('variant_id', keep='first',
                                                                    inplace=True)
                lfn_ni_or_nik_specific_threshold_df = (
                    lfn_ni_or_nik_specific_threshold_df[
                        ['run_id', 'marker_id', 'variant_id', 'read_count', 'N_i',
                         'lfn_variant_threshold']]).drop_duplicates(
                    inplace=False)
            else:  # optimize lfn variant replicate
                N_ik_df = nijk_df_i_obj.get_N_ik_df()
                lfn_ni_or_nik_specific_threshold_df = nijk_run_marker_delete_df.merge(N_ik_df,
                                                                                      on=['run_id',
                                                                                          'marker_id',
                                                                                          'variant_id',
                                                                                          'replicate'])
                lfn_ni_or_nik_specific_threshold_df[
                    'lfn_variant_replicate_threshold'] = lfn_ni_or_nik_specific_threshold_df.read_count \
                                                         / lfn_ni_or_nik_specific_threshold_df.N_ik
                lfn_ni_or_nik_specific_threshold_df.sort_values(
                    by='lfn_variant_replicate_threshold',
                    ascending=False, inplace=True)
                lfn_ni_or_nik_specific_threshold_df.drop_duplicates(['variant_id', 'replicate'],
                                                                    keep='first', inplace=True)
                lfn_ni_or_nik_specific_threshold_df = (lfn_ni_or_nik_specific_threshold_df[
                    ['run_id', 'marker_id', 'variant_id', 'replicate', 'read_count', 'N_ik',
                     'lfn_variant_replicate_threshold']]).drop_duplicates(inplace=False)

            out_optimize2_df = pandas.concat(
                [out_optimize2_df, lfn_ni_or_nik_specific_threshold_df], axis=0)

        return out_optimize_df, out_optimize2_df

    # def to_tsv(self, optimize_path, engine):
    #
    #     ##########################################################################################
    #     #
    #     # out_optimize_df: Format and write
    #     #
    #     ##########################################################################################
    #
    #     out_optimize_df.marker_id = NameIdConverter(out_optimize_df.marker_id,
    #                                                 engine=engine).to_names(Marker)
    #     out_optimize_df.run_id = NameIdConverter(out_optimize_df.run_id, engine=engine).to_names(
    #         Run)
    #     out_optimize_df.rename({'run_id': 'run', 'marker_id': 'marker'}, axis=1, inplace=True)
    #     out_optimize_df.to_csv(output_file_optimize_lfn_tsv, header=True, sep='\t', index=False)
    #
    #     ##########################################################################################
    #     #
    #     # out_optimize_df: Format and write
    #     #
    #     ##########################################################################################
    #
    #     out_optimize2_df.marker_id = NameIdConverter(out_optimize2_df.marker_id,
    #                                                  engine=engine).to_names(Marker)
    #     out_optimize2_df.run_id = NameIdConverter(out_optimize2_df.run_id, engine=engine).to_names(
    #         Run)
    #     out_optimize2_df['action'] = 'delete'
    #     out_optimize2_df['sequence'] = NameIdConverter(out_optimize2_df.variant_id,
    #                                                    engine=engine).variant_id_to_sequence()
    #     out_optimize2_df.rename({'run_id': 'run', 'marker_id': 'marker', 'variant_id': 'variant',
    #                              'read_count': 'read_count_max'}, axis=1, inplace=True)
    #     out_optimize2_df = out_optimize2_df[
    #         ['run', 'marker', 'variant', 'action', 'read_count_max', 'N_i', 'lfn_variant_threshold',
    #          'sequence']]
    #
    #     out_optimize2_df.to_csv(
    #         output_file_lfn_variant_specific_threshold_tsv, header=True, sep='\t', index=False)

    def run_various_filters(self, is_optimize_lfn_variant_replicate,
                            lfn_biosample_replicate_threshold,
                            lfn_read_count_threshold, min_replicate_number,
                            lfn_variant_or_variant_replicate_threshold):

        lfn_filter_runner = FilterLFNrunner(self.nijk_df)

        ############################################################################################
        #
        # Filter lfn_variant
        #
        ############################################################################################

        if not is_optimize_lfn_variant_replicate:  # optimize lfn variant replicate

            lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_i', threshold=lfn_variant_or_variant_replicate_threshold)

        else:  # optimize lfn variant replicate

            lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_ik', threshold=lfn_variant_or_variant_replicate_threshold)

        ############################################################################################
        #
        # Filter lfn_biosample_replicate
        #
        ############################################################################################

        lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_jk', threshold=lfn_biosample_replicate_threshold)

        ############################################################################################
        #
        # Filter absolute read count
        #
        ############################################################################################

        lfn_filter_runner.mark_delete_lfn_absolute_read_count(lfn_read_count_threshold)

        ############################################################################################
        #
        # mark_delete_lfn_do_not_pass_all_filters
        #
        ############################################################################################

        lfn_filter_runner.mark_delete_lfn_do_not_pass_all_filters()

        variant_read_count_remained_df = lfn_filter_runner.variant_read_count_filter_delete_df

        variant_read_count_remained_df = variant_read_count_remained_df.loc[
            (variant_read_count_remained_df.filter_id == 8) &
            (variant_read_count_remained_df.filter_delete == 0)]
        variant_read_count_remained_df.drop(['filter_id', 'filter_delete'], axis=1, inplace=True)
        del (lfn_filter_runner)

        ############################################################################################
        #
        # FilterMinReplicateNumberRunner
        #
        ############################################################################################

        variant_read_count_remained_df = FilterMinReplicateNumberRunner(
            variant_read_count_remained_df).get_variant_read_count_delete_df(min_replicate_number)

        variant_read_count_remained_df = variant_read_count_remained_df.loc[
            variant_read_count_remained_df.filter_delete == 0]
        variant_read_count_remained_df.drop('filter_delete', axis=1, inplace=True)

        ############################################################################################
        #
        # Count keep
        #
        ############################################################################################

        # variant_read_count_remained_df.drop_duplicates(inplace=True)

        # variant_read_count_remained_keep_df = variant_read_count_remained_df.merge(variant_keep_df,
        #                                                                            on=['run_id', 'marker_id',
        #                                                                                'biosample_id', 'variant_id'])
        # variant_read_count_remained_keep_df = variant_read_count_remained_keep_df[
        #     ['run_id', 'marker_id', 'variant_id', 'biosample_id']].drop_duplicates()
        # count_keep = variant_read_count_remained_keep_df.shape[0]

        return variant_read_count_remained_df[['run_id', 'marker_id', 'biosample_id', 'variant_id']].drop_duplicates(inplace=False)
