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
                "Searching upper limit of lfn_variant_cutoff or lfn_variant_replicate_cutoff")

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
                    .rename(columns={'lfn_ni_nik_cutoff': 'lfn_variant_cutoff'})
            else:  # optimize lfn variant replicate
                out_optimize_run_marker_df = out_optimize_run_marker_df \
                    .rename(columns={'lfn_ni_nik_cutoff': 'lfn_variant_replicate_cutoff'})

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
            # Variant delete-specific cutoffs
            #
            ########################################################################################

            nijk_run_marker_delete_df = nijk_run_marker_df.merge(
                delete_run_marker_biosample_variant_df, on=['run_id', 'marker_id', 'biosample_id',
                                                            'variant_id'])
            nijk_df_i_obj = VariantReadCountLikeDF(nijk_run_marker_df)

            if lfn_nik_cutoff is None:  # optimize lfn variant
                N_i_df = nijk_df_i_obj.get_N_i_df()

                lfn_ni_or_nik_specific_cutoff_df = nijk_run_marker_delete_df.merge(
                    N_i_df, on=['run_id', 'marker_id', 'variant_id'])
                lfn_ni_or_nik_specific_cutoff_df[
                    'lfn_variant_cutoff'] = lfn_ni_or_nik_specific_cutoff_df.read_count \
                                               / lfn_ni_or_nik_specific_cutoff_df.N_i
                lfn_ni_or_nik_specific_cutoff_df.sort_values(by='lfn_variant_cutoff',
                                                                ascending=False, inplace=True)
                lfn_ni_or_nik_specific_cutoff_df.drop_duplicates('variant_id', keep='first',
                                                                    inplace=True)
                lfn_ni_or_nik_specific_cutoff_df = (
                    lfn_ni_or_nik_specific_cutoff_df[
                        ['run_id', 'marker_id', 'variant_id', 'read_count', 'N_i',
                         'lfn_variant_cutoff']]).drop_duplicates(
                    inplace=False)
            else:  # optimize lfn variant replicate
                N_ik_df = nijk_df_i_obj.get_N_ik_df()
                lfn_ni_or_nik_specific_cutoff_df = nijk_run_marker_delete_df.merge(N_ik_df,
                                                                                      on=['run_id',
                                                                                          'marker_id',
                                                                                          'variant_id',
                                                                                          'replicate'])
                lfn_ni_or_nik_specific_cutoff_df[
                    'lfn_variant_replicate_cutoff'] = lfn_ni_or_nik_specific_cutoff_df.read_count \
                                                         / lfn_ni_or_nik_specific_cutoff_df.N_ik
                lfn_ni_or_nik_specific_cutoff_df.sort_values(
                    by='lfn_variant_replicate_cutoff',
                    ascending=False, inplace=True)
                lfn_ni_or_nik_specific_cutoff_df.drop_duplicates(['variant_id', 'replicate'],
                                                                    keep='first', inplace=True)
                lfn_ni_or_nik_specific_cutoff_df = (lfn_ni_or_nik_specific_cutoff_df[
                    ['run_id', 'marker_id', 'variant_id', 'replicate', 'read_count', 'N_ik',
                     'lfn_variant_replicate_cutoff']]).drop_duplicates(inplace=False)

            out_optimize2_df = pandas.concat(
                [out_optimize2_df, lfn_ni_or_nik_specific_cutoff_df], axis=0)

        return out_optimize_df, out_optimize2_df
