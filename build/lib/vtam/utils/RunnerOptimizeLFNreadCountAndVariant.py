import pandas
from vtam.utils.constants import lfn_ni_njk_cutoff_global_max, lfn_ni_njk_cutoff_lst_size, \
    lfn_nijk_cutoff_global_max, \
    lfn_nijk_cutoff_lst_size

from vtam.utils.RunnerOptimizeLFNreadCountAndVariantRunMarker import \
    RunnerOptimizeLFNreadCountAndVariantRunMarker


class RunnerOptimizeLFNreadCountAndVariant:

    """This the Runner for Optimize LFN readcount and variant/variantReplicate
    in the presence of several run-marker combinations"""

    def __init__(self, nijk_df, known_occurrences_df):

        self.nijk_df = nijk_df
        self.known_occurrences_df = known_occurrences_df

    def get_optimize_df(self, lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff,
                        min_replicate_number):

        ############################################################################################
        #
        # Group and run_name this genetic_code by run_name/marker_name combination
        #  Loop by run_name/marker_name
        #
        ############################################################################################

        out_optimize_df = pandas.DataFrame()
        out_optimize2_df = pandas.DataFrame()

        lfn_nijk_cutoff_lst = RunnerOptimizeLFNreadCountAndVariantRunMarker.get_lfn_nijk_cutoff_lst(lfn_nijk_cutoff, lfn_nijk_cutoff_global_max, lfn_nijk_cutoff_lst_size)
        if lfn_nik_cutoff is None:
            lfn_ni_nik_cutoff_lst = RunnerOptimizeLFNreadCountAndVariantRunMarker.get_lfn_ni_nik_cutoff_lst(lfn_ni_cutoff, lfn_ni_njk_cutoff_global_max, lfn_ni_njk_cutoff_lst_size)
        else:
            lfn_ni_nik_cutoff_lst = RunnerOptimizeLFNreadCountAndVariantRunMarker.get_lfn_ni_nik_cutoff_lst(
                lfn_nik_cutoff, lfn_ni_njk_cutoff_global_max, lfn_ni_njk_cutoff_lst_size)

        for row in self.known_occurrences_df[['run_id', 'marker_id']].drop_duplicates().itertuples():

            run_id = row.run_id
            marker_id = row.marker_id
            known_occurrs_run_marker_df = self.known_occurrences_df.loc[
                (self.known_occurrences_df.run_id == run_id) & (
                            self.known_occurrences_df.marker_id == marker_id),]
            nijk_run_marker_df = self.nijk_df.loc[
                (self.nijk_df.run_id == run_id) & (self.nijk_df.marker_id == marker_id),]
            nijk_run_marker_df = nijk_run_marker_df[
                ['run_id', 'marker_id', 'sample_id', 'replicate', 'variant_id',
                 'read_count']].drop_duplicates(inplace=False)

            optim_run_marker_obj = RunnerOptimizeLFNreadCountAndVariantRunMarker(
                nijk_run_marker_df, known_occurrs_run_marker_df, lfn_nijk_cutoff_lst, lfn_ni_nik_cutoff_lst)
            out_optimize_run_marker_df = optim_run_marker_obj.get_df_optim_lfn_readcount_variant_replicate_cutoff(
                lfn_ni_cutoff=lfn_ni_cutoff, lfn_nik_cutoff=lfn_nik_cutoff,
                lfn_njk_cutoff=lfn_njk_cutoff, lfn_nijk_cutoff=lfn_nijk_cutoff,
                min_replicate_number=min_replicate_number)

            ########################################################################################
            #
            # Prepare output of this run-marker
            #
            ########################################################################################

            # From list of dics to variant_read_count_input_df
            # out_optimize_run_marker_df = pandas.DataFrame(out_lfn_variant_list)
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
            out_optimize_run_marker_df.lfn_ni_nik_cutoff = round(out_optimize_run_marker_df.lfn_ni_nik_cutoff, 3)
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

            lfn_ni_or_nik_specific_cutoff_df = optim_run_marker_obj.get_df_variant_specific_cutoffs(lfn_nik_cutoff)

            out_optimize2_df = pandas.concat(
                [out_optimize2_df, lfn_ni_or_nik_specific_cutoff_df], axis=0)

        return out_optimize_df, out_optimize2_df
