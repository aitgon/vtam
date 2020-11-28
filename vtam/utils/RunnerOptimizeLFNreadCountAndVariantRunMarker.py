import pandas
import numpy
from vtam.utils.VTAMexception import VTAMexception

from vtam.utils.RunnerFilterLFNreplicateRemain import RunnerFilterLFNreplicateRemain
from vtam.utils.DataframeVariantReadCountLike import DataframeVariantReadCountLike


class RunnerOptimizeLFNreadCountAndVariantRunMarker:

    """This the Runner for Optimize LFN readcount and variant/variantReplicate
    in the presence of one run-marker combination"""

    def __init__(self, nijk_df, known_occurrences_df, lfn_nijk_cutoff_lst, lfn_ni_nik_cutoff_lst):

        if not (len(nijk_df.marker_id.unique()) == 1 or len(nijk_df.run_id.unique()) == 1):
            raise (VTAMexception(
                "Internal error. This variant_read_count df requires a single run_id and marker_id"))

        if not (len(known_occurrences_df.marker_id.unique()) == 1 or len(known_occurrences_df.run_id.unique()) == 1):
            raise (VTAMexception(
                "Internal error. This variant_read_count df requires a single run_id and marker_id"))

        self.nijk_df = nijk_df

        self.known_occurrences_df = known_occurrences_df

        self.lfn_nijk_cutoff_lst = lfn_nijk_cutoff_lst
        self.lfn_ni_nik_cutoff_lst = lfn_ni_nik_cutoff_lst

    @classmethod
    def get_lfn_nijk_cutoff_lst(cls, start: object, stop: object, nb_points: object) -> object:

        return [*range(start, stop, round(int((stop - start + 1) / nb_points), -1))]

    @classmethod
    def get_lfn_ni_nik_cutoff_lst(cls, start, stop, nb_points):

        return [round(x, 3) for x in numpy.arange(start, stop, (stop - start + 0.001)/nb_points)]

    def get_count_keep_max(self):

        count_delete_max = len(self.known_occurrences_df.loc[self.known_occurrences_df.action == 'keep'].variant_id.unique())
        return count_delete_max

    def get_count_delete_max(self):

        count_delete_max = len(self.known_occurrences_df.loc[self.known_occurrences_df.action == 'delete'].variant_id.unique())
        return count_delete_max

    def get_lst_one_par_lfn_nijk_cutoff(self, lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff, min_replicate_number):

        """Loops through self.lfn_nijk_cutoff_lst and keeps only values while count_keep < count_keep_max"""

        out_lfn_nijk_cutoff_lst = []

        count_keep_max = self.get_count_keep_max()

        for lfn_nijk_cutoff_item in self.lfn_nijk_cutoff_lst:

            count_keep, count_delete = RunnerFilterLFNreplicateRemain(
                nijk_df=self.nijk_df, lfn_ni_cutoff=lfn_ni_cutoff, lfn_nik_cutoff=lfn_nik_cutoff, lfn_njk_cutoff=lfn_njk_cutoff, lfn_nijk_cutoff=lfn_nijk_cutoff_item,
                           min_replicate_number=min_replicate_number)\
                .count_keep_delete(known_occurrences_df=self.known_occurrences_df)

            if count_keep < count_keep_max:
                break  # stops when count_keep decreases below count_keep_max

            out_lfn_nijk_cutoff_lst.append(lfn_nijk_cutoff_item)

        return out_lfn_nijk_cutoff_lst

    def get_lst_one_par_lfn_ni_nik_cutoff(self, lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff, min_replicate_number):

        """Loops between default lfn_nijk_cutoff and vtam.utils.constants.lfn_nijk_cutoff_global_max to obtain count_keep_max"""

        out_lfn_ni_nik_cutoff_lst = []

        # lfn_ni_nik_cutoff_lst = self.get_lfn_ni_nik_cutoff_lst(start=lfn_ni_nik_cutoff, stop=lfn_ni_njk_cutoff_global_max, nb_points=lfn_ni_njk_cutoff_lst_size)
        count_keep_max = self.get_count_keep_max()
        for lfn_ni_nik_cutoff_item in self.lfn_ni_nik_cutoff_lst:

            if lfn_nik_cutoff is None:
                count_keep, count_delete = RunnerFilterLFNreplicateRemain(
                    nijk_df=self.nijk_df, lfn_ni_cutoff=lfn_ni_nik_cutoff_item, lfn_nik_cutoff=lfn_nik_cutoff, lfn_njk_cutoff=lfn_njk_cutoff, lfn_nijk_cutoff=lfn_nijk_cutoff,
                               min_replicate_number=min_replicate_number)\
                    .count_keep_delete(known_occurrences_df=self.known_occurrences_df)
            else:
                count_keep, count_delete = RunnerFilterLFNreplicateRemain(
                    nijk_df=self.nijk_df, lfn_ni_cutoff=lfn_ni_cutoff, lfn_nik_cutoff=lfn_ni_nik_cutoff_item, lfn_njk_cutoff=lfn_njk_cutoff, lfn_nijk_cutoff=lfn_nijk_cutoff,
                               min_replicate_number=min_replicate_number)\
                    .count_keep_delete(known_occurrences_df=self.known_occurrences_df)


            if count_keep < count_keep_max:
                break  # stops when count_keep decreases below count_keep_max

            out_lfn_ni_nik_cutoff_lst.append(lfn_ni_nik_cutoff_item)

        return out_lfn_ni_nik_cutoff_lst

    def get_df_optim_lfn_readcount_variant_replicate_cutoff(self, lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff, min_replicate_number):

        """Two parameter loop for lfn_nijk_cutoff and lfn_ni_cutoff/lfn_nik_cutoff to get keep_nb, delete_nb with the two parameters"""

        out_two_pars_lst = []
        # lfn_ni_nik_cutoff = lfn_ni_cutoff  # current max
        # if not (lfn_nik_cutoff is None):
        #     lfn_ni_nik_cutoff = lfn_nik_cutoff  # current max

        lfn_nijk_cutoff_lst = self.get_lst_one_par_lfn_nijk_cutoff(lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff, min_replicate_number)
        lfn_ni_nik_cutoff_lst = self.get_lst_one_par_lfn_ni_nik_cutoff(lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff, min_replicate_number)

        count_keep_max = self.get_count_keep_max()
        # loop over lfn_nijk_cutoff
        for lfn_nijk_cutoff_item in lfn_nijk_cutoff_lst:
            # loop over lfn_ni_nik_cutoff: 0.001, 0.002, ...
            for lfn_ni_nik_cutoff_item in lfn_ni_nik_cutoff_lst:

                if lfn_nik_cutoff is None:
                    count_keep, count_delete = RunnerFilterLFNreplicateRemain(
                        nijk_df=self.nijk_df, lfn_ni_cutoff=lfn_ni_nik_cutoff_item,
                        lfn_nik_cutoff=lfn_nik_cutoff, lfn_njk_cutoff=lfn_njk_cutoff,
                        lfn_nijk_cutoff=lfn_nijk_cutoff_item,
                        min_replicate_number=min_replicate_number) \
                        .count_keep_delete(known_occurrences_df=self.known_occurrences_df)
                else:
                    count_keep, count_delete = RunnerFilterLFNreplicateRemain(
                        nijk_df=self.nijk_df, lfn_ni_cutoff=lfn_ni_cutoff,
                        lfn_nik_cutoff=lfn_ni_nik_cutoff_item, lfn_njk_cutoff=lfn_njk_cutoff,
                        lfn_nijk_cutoff=lfn_nijk_cutoff_item,
                        min_replicate_number=min_replicate_number) \
                        .count_keep_delete(known_occurrences_df=self.known_occurrences_df)

                ################################################################################
                #
                # Store results
                #
                ################################################################################

                if count_keep >= count_keep_max:  # Store results if count_keep maximal
                    out_lfn_variant_row_dic = {
                        "lfn_ni_nik_cutoff": lfn_ni_nik_cutoff_item,
                        "lfn_nijk_cutoff": lfn_nijk_cutoff_item,
                        "occurrence_nb_keep": count_keep, "occurrence_nb_delete": count_delete}
                    out_two_pars_lst.append(out_lfn_variant_row_dic)
                else:
                    break  # stops when count_keep decreases below count_keep_max

        out_two_pars_df = pandas.DataFrame(out_two_pars_lst)
        column_names = ['occurrence_nb_keep', 'occurrence_nb_delete', 'lfn_nijk_cutoff',
                        'lfn_ni_nik_cutoff']
        out_two_pars_df = out_two_pars_df[column_names]
        out_two_pars_df.sort_values(by=column_names,
                                               ascending=[False, True, False, False],
                                               inplace=True)
        return out_two_pars_df

    def get_df_variant_specific_cutoffs(self, lfn_nik_cutoff):

        """Two parameter loop for lfn_nijk_cutoff and lfn_ni_cutoff/lfn_nik_cutoff to get keep_nb, delete_nb with the two parameters"""

        delete_run_marker_sample_variant_df = self.known_occurrences_df.loc[
            self.known_occurrences_df.action == 'delete', ['run_id', 'marker_id', 'sample_id',
                                                           'variant_id']]

        nijk_run_marker_delete_df = self.nijk_df.merge(
            delete_run_marker_sample_variant_df, on=['run_id', 'marker_id', 'sample_id',
                                                        'variant_id'])
        nijk_df_i_obj = DataframeVariantReadCountLike(self.nijk_df)

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

        return lfn_ni_or_nik_specific_cutoff_df