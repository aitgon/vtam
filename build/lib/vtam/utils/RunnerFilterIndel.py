class RunnerFilterIndel:

    def __init__(self, variant_read_count_df):

        self.variant_read_count_df = variant_read_count_df

    def get_variant_read_count_delete_df(self, variant_df, skip_filter_indel):
        """
        filter chimera
        """

        variant_read_count_delete_df = self.variant_read_count_df.copy()

        if skip_filter_indel:

            variant_read_count_delete_df['filter_delete'] = False

        else:

            variant_read_count_delete_df = self.variant_read_count_df.copy()
            variant_read_count_delete_df['filter_delete'] = False
            #
            df = variant_df.copy()
            df['sequence_length_module_3'] = variant_df.sequence.apply(
                lambda x: len(x) % 3)  # compute module for each variant
            #  most common remaining of modulo 3
            majority_sequence_length_module_3 = df.sequence_length_module_3.mode()
            # select id of variant that do not pass on a list
            df = df.loc[df['sequence_length_module_3'] !=
                        majority_sequence_length_module_3.values[0]]
            #
            for id in df.index.tolist():
                variant_read_count_delete_df.loc[variant_read_count_delete_df['variant_id']
                                                 == id, 'filter_delete'] = True

        return variant_read_count_delete_df
