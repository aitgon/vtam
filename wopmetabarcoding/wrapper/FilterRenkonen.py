









def f11_renkonen(self, number_of_replicate, renkonen_tail):
    """
    Function calculating distance between sample replicate of the sample to determine if the sample replicate
    must be deleted or not
    :param engine: engine of the database
    :param replicate_model: model of the replicate table
    :param variant_biosample_replicate_df: data frame containing the information
    :param number_of_replicate: Number of replicate by sample
    :return: None
    """
    # TODO: Must be updated for new FilterNonLFNRunner class
    return
    this_filter_name = inspect.stack()[0][3]
    logger.debug(
        "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
    ########################################################
    # proportion of the reads of variant i per replicate j (Ni,j=1/Nj=1)
    ########################################################
    variant_read_proportion_per_replicate_df = self.variant_read_count_df[
        ['biosample_id', 'replicate_id', 'read_count']].groupby(
        by=['biosample_id', 'replicate_id']).sum().reset_index()
    # Merge the column with the total reads by sample replicates for calculate the ratio
    variant_read_proportion_per_replicate_df = self.variant_read_count_df.merge(
        variant_read_proportion_per_replicate_df, left_on=['biosample_id', 'replicate_id'],
        right_on=['biosample_id', 'replicate_id'])
    variant_read_proportion_per_replicate_df.columns = ['variant_id', 'biosample_id', 'replicate_id',
                                                        'read_count_per_variant_per_biosample_per_replicate',
                                                        'read_count_per_biosample_replicate']
    variant_read_proportion_per_replicate_df[
        'read_proportion_of_variant_in_replicate'] = variant_read_proportion_per_replicate_df.read_count_per_variant_per_biosample_per_replicate / variant_read_proportion_per_replicate_df.read_count_per_biosample_replicate
    # for biosample_id in self.variant_biosample_replicate_df.biosample_id.unique():
    # replicate_combinatorics = itertools.permutations(self.variant_biosample_replicate_df.replicate_id.unique().tolist(), 2)
    biosample_id = 1
    replicate_id1 = 1
    replicate_id2 = 2
    variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
        variant_read_proportion_per_replicate_df.biosample_id == biosample_id]
    ########################################################
    # 2. Calculate renkonen distance index (D) for all pairs of replicates of the same sample
    ########################################################
    variant_read_proportion_per_replicate1_per_biosample_df = \
    variant_read_proportion_per_replicate_per_biosample_df.loc[
        variant_read_proportion_per_replicate_per_biosample_df.replicate_id == replicate_id1, ['variant_id',
                                                                                               'replicate_id',
                                                                                               'read_proportion_of_variant_in_replicate']]
    variant_read_proportion_per_replicate2_per_biosample_df = \
    variant_read_proportion_per_replicate_per_biosample_df.loc[
        variant_read_proportion_per_replicate_per_biosample_df.replicate_id == replicate_id2, ['variant_id',
                                                                                               'replicate_id',
                                                                                               'read_proportion_of_variant_in_replicate']]
    variant_read_proportion_per_replicate_1_2 = variant_read_proportion_per_replicate1_per_biosample_df.merge(
        variant_read_proportion_per_replicate2_per_biosample_df, on='variant_id')
    variant_read_proportion_per_replicate_1_2.columns = ['variant_id', 'replicate_id1',
                                                         'read_proportion_of_variant_in_replicate1',
                                                         'replicate_id2',
                                                         'read_proportion_of_variant_in_replicate_2']
    variant_read_proportion_per_replicate_1_2 = variant_read_proportion_per_replicate_1_2[
        ['variant_id', 'read_proportion_of_variant_in_replicate1', 'read_proportion_of_variant_in_replicate2']]
    variant_read_proportion_per_replicate_1_2['min_read_proportion'] = variant_read_proportion_per_replicate_1_2[
        ['read_proportion_of_variant_in_replicate1', 'read_proportion_of_variant_in_replicate2']].apply(
        lambda row: row.min(), axis=1)
    #
    columns_name = ['repl_i', 'repl_j', 'distance']
    df_read_count_per_sample_replicate = self.variant_read_count_df.groupby(by=['sample_replicate'])['count'].sum()
    df_read_count_per_sample_replicate = df_read_count_per_sample_replicate.to_frame()
    df_read_count_per_sample_replicate.columns = ['replicate_count']
    df_read_count_per_sample_replicate = self.variant_read_count_df.merge(df_read_count_per_sample_replicate,
                                                                          left_on='sample_replicate', right_index=True)
    df_read_count_per_sample_replicate['proportion'] = df_read_count_per_sample_replicate['count'] / \
                                                       df_read_count_per_sample_replicate['replicate_count']
    # df_replicate = df_read_count_per_sample_replicate.groupby(by=['biosample'])['sample_replicate'].to_frame()
    samples = df_read_count_per_sample_replicate['biosample']
    samples = list(set(samples.tolist()))
    for sample in samples:
        df_permutation_distance = pandas.DataFrame(columns=columns_name)
        df_replicate = df_read_count_per_sample_replicate.loc[df_read_count_per_sample_replicate['biosample'] == sample]
        replicates = list(set(df_replicate['sample_replicate'].tolist()))
        for combi in itertools.permutations(replicates, 2):
            combi = list(combi)
            df_repli = df_replicate.loc[df_replicate['sample_replicate'] == combi[0]]
            data_repli = df_repli[['variant_seq', 'sample_replicate', 'proportion']]
            df_replj = df_replicate.loc[df_replicate['sample_replicate'] == combi[1]]
            data_replj = df_replj[['variant_seq', 'sample_replicate', 'proportion']]
            df_replij = data_repli.append(data_replj)
            group_repl = df_replij.groupby(by=['variant_seq'])['proportion'].min()
            distance = 1 - group_repl.sum()
            query = [combi[0], combi[1], distance]
            df_permutation_distance.loc[len(df_permutation_distance)] = query
        # df_calc = df_permutation_distance.loc[df_permutation_distance['repl_i'] == combi[0]]
        indices_to_drop = list(
            df_permutation_distance.loc[df_permutation_distance.distance > renkonen_tail].index.tolist()
        )
        df_permutation_distance.drop(indices_to_drop, inplace=True)
        repl_list = list(set(df_permutation_distance['repl_i'].tolist()))
        for repl_i in repl_list:
            df_calc = df_permutation_distance.loc[df_permutation_distance['repl_i'] == repl_i]
            if len(df_calc) > ((number_of_replicate - 1) / 2):
                index_to_drop = self.variant_read_count_df.loc[
                    self.variant_read_count_df['sample_replicate'] == repl_i].index.tolist()
                self.passsed_variant_ids = sorted(list(set(index_to_drop + self.passsed_variant_ids)))