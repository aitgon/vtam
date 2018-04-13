from sqlalchemy import select, delete
import pandas, subprocess


def lfn1_per_replicate(df, lfn_per_replicate_threshold):
    df2 = df.groupby(by=['sample_replicate']).sum()
    df2 = df.merge(df2, left_on='sample_replicate', right_index=True)
    df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_sample_replicate']
    df2['low_frequence_noice_per_replicate'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_sample_replicate
    return list(df2.ix[df2.low_frequence_noice_per_replicate < lfn_per_replicate_threshold].index) # select rare reads to discard later

def lfn2_per_variant(df, lfn_per_variant_threshold):
    df2 = df.groupby(by=['sample_replicate']).sum()
    df2 = df.merge(df2, left_on='sample_replicate', right_index=True)
    df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_variant']
    df2['low_frequence_noice_per_variant'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_variant
    return list(df2.ix[df2.low_frequence_noice_per_variant < lfn_per_variant_threshold].index) # select rare reads to discard later

def lfn2_per_replicate_series(df, lfn_per_replicate_series_threshold):
    df2 = df.groupby(by=['replicate']).sum()
    df2 = df.merge(df2, left_on='replicate', right_index=True)
    df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_replicate_series']
    df2['low_frequence_noice_per_replicate_series'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_replicate_series
    return list(df2.ix[df2.low_frequence_noice_per_replicate_series < lfn_per_replicate_series_threshold].index) # select rare reads to discard later

def lfn3_read_count(df, lfn_read_count_threshold):
    return list(df.ix[df['count'] <= lfn_read_count_threshold].index) # select rare reads to discard later

def lfn4_per_variant_with_cutoff(df, cutoff_tsv):
    cutoff_df = pandas.read_csv(cutoff_tsv, sep="\t", header=0)
    cutoff_df.columns = ['sequence', 'cutoff']
    #
    df2 = df.groupby(by=['sample_replicate']).sum()
    df2 = df.merge(df2, left_on='sample_replicate', right_index=True)
    df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_variant']
    df2['low_frequence_noice_per_variant'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_variant
    #
    # merge with cutoff
    df2 = df.merge(cutoff_df, left_on="sequence", right_on="sequence")
    return list(df2.ix[df2.low_frequence_noice_per_variant < df2.cutoff].index) # select rare reads to discard later

# def lfn_per_replicate(engine, replicate_model, variant_model, marker_id, data_frame):
#     """
#     Filter out Low Frequency Noise which eliminate a variant from sample-replicate if Nvar-repl/Nrepl < lfn_repl
#     :param engine: Engine of the database
#     :param replicate_model: Model of the replicate table
#     :param variant_model: Model of the variant table
#     :param marker_id: Name of the marker_id
#     :param data_frame: Dataframe containing the content of the sample_count tsv
#     :return: list failed_variants
#     """
#     # Selecting lines of the replicate table
#     replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
#         replicate_model.marker_id == marker_id)
#     replicate_obj = engine.execute(replicate_select)
#     failed_variants = {}
#     for replicate in replicate_obj:
#         failed_variants_list = []
#         sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
#         # Selecting lines of the variant table
#         variant_select = select([variant_model.variant_id, variant_model.sequence]) \
#             .where(variant_model.marker == marker_id)
#         variant_obj = engine.execute(variant_select)
#         # Creating a data frame with only the columns containing the sample replicate
#         data_replicate = data_frame.loc[data_frame['sample_replicate'] == sample_replicate]
#         # Checking if the data frame is not empty
#         if data_replicate.empty is False:
#             replicate_count_series = data_replicate['count']
#             replicate_count = replicate_count_series.sum()
#         for variant in variant_obj:
#             data_variant = data_replicate.loc[data_replicate['sequence'] == variant.sequence]
#             if data_variant.empty is False and data_replicate.empty is False:
#                 variant_replicate_count_series = data_variant['count']
#                 variant_replicate_count = variant_replicate_count_series.sum()
#                 lfn_value = variant_replicate_count / replicate_count
#                 if lfn_value < 0.001 and variant.variant_id not in failed_variants_list:
#                     failed_variants_list.append(variant.variant_id)
#         failed_variants[sample_replicate] = failed_variants_list
#     return failed_variants


# def lfn_per_variant_tmp(engine, replicate_model, variant_model, marker_id, data_frame, series_replicate):
#         # Selecting lines of the variant table
#     variant_select = select([variant_model.variant_id, variant_model.sequence]) \
#             .where(variant_model.marker == marker_id)
#     variant_obj = engine.execute(variant_select)
#     failed_variants = {}
#     for variant in variant_obj:
#         # Selecting lines of the replicate table
#         replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
#                 replicate_model.marker_id == marker_id)
#         replicate_obj = engine.execute(replicate_select)
#         failed_variants_list = []
#         data_variant = data_frame.loc[data_frame['sequence'] == variant.sequence]
#         if data_variant.empty is False:
#             variant_count_series = data_variant['count']
#             variant_count = variant_count_series.sum()
#         for replicate in replicate_obj:
#             sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
#             if series_replicate is True:
#                 data_replicate = data_frame.loc[data_frame['replicate'] == replicate.name]
#                 if data_replicate.empty is False:
#                     replicate_series = data_replicate['count']
#                     replicate_count = replicate_series.sum()
#             data_variantreplicate = data_variant.loc[data_variant['sample_replicate'] == sample_replicate]
#             if data_variantreplicate.empty is False and data_variant.empty is False:
#                 variantreplicate_count_series = data_variantreplicate['count']
#                 variantreplicate_count = variantreplicate_count_series.sum()
#                 variant_count = variant_count_series.sum()
#                 if series_replicate is True:
#                     lfn_value = variantreplicate_count / replicate_count
#                 else:
#                     lfn_value = variantreplicate_count / variant_count
#                 if lfn_value < 0.001 and variant.variant_id not in failed_variants_list:
#                     failed_variants_list.append(variant.sequence)
#             failed_variants[sample_replicate] = failed_variants_list
#     return failed_variants
#
#
# def lfn_per_readcounts(engine, replicate_model, variant_model, marker_id, custom_filter, data_frame):
#     """
#
#     :param engine:
#     :param replicate_model:
#     :param variant_model:
#     :param marker_id:
#     :param custom_filter:
#     :param data_frame:
#     :return:
#     """
#     replicate_select = select(
#         [replicate_model.biosample_name, replicate_model.name]).where(replicate_model.marker_id == marker_id)
#     replicate_obj = engine.execute(replicate_select)
#     failed_variants = {}
#     for replicate in replicate_obj:
#         sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
#         # print(sample_replicate)
#         failed_variants_list = []
#         variant_select = select([variant_model.variant_id, variant_model.sequence]) \
#             .where(variant_model.marker == marker_id)
#         variant_obj = engine.execute(variant_select)
#         data_replicate = data_frame.loc[data_frame['sample_replicate'] == sample_replicate]
#         for variant in variant_obj:
#             data_variant = data_replicate.loc[data_replicate['sequence'] == variant.sequence]
#             if data_variant.empty is False and data_replicate.empty is False:
#                 variant_replicate_count_series = data_variant['count']
#                 variant_replicate_count = variant_replicate_count_series.sum()
#                 if sample_replicate == "14Mon06_repl1":
#                     print(variant.sequence)
#                     print(variant_replicate_count)
#                 if variant_replicate_count < custom_filter and variant.variant_id not in failed_variants_list:
#                     if sample_replicate == "14Mon06_repl1":
#                         print("Denied: " + variant.sequence)
#                     failed_variants_list.append(variant.sequence)
#         failed_variants[sample_replicate] = failed_variants_list
#         # print(failed_variants_list)
#     return failed_variants


def create_cutoff_table(cutoff_file_tsv):
    data_cutoff = pandas.read_csv(cutoff_file_tsv, sep='\t')
    return data_cutoff


def lfn_per_cutoff(engine, replicate_model, variant_model, marker_id, data_frame, cutoff_tsv, series_replicate):
    data_cutoff = create_cutoff_table(cutoff_tsv)
        # Selecting lines of the variant table
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
            .where(variant_model.marker == marker_id)
    variant_obj = engine.execute(variant_select)
    failed_variants = {}
    for variant in variant_obj:
        # Selecting lines of the replicate table
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
                replicate_model.marker_id == marker_id)
        replicate_obj = engine.execute(replicate_select)
        failed_variants_list = []
        data_variant = data_frame.loc[data_frame['sequence'] == variant.sequence]
        if data_variant.empty is False:
            variant_count_series = data_variant['count']
            variant_count = variant_count_series.sum()
        for replicate in replicate_obj:
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            if series_replicate is True:
                data_replicate = data_frame.loc[data_frame['replicate'] == replicate.name]
                if data_replicate.empty is False:
                    replicate_series = data_replicate['count']
                    replicate_count = replicate_series.sum()
            data_variantreplicate = data_variant.loc[data_variant['sample_replicate'] == sample_replicate]
            if data_variantreplicate.empty is False and data_variant.empty is False:
                variantreplicate_count_series = data_variantreplicate['count']
                variantreplicate_count = variantreplicate_count_series.sum()
                variant_count = variant_count_series.sum()
                data_cutoff_filter = data_cutoff.loc[data_cutoff['sequence'] == variant.sequence]
                cutoff_series = data_cutoff_filter['value']
                cutoff = cutoff_series.sum()
                if series_replicate is True:
                    lfn_value = variantreplicate_count / replicate_count
                else:
                    lfn_value = variantreplicate_count / variant_count
                if lfn_value < cutoff and variant.variant_id not in failed_variants_list:
                    failed_variants_list.append(variant.sequence)
            failed_variants[sample_replicate] = failed_variants_list
    return failed_variants


# def min_repln(engine, sample_count_tsv, variant_list, variant_model, data_frame ,value):
#     succeed_variants = []
#     for variant in variant_list:
#         variant_select = select([variant_model.variant_id, variant_model.sequence]).where(variant_model.variant_id == variant)
#         variant_obj = engine.execute(variant_select)
#         for table_variants in variant_obj:
#             print(table_variants)
#             with open(sample_count_tsv, 'r') as fin:
#                 count = 0
#                 for line in fin:
#                     if table_variants.sequence in line:
#                         count += 1
#                 if count >= value:
#                     print(table_variants.variant_id)
#                     print(count)
#                     succeed_variants.append(variant)
#     return succeed_variants


def min_repln(engine, variant_model, replicate_model, marker_id, data_frame, min_count):
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_id)
    variant_obj = engine.execute(variant_select)
    for variant in variant_obj:
        data_variant = data_frame.loc[data_frame['sequence'] == variant.sequence]
        # Selecting lines of the replicate table
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_id == marker_id)
        replicate_obj = engine.execute(replicate_select)
        for replicate in replicate_obj:
            data_sample = data_variant.loc[data_variant['sample'] == replicate.biosample_name]
            if data_variant.empty is False and data_sample.empty is False:
                repl_count = len(data_sample.index)
                if repl_count <= min_count:
                    data_to_drop = data_frame.loc[
                        (data_frame['sequence'] == variant.sequence) &
                        (data_frame['sample'] == replicate.biosample_name)
                        ]
                    indexs_to_drop = data_to_drop.index.tolist()
                    data_frame.drop(indexs_to_drop, inplace=True)


def min_replp(engine, variant_model, replicate_model, marker_id, data_frame, replicate_number):
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_id)
    variant_obj = engine.execute(variant_select)
    for variant in variant_obj:
        data_variant = data_frame.loc[data_frame['sequence'] == variant.sequence]
        # Selecting lines of the replicate table
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_id == marker_id)
        replicate_obj = engine.execute(replicate_select)
        for replicate in replicate_obj:
            data_sample = data_variant.loc[data_variant['sample'] == replicate.biosample_name]
            repl_count = len(data_sample.index)
            if repl_count <= ((1/3)*replicate_number):
                data_to_drop = data_frame.loc[
                    (data_frame['sequence'] == variant.sequence) &
                    (data_frame['sample'] == replicate.biosample_name)
                    ]
                indexs_to_drop = data_to_drop.index.tolist()
                data_frame.drop(indexs_to_drop, inplace=True)


def delete_filtered_variants(engine, replicate_model, marker_id, data_frame, filter1, filter2, filter3, filter4):
    """
    Function which delete variants which don't pass lfn filters
    :param engine:
    :param variant_model:
    :param failed_variants:
    :return:
    """
    replicate_select = select([replicate_model.name, replicate_model.biosample_name])\
        .where(replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(replicate_select)
    for replicate in replicate_obj:
        sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
        filter1_list = filter1.get(sample_replicate)
        filter2_list = filter2.get(sample_replicate)
        filter3_list = filter3.get(sample_replicate)
        filter4_list = filter4.get(sample_replicate)
        variant_list = filter1_list + filter2_list + filter3_list + filter4_list
        data_to_drop = data_frame.loc[(data_frame['sample_replicate'] == sample_replicate) & data_frame['sequence'].isin(variant_list)]
        indexs_to_drop = data_to_drop.index.tolist()
        data_frame.drop(indexs_to_drop, inplace=True)
            # data_frame = data_frame.loc[
            #     (data_frame['sample_replicate'] != sample_replicate) &
            #     (data_frame['sequence'] != variant)
            # ]


def filter_fasta(engine, variant_model, variant_list, sample_fasta, chimera):
    with open(sample_fasta, 'w') as fout:
        for variant in variant_list:
            variant_select = select([variant_model.variant_id, variant_model.sequence]).where(variant_model.sequence == variant)
            variant_obj = engine.execute(variant_select)
            variant_request = list(variant_obj.fetchone())
            print("ok!")
            if chimera is True:
                fout.write(
                    ">{};size={};\n".format(variant_request[0], str(len(variant_request[1])))
                )
            else:
                fout.write(">{}\n".format(variant_request[0]))
            fout.write(variant_request[1] + '\n')


def pcr_error(engine, replicate_model, variant_model, data_frame,marker_id, var_prop, pcr_error_by_sample):
    """

    :param engine:
    :param replicate_model:
    :param variant_model:
    :param data_frame:
    :param marker_id:
    :param pcr_error_by_sample:
    :return:
    """
    select_replicate = select([replicate_model.biosample_name, replicate_model.name])\
        .where(replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(select_replicate)

    for replicate in replicate_obj:

        if pcr_error_by_sample is True:
            data_variant = data_frame.loc[data_frame['sample'] == replicate.biosample_name]
            sample_fasta_name = 'data/output/Filter/{}_{}.fasta'.format(replicate.biosample_name, marker_id)

        else:
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            data_variant = data_frame.loc[data_frame['sample_replicate'] == sample_replicate]
            sample_fasta_name = 'data/output/Filter/{}_{}.fasta'.format(sample_replicate, marker_id)

        variant_list_series = data_variant['sequence']
        if variant_list_series.empty is False:
            variant_list = sorted(set(list(variant_list_series)))
            filter_fasta(engine, variant_model, variant_list,sample_fasta_name, False)
            shortest_sequence = min(variant_list)
            L = len(shortest_sequence)
            id_raw = (L - 1)/ L
            id_rounded = round(id_raw, 2)
            sample_tsv_name = sample_fasta_name.replace('.fasta', '.tsv')
            subprocess.call(
                'vsearch --usearch_global ' + sample_fasta_name + ' --db ' + sample_fasta_name
                + ' --id ' + str(id_rounded) + ' --maxrejects 0 --maxaccepts 0 --userout ' + sample_tsv_name +
                ' --userfields query+target+alnlen+ids+mism+gaps --self', shell=True)

            with open(sample_tsv_name, 'r') as fin:

                for line in fin:
                    line = line.strip().split('\t')
                    query = line[0]
                    target = line[1]
                    mism = line[4]
                    gaps = line[5]

                    if (int(mism) + int(gaps)) == 0:
                        query_variant_select = select([variant_model.sequence]).where(variant_model.variant_id == query)
                        query_variant_obj = engine.execute(query_variant_select)
                        query_variant_request = list(query_variant_obj.fetchone())
                        query_variant_sequence = query_variant_request[0]
                        target_variant_select = select([variant_model.sequence]).where(variant_model.variant_id == target)
                        target_variant_obj = engine.execute(target_variant_select)
                        target_variant_request = list(target_variant_obj.fetchone())
                        target_variant_sequence = target_variant_request[0]
                        data_query_variant = data_variant.loc[data_variant['sequence'] == query_variant_sequence]
                        data_query_variant_series = data_query_variant['count']
                        query_variant_count = data_query_variant_series.sum()
                        data_target_variant = data_variant.loc[data_variant['sequence'] == target_variant_sequence]
                        data_target_variant_series = data_target_variant['count']
                        target_variant_count = data_target_variant_series.sum()
                        count_ratio = query_variant_count / target_variant_count
                        print(count_ratio)
                        if count_ratio < var_prop:
                            if pcr_error_by_sample is True:
                                data_to_drop = data_frame.loc[
                                    (data_frame['sequence'] == query_variant_sequence) &
                                    (data_frame['sample'] == replicate.biosample_name)
                                ]
                                indexs_to_drop = data_to_drop.index.tolist()
                                data_frame.drop(indexs_to_drop, inplace=True)
                            else:
                                data_to_drop = data_frame.loc[
                                    (data_frame['sequence'] == query_variant_sequence) &
                                    (data_frame['sample'] == sample_replicate)
                                    ]
                                indexs_to_drop = data_to_drop.index.tolist()
                                data_frame.drop(indexs_to_drop, inplace=True)


def chimera(engine, replicate_model, variant_model, data_frame, marker_id, chimera_by):
    """

    :param engine: Engine of the database
    :param replicate_model: Model of the replicate table
    :param variant_model:
    :param data_frame:
    :param marker_id:
    :param var_prop:
    :param chimera_by:
    :return:
    """
    select_replicate = select([replicate_model.biosample_name, replicate_model.name]) \
        .where(replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(select_replicate)

    for replicate in replicate_obj:

        if chimera_by == "sample":
            data_variant = data_frame.loc[data_frame['sample'] == replicate.biosample_name]
            repl_fasta_name = 'data/output/Filter/{}_{}_repl.fasta'.format(replicate.biosample_name, marker_id)

        elif chimera_by == "sample_replicate":
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            data_variant = data_frame.loc[data_frame['sample_replicate'] == sample_replicate]
            repl_fasta_name = 'data/output/Filter/{}_{}_repl.fasta'.format(sample_replicate, marker_id)

        sorted_repl_fasta_name = repl_fasta_name.replace('.fasta', '_sorted.fasta')
        variant_list_series = data_variant['sequence']
        if variant_list_series.empty is False:
            variant_list_series = data_variant['sequence']
            variant_list = sorted(set(list(variant_list_series)))
            print(variant_list)
            filter_fasta(engine, variant_model, variant_list, repl_fasta_name, True)
            subprocess.call('vsearch --sortbysize ' + repl_fasta_name + ' --output ' + sorted_repl_fasta_name, shell = True)
            # borderline_repl_filename = sorted_repl_fasta_name.replace('sorted', 'borderline')
            # nonchimera_repl_filename = sorted_repl_fasta_name.replace('sorted', 'correct')
            # chimera_repl_filename = sorted_repl_fasta_name.replace('sorted', 'chimera')
            # subprocess.call(
            #     'vsearch –uchime_denovo ' + sorted_repl_fasta_name + ' --borderline ' + borderline_repl_filename +
            #     ' --nonchimeras '+ nonchimera_repl_filename + ' --chimeras ' + chimera_repl_filename, shell=True
            # )
