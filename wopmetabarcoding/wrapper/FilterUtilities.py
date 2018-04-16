from sqlalchemy import select, delete
from wopmetabarcoding.utils.VSearch import VSearch1, Vsearch2, Vsearch3
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
    df2 = df2.merge(cutoff_df, left_on="sequence", right_on="sequence")
    return list(df2.ix[df2.low_frequence_noice_per_variant < df2.cutoff].index) # select rare reads to discard later


def lfn4_per_replicate_series_with_cutoff(df, cutoff_tsv):
    cutoff_df = pandas.read_csv(cutoff_tsv, sep="\t", header=0)
    cutoff_df.columns = ['sequence', 'cutoff']
    #
    df2 = df.groupby(by=['replicate']).sum()
    df2 = df.merge(df2, left_on='replicate', right_index=True)
    df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_replicate_series']
    df2['low_frequence_noice_per_replicate_series'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_replicate_series
    #
    # merge with cutoff
    df2 = df2.merge(cutoff_df, left_on="sequence", right_on="sequence")
    return list(df2.ix[df2.low_frequence_noice_per_replicate_series < df2.cutoff].index) # select rare reads to discard later


def delete_filtered_variants(df, filter1, filter2, filter3, filter4):
    """
    Function deleting all variants which didn't pass the LFNs filters
    :param df:
    :param filter1:
    :param filter2:
    :param filter3:
    :param filter4:
    :return:
    """
    failed_variants_lfn_filters = filter1 + filter2 + filter3 + filter4
    failed_variants_lfn_filters = sorted(set(failed_variants_lfn_filters))
    df.drop(failed_variants_lfn_filters, inplace=True)


def min_repln(engine, variant_model, replicate_model, marker_id, data_frame, min_count):
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker_id == marker_id)
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
        .where(variant_model.marker_id == marker_id)
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


def filter_fasta(engine, variant_model, variant_list, sample_fasta, chimera):
    with open(sample_fasta, 'w') as fout:
        for variant in variant_list:
            variant_select = select([variant_model.variant_id, variant_model.sequence]).where(variant_model.sequence == variant)
            variant_obj = engine.execute(variant_select)
            variant_request = list(variant_obj.fetchone())
            if chimera is True:
                fout.write(
                    ">{};size={};\n".format(variant_request[0], str(len(variant_request[1])))
                )
            else:
                fout.write(">{}\n".format(variant_request[0]))
            fout.write(variant_request[1] + '\n')


def selecting_target_query(engine, variant_model, target_id, query_id):
    query_variant_select = select([variant_model.sequence]).where(variant_model.variant_id == query_id)
    query_variant_obj = engine.execute(query_variant_select)
    query_variant_request = list(query_variant_obj.fetchone())
    query_variant_sequence = query_variant_request[0]
    target_variant_select = select([variant_model.sequence]).where(variant_model.variant_id == target_id)
    target_variant_obj = engine.execute(target_variant_select)
    target_variant_request = list(target_variant_obj.fetchone())
    target_variant_sequence = target_variant_request[0]
    return [query_variant_sequence, target_variant_sequence]


def pcr_error(engine, replicate_model, variant_model, df, marker_id, var_prop, pcr_error_by_sample):
    """

    :param engine:
    :param replicate_model:
    :param variant_model:
    :param df:
    :param marker_id:
    :param pcr_error_by_sample:
    :return:
    """
    select_replicate = select([replicate_model.biosample_name, replicate_model.name])\
        .where(replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(select_replicate)

    for replicate in replicate_obj:

        if pcr_error_by_sample is True:
            data_variant = df.loc[df['sample'] == replicate.biosample_name]
            sample_fasta_name = 'data/output/Filter/{}_{}.fasta'.format(replicate.biosample_name, marker_id)

        else:
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            data_variant = df.loc[df['sample_replicate'] == sample_replicate]
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
            vsearch_usearch_global_args = {'db': sample_fasta_name,
                              'usearch_global': sample_fasta_name,
                              'id': id_rounded,
                              'maxrejects': 0,
                              'maxaccepts': 0,
                              'userfields': "query+target+alnlen+ids+mism+gaps --self'",
                              'userout': sample_tsv_name,
                              }
            vsearch_usearch_global = VSearch1(**vsearch_usearch_global_args)
            vsearch_usearch_global.run()
            del vsearch_usearch_global

            with open(sample_tsv_name, 'r') as fin:

                for line in fin:
                    line = line.strip().split('\t')
                    query = line[0]
                    target = line[1]
                    mism = line[4]
                    gaps = line[5]

                    if (int(mism) + int(gaps)) == 0:
                        sequence_target_query = selecting_target_query(engine, variant_model, query, target)
                        query_sequence = sequence_target_query[0]
                        target_sequence = sequence_target_query[1]
                        data_query_variant = data_variant.loc[data_variant['sequence'] == query_sequence]
                        query_count = data_query_variant['count'].sum()
                        data_target_variant = data_variant.loc[data_variant['sequence'] == target_sequence]
                        target_count = data_target_variant['count'].sum()
                        count_ratio = query_count / target_count

                        if count_ratio < var_prop:
                            if pcr_error_by_sample is True:
                                data_to_drop = df.loc[
                                    (df['sequence'] == query_sequence) &
                                    (df['sample'] == replicate.biosample_name)
                                    ]
                                indexs_to_drop = data_to_drop.index.tolist()
                                df.drop(indexs_to_drop, inplace=True)


def chimera(engine, replicate_model, variant_model, df, marker_id, chimera_by):
    """

    :param engine: Engine of the database
    :param replicate_model: Model of the replicate table
    :param variant_model:
    :param df:
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
            df2 = df.loc[df['sample'] == replicate.biosample_name]
            repl_fasta_name = 'data/output/Filter/{}_{}_repl.fasta'.format(replicate.biosample_name, marker_id)

        elif chimera_by == "sample_replicate":
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            df2 = df.loc[df['sample_replicate'] == sample_replicate]
            repl_fasta_name = 'data/output/Filter/{}_{}_repl.fasta'.format(sample_replicate, marker_id)

        sorted_repl_fasta_name = repl_fasta_name.replace('.fasta', '_sorted.fasta')
        variant_list_series = df2['sequence']
        if variant_list_series.empty is False:
            variant_list_series = df2['sequence']
            variant_list = sorted(set(list(variant_list_series)))
            filter_fasta(engine, variant_model, variant_list, repl_fasta_name, True)
            vsearch_sortbysize_args = {"sortbysize": repl_fasta_name, "output": sorted_repl_fasta_name}
            vsearch_sortbysize = Vsearch2(**vsearch_sortbysize_args)
            vsearch_sortbysize.run()
            del vsearch_sortbysize
            borderline_repl_filename = sorted_repl_fasta_name.replace('sorted', 'borderline')
            nonchimera_repl_filename = sorted_repl_fasta_name.replace('sorted', 'correct')
            chimera_repl_filename = sorted_repl_fasta_name.replace('sorted', 'chimera')
            vsearch_chimera_args = {
                "uchime_denovo": sorted_repl_fasta_name, "borderline": borderline_repl_filename,
                "nonchimeras": nonchimera_repl_filename, "chimeras": chimera_repl_filename
            }
            vsearch_chimera = Vsearch3(**vsearch_chimera_args)
            vsearch_chimera.run()
            del vsearch_chimera

            with open(chimera_repl_filename, 'r') as fin:
                for line in fin:

                    if ">" in line:
                        line = line.replace('>', '')
                        variant_id = line.strip().split(';')[0]
                        variant_select = select([variant_model.sequence]).where(variant_model.variant_id == variant_id)
                        variant_obj = engine.execute(variant_select)
                        variant = variant_obj.fetchone()
                        df_variant = df2.loc[df2['sequence'] == variant.sequence]

                        if df_variant.empty is False:
                            indexs_to_drop = df_variant.index.tolist()
                            df.drop(indexs_to_drop, inplace=True)
                        else:
                            pass



