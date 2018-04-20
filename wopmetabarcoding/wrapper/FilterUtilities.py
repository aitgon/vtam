from sqlalchemy import select
from wopmetabarcoding.utils.VSearch import VSearch1, Vsearch2, Vsearch3
import pandas, itertools
from Bio import SeqIO
import os
from wopmetabarcoding.utils.constants import tempdir


def lfn1_per_replicate(df, lfn_per_replicate_threshold):
    """
    Function calculating the Low Frequency Noise per sample replicate
    :param df: dataframe containing the information
    :param lfn_per_replicate_threshold: threshold defined by the user
    :return: List of the index which don't pass the filter
    """
    df2 = df.groupby(by=['sample_replicate']).sum()
    df2 = df.merge(df2, left_on='sample_replicate', right_index=True)
    df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_sample_replicate']
    df2['low_frequence_noice_per_replicate'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_sample_replicate
    return list(df2.ix[df2.low_frequence_noice_per_replicate < lfn_per_replicate_threshold].index) # select rare reads to discard later


def lfn2_per_variant(df, lfn_per_variant_threshold):
    """
    Function calculating the Low Frequency Noise per variant
    :param df: dataframe containing the information
    :param lfn_per_variant_threshold: threshold defined by the user
    :return: List of the index which don't pass the filter
    """
    df2 = df.groupby(by=['sample_replicate']).sum()
    df2 = df.merge(df2, left_on='sample_replicate', right_index=True)
    df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_variant']
    df2['low_frequence_noice_per_variant'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_variant
    return list(df2.ix[df2.low_frequence_noice_per_variant < lfn_per_variant_threshold].index) # select rare reads to discard later


def lfn2_per_replicate_series(df, lfn_per_replicate_series_threshold):
    """
    Function calculating the Low Frequency Noise per replicate series
    :param df: dataframe containing the information
    :param lfn_per_replicate_series_threshold: threshold defined by the user
    :return: List of the index which don't pass the filter
    """
    df2 = df.groupby(by=['replicate']).sum()
    df2 = df.merge(df2, left_on='replicate', right_index=True)
    df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_replicate_series']
    df2['low_frequence_noice_per_replicate_series'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_replicate_series
    return list(df2.ix[df2.low_frequence_noice_per_replicate_series < lfn_per_replicate_series_threshold].index) # select rare reads to discard later


def lfn3_read_count(df, lfn_read_count_threshold):
    """
    Function calculating the Low Frequency Noise per users defined minimal readcount
    :param df: dataframe containing the information
    :param lfn_read_count_threshold: threshold defined by the user
    :return: List of the index which don't pass the filter
    """
    return list(df.ix[df['count'] <= lfn_read_count_threshold].index) # select rare reads to discard later


def lfn4_per_variant_with_cutoff(df, cutoff_tsv):
    """
    Function calculating the Low Frequency Noise per variant against cutoff
    :param df: dataframe containing the information
    :param cutoff_tsv: file containing the cutoffs for each variant
    :return: List of the index which don't pass the filter
    """
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
    """
        Function calculating the Low Frequency Noise per replicate series against cutoff
        :param df: dataframe containing the information
        :param cutoff_tsv: file containing the cutoffs for each variant
        :return: List of the index which don't pass the filter
        """
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
    :param df: dataframe containing the information
    :param filter1: List of the index which don't pass the filter 1
    :param filter2: List of the index which don't pass the filter 2
    :param filter3: List of the index which don't pass the filter 3
    :param filter4: List of the index which don't pass the filter 4
    :return: None
    """
    failed_variants_lfn_filters = filter1 + filter2 + filter3 + filter4
    failed_variants_lfn_filters = sorted(set(failed_variants_lfn_filters))
    df.drop(failed_variants_lfn_filters, inplace=True)


def min_repln(engine, replicate_model, marker_id, df, min_count):
    """
    Filter watching if a variant is present the minimal number of sample replicates of a sample
    :param engine: engine of the database
    :param replicate_model: model of the replicate table
    :param marker_id: id of the marker
    :param df: dataframe containing the information
    :param min_count: minimal number of replicates in which the variant must be present
    :return: None
    """
    replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
                replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(replicate_select)
    for replicate in replicate_obj:
        df2 = df.loc[df['sample'] == replicate.biosample_name]
        df3 = df2['sequence'].value_counts().to_frame()
        df3.columns = ['sequence_count']
        df2 = df2.merge(df3, left_on='sequence', right_index=True)
        index_to_drop = list(df2.ix[df2.sequence_count < min_count].index)
        df.drop(index_to_drop, inplace=True)


def min_replp(engine, replicate_model, marker_id, df, replicate_count):
    """
    Filter watching if a variant is present the more more than 1/3 of the number of sample replicates of a sample
    :param engine: engine of the database
    :param replicate_model: model of the replicate table
    :param marker_id: id of the marker
    :param df: data frame containing the information
    :param replicate_count: number of sample_replicate in each
    :return: None
    """
    replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
                replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(replicate_select)
    for replicate in replicate_obj:
        df2 = df.loc[df['sample'] == replicate.biosample_name]
        df3 = df2['sequence'].value_counts().to_frame()
        df3.columns = ['sequence_count']
        df2 = df2.merge(df3, left_on='sequence', right_index=True)
        index_to_drop = list(df2.ix[df2.sequence_count < ((1/3) * replicate_count)].index)
        df.drop(index_to_drop, inplace=True)


def filter_fasta(engine, variant_model, variant_list, sample_fasta, chimera):
    """
    Function used to create fasta file from the dataframe
    :param engine: engine of the database
    :param variant_model: model of the variant table
    :param variant_list: list of the variant present in the data frame
    :param sample_fasta: file name choose to name the fasta file
    :param chimera: choose if the function is used for the pcr_error or the chimera step
    :return: None
    """
    with open(sample_fasta, 'w') as fout:
        for variant in variant_list:
            variant_select = select([variant_model.variant_id, variant_model.sequence]).where(variant_model.sequence == variant)
            variant_obj = engine.execute(variant_select)
            variant_request = list(variant_obj.fetchone())
            if chimera is True:
                fout.write(
                    ">{};size={};\n".format(variant_request[1], str(len(variant_request[1])))
                )
            else:
                fout.write(">{}\n".format(variant_request[1]))
            fout.write(variant_request[1] + '\n')


def filter_fastabis(variant_list, sample_fasta, chimera):
    """
    Function used to create fasta file from the dataframe
    :param variant_list: list of the variant present in the data frame
    :param sample_fasta: file name choose to name the fasta file
    :param chimera: choose if the function is used for the pcr_error or the chimera step
    :return: None
    """
    with open(sample_fasta, 'w') as fout:
        for variant in variant_list:
            if chimera is True:
                fout.write(
                    ">{};size={};\n".format(variant, str(len(variant)))
                )
            else:
                fout.write(">{}\n".format(variant))
            fout.write(variant + '\n')


def pcr_error(engine, replicate_model, variant_model, df, marker_id, var_prop, pcr_error_by_sample):
    """
    Function used to eliminate variants from the data frame in which a pcr error is spotted
    :param engine: engine of the database
    :param replicate_model: model of the replicate table
    :param variant_model: model of the variants table
    :param df: data frame containing the information
    :param marker_id: id of the marker
    :param var_prop: threshold choose by the user
    :param pcr_error_by_sample: boolean used to choose if the analysis is made by sample or sample_replicate
    :return: None
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
                              'id': str(id_rounded),
                              'maxrejects': 0,
                              'maxaccepts': 0,
                              'userout': sample_tsv_name,
                              'userfields': "query+target+alnlen+ids+mism+gaps",
                                           }
            vsearch_usearch_global = VSearch1(**vsearch_usearch_global_args)
            vsearch_usearch_global.run()
            column_names = ['query', 'target', 'alnlen', 'ids', 'mism', 'gaps']
            df2 = pandas.read_csv(sample_tsv_name, sep='\t', names = column_names)
            false_df2 = list(df2.ix[(df2.mism + df2.gaps) != 1].index)
            df2.drop(false_df2, inplace=True)
            df3 = df2[['query', 'target']]
            df4 = data_variant[['sequence', 'count']]
            df3 = df3.merge(df4, left_on=['query'], right_on=['sequence'])
            df3 = df3.merge(df4, left_on=['target'], right_on=['sequence'])
            df3['count_ratio'] = df3.count_x / df3.count_y
            df3 = df3[['query', 'count_ratio']]
            df3.index = df3['query']
            data_variant = data_variant.merge(df3, left_on='sequence', right_index=True)
            index_to_drop = list(data_variant.ix[data_variant.count_ratio < var_prop].index)
            df.drop(index_to_drop, inplace=True)


def chimera(replicate_obj_list, df, marker_id, chimera_by):
    """
    Function using Vsearch to identify chimera or borderline variant
    :param replicate_obj_list: list of the sample and the corresponding sample replicate
    :param df: data frame containing the information
    :param marker_id: id of the marker
    :param chimera_by: parameter to choose between sample, sample replicate or all
    :return: None
    """

    for replicate in replicate_obj_list:

        if chimera_by == "sample":
            df2 = df.loc[df['sample'] == replicate.biosample_name]
            filename = '{}_{}_repl.fasta'.format(replicate.biosample_name, marker_id)
            repl_fasta_name = os.path.join(tempdir, filename)

        elif chimera_by == "sample_replicate":
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            df2 = df.loc[df['sample_replicate'] == sample_replicate]
            # repl_fasta_name = 'data/output/Filter/{}_{}_repl.fasta'.format(sample_replicate, marker_id)
            filename = '{}_{}_repl.fasta'.format(sample_replicate, marker_id)
            repl_fasta_name = os.path.join(tempdir, filename)

        sorted_repl_fasta_name = repl_fasta_name.replace('.fasta', '_sorted.fasta')
        variant_list_series = df2['sequence']
        if variant_list_series.empty is False:
            variant_list_series = df2['sequence']
            variant_list = sorted(set(list(variant_list_series)))
            filter_fastabis(variant_list, repl_fasta_name, True)
            vsearch_sortbysize_args = {"sortbysize": repl_fasta_name, "output": sorted_repl_fasta_name}
            vsearch_sortbysize = Vsearch2(**vsearch_sortbysize_args)
            vsearch_sortbysize.run()
            del vsearch_sortbysize
            borderline_repl_filename = sorted_repl_fasta_name.replace('sorted', 'borderline')
            nonchimera_repl_filename = sorted_repl_fasta_name.replace('sorted', 'correct')
            chimera_repl_filename = sorted_repl_fasta_name.replace('sorted', 'chimera')
            vsearch_chimera_args = {
                "uchime_denovo": sorted_repl_fasta_name,
                "borderline": borderline_repl_filename,
                "nonchimeras": nonchimera_repl_filename,
                "chimeras": chimera_repl_filename
            }
            vsearch_chimera = Vsearch3(**vsearch_chimera_args)
            vsearch_chimera.run()
            del vsearch_chimera
            with open(chimera_repl_filename, 'r') as fin:
                for line in fin:

                    if ">" in line:
                        line = line.replace('>', '')
                        variant_sequence = line.strip().split(';')[0]
                        df_variant = df2.loc[df2['sequence'] == variant_sequence]

                        if df_variant.empty is False:
                            indexs_to_drop = df_variant.index.tolist()
                            df.drop(indexs_to_drop, inplace=True)
                        else:
                            pass
            borderline_variants = []
            for record in SeqIO.parse(borderline_repl_filename, 'fasta'):
                borderline_variants.append(record.description.strip().split(';')[0])
            df['is_borderline'] = (df['sequence'] in borderline_variants)
            print(df)


def renkonen(engine, replicate_model, df, marker_id, number_of_replicate):
    """
    Function calculating distance between sample replicate of the sample to determine if the sample replicate
    must be deleted or not
    :param engine: engine of the database
    :param replicate_model: model of the replicate table
    :param df: data frame containing the information
    :param marker_id: id of the marker
    :param number_of_replicate: Number of replicate by sample
    :return: None
    """
    df2 = df.groupby(by=['sample_replicate']).sum()
    df2 = df.merge(df2, left_on='sample_replicate', right_index=True)
    select_replicate = select([replicate_model.biosample_name, replicate_model.name]) \
        .where(replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(select_replicate)
    for replicate in replicate_obj:
        df_replicate = df2.loc[df['sample'] == replicate.biosample_name]
        replicates = df_replicate['sample_replicate'].tolist()
        for combi in itertools.combinations(replicates, 2):
            combi = list(combi)
            df_repl1 = df2.loc[df2['sample_replicate'] == combi[0]]
            df_repl2 = df2.loc[df2['sample_replicate'] == combi[1]]
            repl1_variant = df_repl1[['sequence', 'sample_replicate', 'count_x', 'count_y']]
            repl2_variant = df_repl2[['sequence', 'sample_replicate', 'count_x', 'count_y']]
            # repl1_2_variant = pandas.concat((repl1_variant, repl2_variant), axis=1)
            repl1_2_variant = repl1_variant.append(repl2_variant)
            repl1_2_variant['ratio'] = repl1_2_variant['count_x'] / repl1_2_variant['count_y']
            group_repl = repl1_2_variant.groupby(by=['sequence']).min()
            ratios = group_repl['ratio']
            distance = 1 - group_repl['ratio'].sum()
            if distance < ((number_of_replicate - 1)/2):
                sample_replicate_rm = df.loc[(df['sample'] == replicate.biosample_name) & (df['sample_replicate'] == combi[0])]
                # index_to_drop = sample_replicate_rm
                # df.drop(index_to_drop, inplace=True)
                index_to_drop = sample_replicate_rm.index.tolist()
                df.drop(index_to_drop, inplace=True)
















