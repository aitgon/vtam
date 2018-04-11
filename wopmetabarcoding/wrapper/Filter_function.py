from sqlalchemy import select, delete
import pandas, subprocess


def lfn_per_replicate(engine, replicate_model, variant_model, marker_id, data_frame):
    """
    Filter out Low Frequency Noise which eliminate a variant from sample-replicate if Nvar-repl/Nrepl < lfn_repl
    :param engine: Engine of the database
    :param replicate_model: Model of the replicate table
    :param variant_model: Model of the variant table
    :param marker_id: Name of the marker
    :param data_frame: Dataframe containing the content of the sample_count tsv
    :return: list failed_variants
    """
    # Selecting lines of the replicate table
    replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
        replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(replicate_select)
    failed_variants = {}
    for replicate in replicate_obj:
        failed_variants_list = []
        sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
        # Selecting lines of the variant table
        variant_select = select([variant_model.variant_id, variant_model.sequence]) \
            .where(variant_model.marker == marker_id)
        variant_obj = engine.execute(variant_select)
        # Creating a data frame with only the columns containing the sample replicate
        data_replicate = data_frame.loc[data_frame['sample_replicate'] == sample_replicate]
        # Checking if the data frame is not empty
        if data_replicate.empty is False:
            replicate_count_series = data_replicate['count']
            replicate_count = replicate_count_series.sum()
        for variant in variant_obj:
            data_variant = data_replicate.loc[data_replicate['sequence'] == variant.sequence]
            if data_variant.empty is False and data_replicate.empty is False:
                variant_replicate_count_series = data_variant['count']
                variant_replicate_count = variant_replicate_count_series.sum()
                lfn_value = variant_replicate_count / replicate_count
                if lfn_value < 0.001 and variant.variant_id not in failed_variants_list:
                    failed_variants_list.append(variant.variant_id)
        failed_variants[sample_replicate] = failed_variants_list
    return failed_variants


def lfn_per_variant(engine, replicate_model, variant_model, marker_id, data_frame, series_replicate):
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
                if series_replicate is True:
                    lfn_value = variantreplicate_count / replicate_count
                else:
                    lfn_value = variantreplicate_count / variant_count
                if lfn_value < 0.001 and variant.variant_id not in failed_variants_list:
                    failed_variants_list.append(variant.variant_id)
            failed_variants[variant.variant_id] = failed_variants_list
    return failed_variants


def lfn_per_readcounts(engine, replicate_model, variant_model, marker_id, custom_filter, data_frame):
    """

    :param engine:
    :param replicate_model:
    :param variant_model:
    :param marker_id:
    :param custom_filter:
    :param data_frame:
    :return:
    """
    replicate_select = select(
        [replicate_model.biosample_name, replicate_model.name]).where(replicate_model.marker_id == marker_id)
    replicate_obj = engine.execute(replicate_select)
    failed_variants = {}
    for replicate in replicate_obj:
        sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
        failed_variants_list = []
        variant_select = select([variant_model.variant_id, variant_model.sequence]) \
            .where(variant_model.marker == marker_id)
        variant_obj = engine.execute(variant_select)
        data_replicate = data_frame.loc[data_frame['sample_replicate'] == sample_replicate]
        for variant in variant_obj:
            data_variant = data_replicate.loc[data_replicate['sequence'] == variant.sequence]
            if data_variant.empty is False and data_replicate.empty is False:
                variant_replicate_count_series = data_variant['count']
                variant_replicate_count = variant_replicate_count_series.sum()
                if variant_replicate_count < custom_filter and variant.variant_id not in failed_variants_list:
                    failed_variants_list.append(variant.variant_id)
        failed_variants[sample_replicate] = failed_variants_list
    return failed_variants


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
                    failed_variants_list.append(variant.variant_id)
            failed_variants[variant.variant_id] = failed_variants_list
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


def min_repln(engine, variant_model, replicate_model, marker_id, data_frame):
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
            if repl_count >= 2 and variant.variant_id:
                print(variant.variant_id, replicate.biosample_name)


def min_replp(engine, variant_model, replicate_model, marker_id, data_frame, replicate_number):
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_id)
    variant_obj = engine.execute(variant_select)
    failed_variant = []
    for variant in variant_obj:
        data_variant = data_frame.loc[data_frame['sequence'] == variant.sequence]
        # Selecting lines of the replicate table
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_id == marker_id)
        replicate_obj = engine.execute(replicate_select)
        for replicate in replicate_obj:
            data_sample = data_variant.loc[data_variant['sample'] == replicate.biosample_name]
            repl_count = len(data_sample.index)
            if repl_count >= ((1/3)*replicate_number) and variant.variant_id not in failed_variant:
                failed_variant.append(variant.variant_id)
    return failed_variant


def delete_filtered_variants(session, variant_model, failed_variants):
    """
    Function which delete variants which don't pass lfn filters
    :param engine:
    :param variant_model:
    :param failed_variants:
    :return:
    """
    if len(failed_variants) != 0:
        for variant in failed_variants:
            stmt = session.query(variant_model).filter(variant_model.variant_id == variant)
            stmt.delete()
    else:
        pass


def pcr_error_fasta(engine, variant_model, variant_list, sample_fasta):
    with open(sample_fasta, 'w') as fout:
        for variant in variant_list:
            variant_select = select([variant_model.variant_id, variant_model.sequence]).where(variant_model.sequence == variant)
            variant_obj = engine.execute(variant_select)
            variant_request = list(variant_obj.fetchone())
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
            pcr_error_fasta(engine, variant_model, variant_list,sample_fasta_name)
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
                        if count_ratio < var_prop:
                            if pcr_error_by_sample is True:
                                data_frame = data_frame.loc[
                                    (data_frame['sequence'] != query_variant_sequence) &
                                    (data_frame['sample'] != replicate.biosample_name)
                                ]
                            else:
                                data_frame = data_frame.loc[
                                    (data_frame['sequence'] != query_variant_sequence) &
                                    (data_frame['sample'] != sample_replicate)
                                    ]
    return data_frame








