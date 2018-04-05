from sqlalchemy import select
import pandas


def lfn_per_replicate(engine, replicate_model, variant_model, marker_name, data_frame):
    """
    Filter out Low Frequency Noise which eliminate a variant from sample-replicate if Nvar-repl/Nrepl < lfn_repl
    :param engine: Engine of the database
    :param replicate_model: Model of the replicate table
    :param variant_model: Model of the variant table
    :param marker_name: Name of the marker
    :param data_frame: Dataframe containing the content of the sample_count tsv
    :return: list failed_variants
    """
    # Selecting lines of the replicate table
    replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
        replicate_model.marker_name == marker_name)
    replicate_obj = engine.execute(replicate_select)
    failed_variants = {}
    for replicate in replicate_obj:
        failed_variants_list = []
        sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
        # Selecting lines of the variant table
        variant_select = select([variant_model.variant_id, variant_model.sequence]) \
            .where(variant_model.marker == marker_name)
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


def lfn_per_variant(engine, replicate_model, variant_model, marker_name, data_frame, replicate_series):
    """
    Filter out Low Frequency Noise which eliminate a variant from sample-replicate if Nvar-repl/Nvar < lfn_var
    :param engine: Engine of the database
    :param replicate_model: Model of the replicate table
    :param variant_model: Model of the variant table
    :param marker_name: Name of the marker
    :param data_frame: Dataframe containing the content of the sample_count tsv
    :return: list failed_variants
    """
    # Selecting lines of the variant table
    print(replicate_series)
    variantreplicate_count = 0
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_name)
    variant_obj = engine.execute(variant_select)
    failed_variants = {}
    for variant in variant_obj:
        # Selecting lines of the replicate table
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_name == marker_name)
        replicate_obj = engine.execute(replicate_select)
        data_variant = data_frame.loc[data_frame['sequence'] == variant.sequence]
        failed_variants_list = []
        if replicate_series is False:
            if data_variant.empty is False:
                variant_count_series = data_variant['count']
                variant_count = variant_count_series.sum()
                print(variant.variant_id)
                print(variant_count)
        for replicate in replicate_obj:
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            if replicate_series is True:
                data_replicate = data_frame.loc[data_frame['replicate'] == replicate.replicate_name]
                if data_replicate.empty is False:
                    replicate_count_series = data_replicate['count']
                    replicate_count = replicate_count_series.sum()
            data_variantreplicate = data_variant.loc[data_variant['sample_replicate'] == sample_replicate]
            if data_variantreplicate.empty is False:
                variantreplicate_count_series = data_variantreplicate['count']
                variantreplicate_count = variantreplicate_count_series.sum()
                print(variantreplicate_count)
            if replicate_series is True:
                print(variantreplicate_count)
                lfn_value = variantreplicate_count / replicate_count
            else:
                print(variantreplicate_count)
                lfn_value = variantreplicate_count / variant_count
            if lfn_value < 0.001 and variant.variant_id not in failed_variants_list:
                failed_variants_list.append(variant.variant_id)
        failed_variants[variant.variant_id] = failed_variants_list
    return failed_variants


def lfn_per_readcounts(engine, replicate_model, variant_model, marker_name, custom_filter, data_frame):
    """

    :param engine:
    :param replicate_model:
    :param variant_model:
    :param marker_name:
    :param custom_filter:
    :param data_frame:
    :return:
    """
    replicate_select = select(
        [replicate_model.biosample_name, replicate_model.name]).where(replicate_model.marker_name == marker_name)
    replicate_obj = engine.execute(replicate_select)
    failed_variants = {}
    for replicate in replicate_obj:
        sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
        failed_variants_list = []
        variant_select = select([variant_model.variant_id, variant_model.sequence]) \
            .where(variant_model.marker == marker_name)
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


def lfn_per_cutoff(engine, replicate_model, variant_model, marker_name, data_frame, replicate_series, cutoff_tsv):
    """
        Filter out Low Frequency Noise which eliminate a variant from sample-replicate if Nvar-repl/Nvar < lfn_var
        :param engine: Engine of the database
        :param replicate_model: Model of the replicate table
        :param variant_model: Model of the variant table
        :param marker_name: Name of the marker
        :param data_frame: Dataframe containing the content of the sample_count tsv
        :return: list failed_variants
        """
    # Selecting lines of the variant table
    variantreplicate_count = 0
    data_cutoff = create_cutoff_table(cutoff_tsv)
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_name)
    variant_obj = engine.execute(variant_select)
    failed_variants = {}
    for variant in variant_obj:
        # Selecting lines of the replicate table
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_name == marker_name)
        replicate_obj = engine.execute(replicate_select)
        data_variant = data_frame.loc[data_frame['sequence'] == variant.sequence]
        failed_variants_list = []
        if replicate_series is False:
            if data_variant.empty is False:
                variant_count_series = data_variant['count']
                variant_count = variant_count_series.sum()
        for replicate in replicate_obj:
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            if replicate_series is True:
                data_replicate = data_frame.loc[data_frame['replicate'] == replicate.replicate_name]
                if data_replicate.empty is False:
                    replicate_count_series = data_replicate['count']
                    replicate_count = replicate_count_series.sum()
            data_variantreplicate = data_variant.loc[data_variant['sample_replicate'] == sample_replicate]
            if data_variantreplicate.empty is False:
                variantreplicate_count_series = data_variantreplicate['count']
                variantreplicate_count = variantreplicate_count_series.sum()
            if replicate_series is True:
                lfn_value = variantreplicate_count / replicate_count
            else:
                lfn_value = variantreplicate_count / variant_count
            data_cutoff_filter = data_cutoff.loc[data_cutoff['sequence'] == variant.sequence]
            cutoff_series = data_cutoff_filter['value']
            cutoff = cutoff_series.sum()
            if lfn_value < cutoff and variant.variant_id not in failed_variants_list:
                failed_variants_list.append(variant.variant_id)
        failed_variants[variant.variant_id] = failed_variants_list
    return failed_variants


def min_repln(engine, sample_count_tsv, variant_list, variant_model, value):
    succeed_variants = []
    for variant in variant_list:
        variant_select = select([variant_model.variant_id, variant_model.sequence]).where(variant_model.variant_id == variant)
        variant_obj = engine.execute(variant_select)
        for table_variants in variant_obj:
            print(table_variants)
            with open(sample_count_tsv, 'r') as fin:
                count = 0
                for line in fin:
                    if table_variants.sequence in line:
                        count += 1
                if count >= value:
                    print(table_variants.variant_id)
                    print(count)
                    succeed_variants.append(variant)
    return succeed_variants
