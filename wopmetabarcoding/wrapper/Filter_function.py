from sqlalchemy import select
import pandas


def filter1(engine, replicate_model, variant_model, marker_name, sample_count_tsv, data_frame):
    """

    :param engine:
    :param replicate_model:
    :param variant_model:
    :param marker_name:
    :param sample_count_tsv:
    :return:
    """
    # Selecting lines of the replicate table
    replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
        replicate_model.marker_name == marker_name)
    replicate_obj = engine.execute(replicate_select)
    succeed_variants = {}
    for replicate in replicate_obj:
        succeed_variants_list = []
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
            data_variant = data_replicate.loc[data_frame['sequence'] == variant.sequence]
            if data_variant.empty is False and data_replicate.empty is False:
                variant_replicate_count_series = data_variant['count']
                variant_replicate_count = variant_replicate_count_series.sum()
                lfn_value = variant_replicate_count / replicate_count
                if lfn_value > 0.001:
                    succeed_variants_list.append(variant.variant_id)
        succeed_variants[sample_replicate] = succeed_variants_list
    return succeed_variants


def filter2(engine, replicate_model, variant_model, marker_name, sample_count_tsv, data_frame):
    # Selecting lines of the variant table
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_name)
    variant_obj = engine.execute(variant_select)
    succeed_variants = {}
    for variant in variant_obj:
        # Selecting lines of the replicate table
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_name == marker_name)
        replicate_obj = engine.execute(replicate_select)
        data_variant = data_frame.loc[data_frame['sequence'] == variant.sequence]
        succeed_variants_list = []
        if data_variant.empty is False:
            variant_replicate_count_series = data_variant['count']
            variant_replicate_count = variant_replicate_count_series.sum()
        for replicate in replicate_obj:
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            data_replicate = data_variant.loc[data_frame['sample_replicate'] == sample_replicate]
            if data_replicate.empty is False and data_variant.empty is False:
                replicate_count_series = data_replicate['count']
                replicate_count = replicate_count_series.sum()
                lfn_value = replicate_count / variant_replicate_count
                print(lfn_value)
                if lfn_value > 0.001:
                    succeed_variants_list.append(variant.variant_id)
        succeed_variants[variant.variant_id] = succeed_variants_list
    return succeed_variants