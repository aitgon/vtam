from sqlalchemy import select


def lfn_per_replicate_filter1(engine, replicate_model, variant_model, marker_name, sample_count_tsv):
    # Selecting lines of the replicate table
    replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
        replicate_model.marker_name == marker_name)
    replicate_obj = engine.execute(replicate_select)
    failed_variants = {}
    for replicate in replicate_obj:
        failed_variants_list = []
        sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
        print(sample_replicate)
        # Selecting lines of the variant table
        variant_select = select([variant_model.variant_id, variant_model.sequence]) \
            .where(variant_model.marker == marker_name)
        variant_obj = engine.execute(variant_select)
        for variant in variant_obj:
            count_variant = 0
            replicate_count = 0
            with open(sample_count_tsv, 'r') as fin:
                for line in fin:
                    if sample_replicate in line and variant.sequence in line:
                        count_variant += int(line.strip().split('\t')[3])
                    if sample_replicate in line:
                        replicate_count += int(line.strip().split('\t')[3])
                if replicate_count != 0:
                    lfn_value = count_variant/replicate_count
                    if lfn_value < 0.001:
                        failed_variants_list.append(variant.variant_id)
        failed_variants[sample_replicate] = failed_variants_list
    return failed_variants_list


def lfn_per_variant_filter2(engine, replicate_model, variant_model, marker_name, sample_count_tsv):
    # Selecting lines of the variant table
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_name)
    variant_obj = engine.execute(variant_select)
    failed_variants = {}
    for variant in variant_obj:
        # Selecting lines of the replicate table
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_name == marker_name)
        replicate_obj = engine.execute(replicate_select)
        failed_variants_list = []
        for replicate in replicate_obj:
            sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
            with open(sample_count_tsv, 'r') as fin:
                count_variant = 0
                variant_count = 0
                for line in fin:
                    if sample_replicate in line and variant.sequence in line:
                        count_variant += int(line.strip().split('\t')[3])
                    if sample_replicate in line:
                        variant_count += int(line.strip().split('\t')[3])
                if variant_count != 0:
                    lfn_value = count_variant / variant_count
                    if lfn_value < 0.001:
                        failed_variants_list.append(variant.variant_id)
        failed_variants[variant.variant_id] = failed_variants_list
    return failed_variants_list



