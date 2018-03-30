from sqlalchemy import select
import sqlite3


def lfn_per_replicate_filter1(engine, replicate_model, variant_model, marker_name, sample_count_tsv):
    replicate_select = select([replicate_model.biosample_name,replicate_model.name]).where(replicate_model.marker_name == marker_name)
    replicate_obj = engine.execute(replicate_select)
    failed_variants = {}
    for replicate in replicate_obj:
        failed_variants_list = []
        variant_select = select([variant_model.variant_id, variant_model.sequence])\
            .where(variant_model.marker == marker_name)
        variant_obj = engine.execute(variant_select)
        sample_replicate = replicate.biosample_name + "_" + replicate.name
        for variant in variant_obj:
            count_variant = 0
            replicate_count = 0
            with open(sample_count_tsv, 'r') as fin:
                for line in fin:
                    variant_sequence = line.strip().split('\t')[0]
                    sample_replicate_name = line.strip().split('\t')[2]
                    if variant.sequence == variant_sequence and sample_replicate == sample_replicate_name:
                        count_variant += int(line.strip().split('\t')[3])
                    if sample_replicate == sample_replicate_name:
                        replicate_count += int(line.strip().split('\t')[3])
                try:
                    lfn_value = count_variant/replicate_count
                    print(sample_replicate)
                    print(lfn_value)
                    if lfn_value < 0.001:
                        failed_variants_list.append(variant.variant_id)
                except ZeroDivisionError:
                    pass
        if len(failed_variants_list) != 0:
            failed_variants[sample_replicate] = failed_variants_list
    return failed_variants


def lfn_per_variant_filter2(engine, replicate_model, variant_model, marker_name, sample_count_tsv):
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_name)
    variant_obj = engine.execute(variant_select)
    failed_variants = {}
    for variant in variant_obj:
        failed_variants_list = []
        replicate_select = select([replicate_model.biosample_name,replicate_model.name]).where(replicate_model.marker_name == marker_name)
        replicate_obj = engine.execute(replicate_select)
        for replicate in replicate_obj:
            count_variant = 0
            variant_count = 0
            sample_replicate = replicate.biosample_name + "_" + replicate.name
            with open(sample_count_tsv, 'r') as fin:
                for line in fin:
                    variant_sequence = line.strip().split('\t')[0]
                    sample_replicate_name = line.strip().split('\t')[2]
                    if variant.sequence == variant_sequence and sample_replicate == sample_replicate_name:
                        count_variant += int(line.strip().split('\t')[3])
                    if variant.sequence == variant_sequence:
                        variant_count += int(line.strip().split('\t')[3])
                try:
                    lfn_value = count_variant/variant_count
                    if lfn_value < 0.001:
                        failed_variants_list.append(variant.variant_id)
                except ZeroDivisionError:
                    pass
        failed_variants[sample_replicate] = failed_variants_list
    return failed_variants


def empirical_filter3(engine, replicatemarker_model, variant_model, marker_name, sample_count_tsv, filter):
    replicate_select = select([replicatemarker_model.name]).where(replicatemarker_model.marker_name == marker_name)
    replicate_obj = engine.execute(replicate_select)
    failed_variants = {}
    for replicate in replicate_obj:
        failed_variants_list = []
        variant_select = select([variant_model.variant_id, variant_model.sequence]) \
            .where(variant_model.marker == marker_name)
        variant_obj = engine.execute(variant_select)
        for variant in variant_obj:
            count_variant = 0
            with open(sample_count_tsv, 'r') as fin:
                for line in fin:
                    variant_sequence = line.strip().split('\t')[0]
                    replicate_name = line.strip().split('\t')[1]
                    if variant.sequence == variant_sequence and replicate.name == replicate_name:
                        count_variant += int(line.strip().split('\t')[3])
                if count_variant < filter and variant.variant_id not in failed_variants_list:
                    failed_variants_list.append(variant.variant_id)
        failed_variants[replicate.name] = failed_variants_list
    return failed_variants


def create_cutoff_table(conn, cutoff_file_tsv):
    with open(cutoff_file_tsv, 'r') as fin:
        for line in fin:
            conn.execute("INSERT INTO read_fasta (variant_id, cutoff) VALUES (?, ?)",
                         (line.strip().split('\t')[0], line.strip().split('\t')[1],))
    conn.commit()


# sc = Specific cutoff
def lfn_per_variant_sc_filter4(engine, replicatemarker_model, variant_model, marker_name, sample_count_tsv, cutoff_tsv):
    table_filename = cutoff_tsv.replace('.tsv', '.sqlite')
    conn = sqlite3.connect(table_filename)
    conn.execute("DROP TABLE IF EXISTS cutoff_table")
    conn.execute("CREATE TABLE  cutoff_table (variant_id VARCHAR PRIMARY KEY , cutoff INT)")
    create_cutoff_table(conn, cutoff_tsv)
    #
    variant_select = select([variant_model.variant_id, variant_model.sequence]) \
        .where(variant_model.marker == marker_name)
    variant_obj = engine.execute(variant_select)
    failed_variants = {}
    for variant in variant_obj:
        #
        failed_variants_list = []
        replicate_select = select([replicatemarker_model.name]).where(replicatemarker_model.marker_name == marker_name)
        replicate_obj = engine.execute(replicate_select)
        #
        cutoff_cursor = conn.execute('SELECT cutoff FROM cutoff WHERE variant_id=?', (variant.variant_id,))
        cutoff_list = list(cutoff_cursor.fetchone())
        cutoff_cursor.close()
        cutoff = cutoff_list[0]
        #
        for replicate in replicate_obj:
            count_variant = 0
            variant_count = 0
            #
            with open(sample_count_tsv, 'r') as fin:
                #
                for line in fin:
                    variant_sequence = line.strip().split('\t')[0]
                    replicate_name = line.strip().split('\t')[1]
                    #
                    if variant.sequence == variant_sequence and replicate.name == replicate_name:
                        count_variant += int(line.strip().split('\t')[3])
                    #
                    if variant.sequence in line:
                        variant_count += int(line.strip().split('\t')[3])
                lfn_value = count_variant / variant_count
                #
                if lfn_value < int(cutoff):
                    failed_variants_list.append(variant.variant_id)
        #
        failed_variants[replicate.name] = failed_variants_list
    conn.close()
    return failed_variants

# def lfn_per_replicate_filter1(engine, replicate_model, variant_model, marker_name, sample_count_tsv):
#     # Selecting lines of the replicate table
#     replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
#         replicate_model.marker_name == marker_name)
#     replicate_obj = engine.execute(replicate_select)
#     failed_variants = {}
#     for replicate in replicate_obj:
#         failed_variants_list = []
#         sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
#         # Selecting lines of the variant table
#         variant_select = select([variant_model.variant_id, variant_model.sequence]) \
#             .where(variant_model.marker == marker_name)
#         variant_obj = engine.execute(variant_select)
#         for variant in variant_obj:
#             count_variant = 0
#             replicate_count = 0
#             with open(sample_count_tsv, 'r') as fin:
#                 for line in fin:
#                     if sample_replicate in line and variant.sequence in line:
#                         count_variant += int(line.strip().split('\t')[3])
#                     if sample_replicate in line:
#                         replicate_count += int(line.strip().split('\t')[3])
#                 if replicate_count != 0:
#                     lfn_value = count_variant/replicate_count
#                     if lfn_value < 0.001:
#                         failed_variants_list.append(variant.variant_id)
#         failed_variants[sample_replicate] = failed_variants_list
#     return failed_variants
#
#
# def lfn_per_variant_filter2(engine, replicate_model, variant_model, marker_name, sample_count_tsv):
#     # Selecting lines of the variant table
#     variant_select = select([variant_model.variant_id, variant_model.sequence]) \
#         .where(variant_model.marker == marker_name)
#     variant_obj = engine.execute(variant_select)
#     failed_variants = {}
#     for variant in variant_obj:
#         # Selecting lines of the replicate table
#         replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
#             replicate_model.marker_name == marker_name)
#         replicate_obj = engine.execute(replicate_select)
#         failed_variants_list = []
#         for replicate in replicate_obj:
#             sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
#             with open(sample_count_tsv, 'r') as fin:
#                 count_variant = 0
#                 variant_count = 0
#                 for line in fin:
#                     if sample_replicate in line and variant.sequence in line:
#                         count_variant += int(line.strip().split('\t')[3])
#                     if sample_replicate in line:
#                         variant_count += int(line.strip().split('\t')[3])
#                 if variant_count != 0:
#                     lfn_value = count_variant / variant_count
#                     if lfn_value < 0.001:
#                         failed_variants_list.append(variant.variant_id)
#         failed_variants[variant.variant_id] = failed_variants_list
#     return failed_variants