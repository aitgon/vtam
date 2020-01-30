##########################################################
#
# Default VTAM parameters_numerical
#
##########################################################

parameters_numerical = {
    'fastq_ascii': 33,
    'fastq_maxee': 1,
    'fastq_maxmergelen': 500,
    'fastq_maxns': 0,
    'fastq_minlen': 50,
    'fastq_minmergelen': 100,
    'fastq_minovlen': 50,
    'fastq_truncqual': 10,
    'genetic_table_number': 5,
    'ltg_rule_threshold': 97,
    'include_prop': 90,
    'lfn_biosample_replicate_threshold': 0.001,
    'lfn_read_count_threshold': 10,
    'lfn_variant_replicate_threshold': None,
    'lfn_variant_threshold': 0.001,
    'min_id': 0.8,
    'min_number_of_taxa': 3,
    'min_replicate_number': 2,
    'minseqlength': 32,
    'overhang': 0,
    'pcr_error_var_prop': 0.1,
    'skip_filter_codon_stop': 0,
    'skip_filter_indel': 0,
    'upper_renkonen_tail': 0.1,
    'threads': 8,
}

##########################################################
#
# Tax_assign parameters_numerical
#
##########################################################

rank_hierarchy =['no rank', 'phylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order',
                 'suborder', 'infraorder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']
rank_hierarchy_asv_table =['phylum', 'class', 'order', 'family', 'genus', 'species']

public_data_dir = "http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/"
url_taxonomy_tsv = "http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/taxonomy.tsv"

identity_list = [100, 99, 97, 95, 90, 85, 80, 75, 70]

##########################################################
#
# FilterLFNreference
#
##########################################################

FilterLFNreference_records = [
    {'id': 2, 'name': 'lfn_per_variant'},
    {'id': 3, 'name': 'lfn_per_varian'
                      't_replicate'},
    {'id': 4, 'name': 'lfn_per_variant_specific'},
    {'id': 5, 'name': 'lfn_per_variant_replicate_specific'},
    {'id': 6, 'name': 'lfn_biosample_replicate'},
    {'id': 7, 'name': 'lfn_read_count'},
    {'id': 8, 'name': 'lfn_did_not_passed_all'},
]
