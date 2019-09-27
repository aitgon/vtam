import os

##########################################################
#
# Default VTAM parameters_numerical
#
##########################################################

parameters_numerical = {
'fastq_minovlen' : 50,
'fastq_maxmergelen' : 300,
'fastq_minmergelen' : 100,
'fastq_minlen' : 50,
'fastq_maxee' : 1,
'fastq_truncqual' : 10,
'fastq_maxns' : 0,
'threads' : 8,
'fastq_ascii' : 33,
'min_id' : 0.8,
'minseqlength' : 32,
'overhang' : 0,
'lfn_variant_threshold' : 0.001,
'lfn_variant_replicate_threshold' : 0.001,
'lfn_biosample_replicate_threshold' : 0.001,
'lfn_read_count_threshold' : 10,
'min_replicate_number' : 2,
'pcr_error_var_prop' : 0.1,
'renkonen_threshold' : 0.5,
'genetic_table_number' : 5,
'identity_threshold' : 97,
'include_prop' : 90,
'min_number_of_taxa' : 3
}

##########################################################
#
# Tax_assign parameters_numerical
#
##########################################################
rank_hierarchy =['no rank', 'phylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order',
                 'suborder', 'infraorder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']
rank_hierarchy_otu_table =['phylum', 'class', 'order', 'family', 'genus', 'species']


wop_dir = os.path.join("{}/Software/repositories/wopmetabarcodin".format(os.environ['HOME']))

public_data_dir = "http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/"

identity_list = [100, 99, 97, 95, 90, 85, 80, 75, 70]

# # Old tax_assign parameters_numerical
# order = [100.0, 97.0, 95.0, 90.0]
# # order = [100.0, 97.0, 95.0, 90.0, 85.0, 80.0]
# taxonomic_levels = {"family": 5, "order": 4, "genus": 3, "species": 2, "subspecies": 1}
