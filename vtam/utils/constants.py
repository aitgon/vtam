import multiprocessing
import yaml

from vtam import __version__

####################################################################################################
#
# Default VTAM parameters_numerical_default
#
####################################################################################################

params_default_str = """################################################################################
# Parameters of the "merge" command
# These parameters are used by the vsearch --fastq_mergepairs tool that underlies the "vtam merge" command
# For a description of these parameters run_name "vsearch --help"
fastq_ascii: 33
fastq_maxee: 1
fastq_maxmergelen: 500
fastq_maxns: 0
fastq_minlen: 50
fastq_minmergelen: 100
fastq_minovlen: 50
fastq_truncqual: 10

################################################################################
# Parameters of the "sortreads" command
# These parameters correspond to the corresponding parametes by cutadapt that underlies the "vtam sortreads" command
# For a description of these parameters run_name "cutadapt --help"
cutadapt_error_rate: 0.1 # -e in cutadapt
cutadapt_minimum_length: 50 # -m in cutadapt
cutadapt_maximum_length: 500 # -M in cutadapt

################################################################################
# Parameters of the "filter" command
# This parameter sets the minimal number of reads of a variant in the whole run_name
global_read_count_cutoff: 2

################################################################################
# Parameters of the "FilterLFN" filter in the "filter" command
# These parameters set the cutoffs for the low frequency noise (LFN) filters
# Occurrence is deleted if N_ijk/N_i < lfn_variant_cutoff
lfn_variant_cutoff: 0.001
# Occurrence is deleted if N_ijk/N_ik < lfn_variant_replicate_cutoff
# This parameter is used if the --lfn_variant_replicate option is set in "vtam filter" or "vtam optimize"
lfn_variant_replicate_cutoff: 0.001
# Occurrence is deleted if N_ijk/N_jk < lfn_ sample lfn_sample_replicate_cutoff
lfn_sample_replicate_cutoff: 0.001
# Occurrence is deleted if N_ijk < lfn_ lfn_read_count_cutoff
lfn_read_count_cutoff: 10

################################################################################
# Parameters of the "FilterMinReplicateNumber" filter in the "filter" command
# Occurrences of a variant in a given sample are retained only if it is present in at least min_replicate_number replicates of the sample
min_replicate_number: 2

################################################################################
# Parameter of the "FilterPCRerror" filter in the "filter" command
# A given variant 1 is eliminated if N_1j/N_2j < pcr_error_var_prop, where variant 2 is identical to variant 1 except a single mismatch
pcr_error_var_prop: 0.1

################################################################################
# Parameter of the "FilterChimera" filter in the "filter" command
# This parameter corresponds to the abskew parameter in the vsearch --uchime3_denovo tool that underlies the vtam FilterChimera
# For a description of this parameter run_name "vsearch --help"
uchime3_denovo_abskew: 16.0

################################################################################
# Parameter of the "FilterRenkonen" filter in the "filter" command
# Quantile renkonen distance to drop more extreme values
# For. a 0.9 value will set the 9th decile of all renkonen distances as cutoff
renkonen_distance_quantile: 0.9

################################################################################
# Parameter of the "FilterIndel" filter in the "filter" command
# If 1, skips this filter for non-coding markers
skip_filter_indel: 0

################################################################################
# Parameter of the "FilterCondonStop" filter in the "filter" command
# If 1, skips this filter for non-coding markers
skip_filter_codon_stop: 0
# Translation table number from NCBI [ link]
# Default NCBI translation table 5: stops: ['TAA', 'UAA', 'TAG', 'UAG']
genetic_code: 5

################################################################################
# Parameter of the "MakeAsvTable" filter in the "filter" command
# Cluster identity value to clusterize sequences
cluster_identity: 0.97

################################################################################
# Parameters of the "taxassign" command
# Blast parameter for the minimum query coverage
qcov_hsp_perc: 80
# The LTG must include include_prop percent of the hits
include_prop: 90
# Minimal number of taxa among the hits to assign LTG when %identity is below ltg_rule_threshold
min_number_of_taxa: 3
ltg_rule_threshold: 97"""

####################################################################################################
#
#  Optimize lfn_readcount_cutoff, lfn_variant_cutoff, lfn_variant_replicate_cutoff, ...
#
####################################################################################################

# 100 110 absolute lfn_nijk_cutoff maximal value, used by the optimized
lfn_nijk_cutoff_global_max = 101
lfn_nijk_cutoff_lst_size = 10
# 0.05 1 absolute maximal lfn_ni_cutoff or lfn_njk_cutoff, used by the optimized
lfn_ni_njk_cutoff_global_max = 0.051
lfn_ni_njk_cutoff_lst_size = 10

def get_params_default_dic():
    params_default_dic = yaml.load(params_default_str, Loader=yaml.SafeLoader)
    params_default_dic['threads'] = multiprocessing.cpu_count()
    return params_default_dic


####################################################################################################
#
#  Header of these information files
#  PairedFastq, MergedFasta, SortedReadFasta
#
####################################################################################################

header_paired_fastq = {'run', 'marker', 'sample', 'replicate'}
header_merged_fasta = {'run', 'marker', 'sample', 'replicate', 'tagfwd', 'primerfwd', 'tagrev',
    'primerrev', 'mergedfasta'}
header_sortedread_fasta = {'run', 'marker', 'sample', 'replicate', 'sortedfasta'}
header_known_occurrences = {'run', 'marker', 'sample', 'mock', 'variant', 'action', 'sequence'}
header_cutoff_specific_variant_replicate = {'run', 'marker', 'variant', 'lfn_variant_replicate_cutoff'}
header_cutoff_specific_variant = {'run', 'marker', 'variant', 'lfn_variant_cutoff'}

####################################################################################################
#
#  Tax_assign parameters_numerical_default
#
####################################################################################################

rank_hierarchy = [
    'no rank',
    'phylum',
    'superclass',
    'class',
    'subclass',
    'infraclass',
    'superorder',
    'order',
    'suborder',
    'infraorder',
    'family',
    'subfamily',
    'genus',
    'subgenus',
    'species',
    'subspecies']
rank_hierarchy_asv_table = [
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species']

####################################################################################################
#
# External data
#
####################################################################################################

fastq_tar_gz_url1 = "https://github.com/aitgon/vtam/releases/download/{}/fastq.tar.gz".format(__version__)
fastq_tar_gz_url2 = "https://github.com/aitgon/vtam/releases/latest/download/fastq.tar.gz"
fastq_tar_gz_url3 = "https://osf.io/9mvzw/download"

sorted_tar_gz_url1 = "https://github.com/aitgon/vtam/releases/download/{}/sorted.tar.gz".format(__version__)
sorted_tar_gz_url2 = "https://github.com/aitgon/vtam/releases/latest/download/sorted.tar.gz"
sorted_tar_gz_url3 = "https://osf.io/dhk5g/download"


taxonomy_tsv_gz_url1 = "https://github.com/aitgon/vtam/releases/download/{}/taxonomy.tsv.gz".format(__version__)
taxonomy_tsv_gz_url2 = "https://github.com/aitgon/vtam/releases/latest/download/taxonomy.tsv.gz"
taxonomy_tsv_gz_url3 = "https://osf.io/uzk87/download"

coi_blast_db_gz_url1 = "https://github.com/aitgon/vtam/releases/download/{}/coi_blast_db.tar.gz".format(__version__)
coi_blast_db_gz_url2 = "https://github.com/aitgon/vtam/releases/latest/download/coi_blast_db.tar.gz"
def get_coi_blast_db_gz_url3(blastdbname):
    if blastdbname == "coi_blast_db_20200420":
        coi_blast_db_gz_url3 = "https://osf.io/kw9ms/download"
    elif blastdbname == "coi_blast_db":
        coi_blast_db_gz_url3 = "https://osf.io/45zfd/download"
    return coi_blast_db_gz_url3

identity_list = [100, 99, 97, 95, 90, 85, 80, 75, 70]

####################################################################################################
#
#  FilterLFNreference
#
####################################################################################################

FilterLFNreference_records = [
    {'filter_id': 2, 'filter_name': 'lfn_per_variant'},
    {'filter_id': 3, 'filter_name': 'lfn_per_varian'
                                    't_replicate'},
    {'filter_id': 4, 'filter_name': 'lfn_per_variant_specific'},
    {'filter_id': 5, 'filter_name': 'lfn_per_variant_replicate_specific'},
    {'filter_id': 6, 'filter_name': 'lfn_sample_replicate'},
    {'filter_id': 7, 'filter_name': 'lfn_read_count'},
    {'filter_id': 8, 'filter_name': 'lfn_did_not_passed_all'},
]
