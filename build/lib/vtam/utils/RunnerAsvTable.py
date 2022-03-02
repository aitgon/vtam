from vtam.models.Sample import Sample
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.utils.NameIdConverter import NameIdConverter
from vtam.utils.SequenceClusterer import SequenceClusterer
from vtam.utils.DataframeVariantReadCountLike import DataframeVariantReadCountLike


class RunnerAsvTable(object):

    def __init__(self, variant_read_count_df, engine, sample_list, cluster_identity, known_occurrences_df=None):

        self.variant_read_count_df = variant_read_count_df
        self.engine = engine
        self.sample_list = sample_list
        self.cluster_identity = cluster_identity
        self.known_occurrences_df = known_occurrences_df

    def to_tsv(self, asvtable_path):

        asvtable_df = self.create_asvtable_df()
        asvtable_df.to_csv(asvtable_path, sep="\t", header=True, index=False)

    def create_asvtable_df(self):

        """Merge asvtable information and reorder columns"""

        asvtable_variant_info_df = self.get_asvtable_variants()
        asvtable_samples_df = self.get_asvtable_samples()

        # Merge variant and sample sides of asvtables
        asvtable_df = asvtable_variant_info_df.merge(asvtable_samples_df, on=['run_id', 'marker_id', 'variant_id'])

        ############################################################################################
        #
        # Label keep variants if known_occurrences not None
        #
        ############################################################################################

        if not (self.known_occurrences_df is None):

            sample_name = NameIdConverter(id_name_or_sequence_list=self.known_occurrences_df.sample_id.tolist(), engine=self.engine).to_names(Sample)
            self.known_occurrences_df.sample_id = ['keep_{}'.format(x) for x in sample_name]
            variant_keep_info_df = self.known_occurrences_df.pivot_table(index=['run_id', 'marker_id', 'variant_id'], columns='sample_id', values='action', aggfunc='first', fill_value=0).reset_index()
            variant_keep_info_df.replace(to_replace='keep', value=1, inplace=True)
            asvtable_df = asvtable_df.merge(variant_keep_info_df, on=['run_id', 'marker_id', 'variant_id'], how='left')
            asvtable_df.fillna(0, inplace=True)
            for label in asvtable_df.columns:
                if label.startswith('keep_'):
                    asvtable_df[label] = asvtable_df[label].astype('int')

        ############################################################################################
        #
        # Rename
        #
        ############################################################################################

        asvtable_df.run_id = NameIdConverter(id_name_or_sequence_list=asvtable_df.run_id.tolist(), engine=self.engine) \
            .to_names(Run)
        asvtable_df.marker_id = NameIdConverter(id_name_or_sequence_list=asvtable_df.marker_id.tolist(), engine=self.engine) \
            .to_names(Marker)
        asvtable_df.rename({'run_id': 'run', 'marker_id': 'marker', 'variant_id': 'variant'}, axis=1, inplace=True)

        ############################################################################################
        #
        # Reorder columns
        #
        ############################################################################################

        column_list = asvtable_df.columns.tolist()
        column_list.remove("clusterid")
        column_list.remove("clustersize")
        column_list.remove("chimera_borderline")
        column_list.remove("sequence")
        column_list = column_list + ['clusterid', 'clustersize', 'chimera_borderline', 'sequence']

        column_list.remove("sequence_length")
        column_list.remove("read_count")
        column_list.insert(3, "sequence_length")
        column_list.insert(4, "read_count")

        asvtable_df = asvtable_df[column_list]

        return asvtable_df

    def get_asvtable_variants(self):

        """This function gets variant-related data such sequence, sequence length and chimera
        borderline status

        Example:
    run_id  marker_id  variant_id  read_count                                           sequence  sequence_length  chimera_borderline  clusterid  clustersize
0        1          1          95         563  ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCT...              175               False       1651            2
1        1          1         271        3471  ACTTTATTTTATTTTTGGTGCTTGATCAGGAATAGTAGGAACTTCT...              175               False       3679            2
2        1          1         603         489  CCTTTATCTTGTATTTGGTGCCTGGGCCGGAATGGTAGGGACCGCC...              175               False        603            1
3        1          1         954       34970  CCTTTATTTTATTTTCGGTATCTGATCAGGTCTCGTAGGATCATCA...              175               False        954            1
4        1          1        1309         755  CTTATATTTTATTTTTGGTGCTTGATCAGGGATAGTGGGAACTTCT...              175               False       1309            2
            """

        asvtable_variant_info_df = DataframeVariantReadCountLike(self.variant_read_count_df).get_N_i_df()
        asvtable_variant_info_df.rename({'N_i': 'read_count'}, axis=1, inplace=True)

        asvtable_variant_info_df['sequence'] = NameIdConverter(
            id_name_or_sequence_list=asvtable_variant_info_df.variant_id.tolist(), engine=self.engine)\
            .variant_id_to_sequence()

        asvtable_variant_info_df['sequence_length'] = asvtable_variant_info_df.sequence.apply(len)

        asvtable_variant_info_df['chimera_borderline'] = NameIdConverter(
            id_name_or_sequence_list=asvtable_variant_info_df.variant_id.tolist(), engine=self.engine)\
            .variant_id_is_chimera_borderline()

        #######################################################################
        #
        # Cluster sequences
        #
        #######################################################################

        seq_clusterer_obj = SequenceClusterer(asvtable_variant_info_df, cluster_identity=self.cluster_identity)
        cluster_count_df = seq_clusterer_obj.compute_clusters()

        asvtable_variant_info_df = asvtable_variant_info_df.merge(cluster_count_df, on='variant_id')

        return asvtable_variant_info_df

    def get_asvtable_samples(self):

        """This function gets sample-related data such as run, marker, variants and read_count per sample

        Example:
sample  run_id  marker_id  variant_id  tpos1_run1  tnegtag_run1  14ben01  14ben02
0               1          1          95         563             0        0        0
1               1          1         271        3471             0        0        0
2               1          1         603         489             0        0        0
3               1          1         954       34970             0        0        0
4               1          1        1309         755             0        0        0
5               1          1        2371        9124             0        0        0
6               1          1        2684           0             0        0     5343
7               1          1        4134           0             0        0      990
            """

        asvtable_2nd_df = DataframeVariantReadCountLike(self.variant_read_count_df).get_N_ij_df()
        asvtable_2nd_df.sample_id = NameIdConverter(id_name_or_sequence_list=asvtable_2nd_df.sample_id.tolist(), engine=self.engine)\
            .to_names(Sample)
        asvtable_2nd_df.rename({'sample_id': 'sample'}, axis=1, inplace=True)
        asvtable_2nd_df = asvtable_2nd_df.pivot_table(index=['run_id', 'marker_id', 'variant_id'], columns='sample',
                            values='N_ij', fill_value=0).reset_index()

        ############################################################################################
        #
        # Set order and fill 0-read count samples with Zeros
        #
        ############################################################################################

        asvtable_2nd_2_df = (asvtable_2nd_df.copy())[['run_id', 'marker_id', 'variant_id']]
        for sample_name in self.sample_list:
            if sample_name in asvtable_2nd_df.columns:  # sample with read counts
                asvtable_2nd_2_df[sample_name] = asvtable_2nd_df[sample_name].tolist()
            else:  # sample without read counts
                asvtable_2nd_2_df[sample_name] = [0] * asvtable_2nd_df.shape[0]

        return asvtable_2nd_2_df
