from vtam.models.Biosample import Biosample
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.utils.NameIdConverter import NameIdConverter
from vtam.utils.VariantReadCountLikeDF import VariantReadCountLikeDF


class AsvTableRunner(object):

    def __init__(self, variant_read_count_df, engine, biosample_list):

        self.variant_read_count_df = variant_read_count_df
        self.engine = engine
        self.biosample_list = biosample_list

    def to_tsv(self, asvtable_path):

        asvtable_df = self.get_asvtable_df()
        asvtable_df.to_csv(asvtable_path, sep="\t", header=True, index=False)

    def get_asvtable_df(self):

        asvtable_variants_df = self.get_asvtable_variants()
        asvtable_biosamples_df = self.get_asvtable_biosamples()

        ############################################################################################
        #
        # Merge variant and biosample sides of asvtables
        #
        ############################################################################################

        asvtable_df = asvtable_variants_df.merge(asvtable_biosamples_df, on=['run_id', 'marker_id', 'variant_id'])
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
        column_list.remove("chimera_borderline")
        column_list.remove("sequence")
        column_list = column_list + ['chimera_borderline', 'sequence']

        column_list.remove("sequence_length")
        column_list.remove("read_count")
        column_list.insert(3, "sequence_length")
        column_list.insert(4, "read_count")
        asvtable_df = asvtable_df[column_list]

        return asvtable_df

    def get_asvtable_variants(self):

        asvtable_1st_df = VariantReadCountLikeDF(self.variant_read_count_df).get_N_i_df()
        asvtable_1st_df.rename({'N_i': 'read_count'}, axis=1, inplace=True)

        asvtable_1st_df['sequence'] = NameIdConverter(
            id_name_or_sequence_list=asvtable_1st_df.variant_id.tolist(), engine=self.engine)\
            .variant_id_to_sequence()

        asvtable_1st_df['sequence_length'] = asvtable_1st_df.sequence.apply(len)

        asvtable_1st_df['chimera_borderline'] = NameIdConverter(
            id_name_or_sequence_list=asvtable_1st_df.variant_id.tolist(), engine=self.engine)\
            .variant_id_is_chimera_borderline()

        return asvtable_1st_df

    def get_asvtable_biosamples(self):

        asvtable_2nd_df = VariantReadCountLikeDF(self.variant_read_count_df).get_N_ij_df()
        asvtable_2nd_df.biosample_id = NameIdConverter(id_name_or_sequence_list=asvtable_2nd_df.biosample_id.tolist(), engine=self.engine)\
            .to_names(Biosample)
        asvtable_2nd_df.rename({'biosample_id': 'biosample'}, axis=1, inplace=True)
        asvtable_2nd_df = asvtable_2nd_df.pivot_table(index=['run_id', 'marker_id', 'variant_id'], columns='biosample',
                            values='N_ij', fill_value=0).reset_index()

        ############################################################################################
        #
        # Set order and fill 0-read count biosamples with Zeros
        #
        ############################################################################################

        asvtable_2nd_2_df = (asvtable_2nd_df.copy())[['run_id', 'marker_id', 'variant_id']]
        for biosample_name in self.biosample_list:
            if biosample_name in asvtable_2nd_df.columns:  # biosample with read counts
                asvtable_2nd_2_df[biosample_name] = asvtable_2nd_df[biosample_name].tolist()
            else:  # biosample without read counts
                asvtable_2nd_2_df[biosample_name] = [0] * asvtable_2nd_df.shape[0]
        return asvtable_2nd_2_df
