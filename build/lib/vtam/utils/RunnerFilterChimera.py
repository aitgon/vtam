import os
import pathlib

from Bio import SeqIO

from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.DataframeVariant import DataframeVariant
from vtam.utils.DataframeVariantReadCountLike import DataframeVariantReadCountLike
from vtam.utils.RunnerVSearch import RunnerVSearch


class RunnerFilterChimera(object):

    def __init__(self, variant_read_count_df):
        """Carries out a chimera analysis"""

        self.variant_read_count_df = variant_read_count_df

    def get_variant_read_count_delete_df(
            self, variant_df, uchime3_denovo_abskew):

        temp_dir = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(temp_dir).mkdir(exist_ok=True)

        filter_output_chimera_df = self.variant_read_count_df.copy()
        filter_output_chimera_df['filter_delete'] = False
        #
        filter_output_borderline_df = self.variant_read_count_df.copy()
        filter_output_borderline_df['filter_delete'] = False

        run_marker_sample_df = self.variant_read_count_df[[
            'run_id', 'marker_id', 'sample_id']].drop_duplicates(inplace=False)
        for row in run_marker_sample_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            sample_id = row.sample_id

            variant_read_count_df = self.variant_read_count_df.loc[(self.variant_read_count_df.run_id == run_id) & (
                self.variant_read_count_df.marker_id == marker_id) & (self.variant_read_count_df.sample_id == sample_id)]

            variant_read_count_df_obj = DataframeVariantReadCountLike(
                variant_read_count_df=variant_read_count_df)
            N_i_df = variant_read_count_df_obj.get_N_i_df()

            variant_size_df = variant_df.merge(
                N_i_df, left_index=True, right_on='variant_id')
            variant_size_df = variant_size_df[[
                'variant_id', 'sequence', 'N_i']]
            variant_size_df.rename(columns={'N_i': 'size'}, inplace=True)
            variant_size_df.set_index('variant_id', inplace=True)

            ###################################################################
            #
            # Sort variants by abundance and write to fasta_path
            #
            ###################################################################

            variant_size_df.sort_values(
                by='size', ascending=False, inplace=True)

            variant_df_utils_obj = DataframeVariant(variant_size_df)

            uchime_fasta_path = os.path.join(
                temp_dir, 'run_{}_marker_{}_sample_{}.fasta' .format(
                    run_id, marker_id, sample_id))
            variant_df_utils_obj.to_fasta(
                fasta_path=uchime_fasta_path, add_column="size")

            ###################################################################
            #
            # Run uchime_denovo
            #
            ###################################################################

            uchime_borderline_fasta_path = os.path.join(
                temp_dir, 'run_{}_marker_{}_sample_{}_borderline.fasta' .format(
                    run_id, marker_id, sample_id))
            uchime_nonchimeras_fasta_path = os.path.join(
                temp_dir, 'run_{}_marker_{}_sample_id_{}_nonchimeras.fasta' .format(
                    run_id, marker_id, sample_id))
            uchime_chimeras_fasta_path = os.path.join(
                temp_dir, 'run_{}_marker_{}_sample_{}_chimeras.fasta' .format(
                    run_id, marker_id, sample_id))

            #
            # Create object and run_name vsearch
            vsearch_parameters = {'uchime3_denovo': uchime_fasta_path,
                                  'borderline': uchime_borderline_fasta_path,
                                  'nonchimeras': uchime_nonchimeras_fasta_path,
                                  'chimeras': uchime_chimeras_fasta_path,
                                  'abskew': uchime3_denovo_abskew,
                                  }
            vsearch_cluster = RunnerVSearch(parameters=vsearch_parameters)
            vsearch_cluster.run()

            ###################################################################
            #
            # 4. Delete variant from replicate/sample if chimeras
            #
            ###################################################################

            Logger.instance().debug(
                "Vsearch uchime chimera tsv_path: {}".format(uchime_chimeras_fasta_path))
            with open(uchime_chimeras_fasta_path, "r") as handle:
                for chimera_seqrecord in SeqIO.parse(handle, "fasta"):
                    variant_id = int(chimera_seqrecord.id.split(';')[0])
                    filter_output_chimera_df.loc[(filter_output_chimera_df['run_id'] == run_id)
                                                 & (filter_output_chimera_df['marker_id'] == marker_id)
                                                 & (filter_output_chimera_df['sample_id'] == sample_id)
                                                 & (filter_output_chimera_df['variant_id'] == variant_id), 'filter_delete'] = True

            Logger.instance().debug("Vsearch uchime chimera borderline tsv_path: {}".format(
                uchime_borderline_fasta_path))
            with open(uchime_borderline_fasta_path, "r") as handle:
                for chimera_seqrecord in SeqIO.parse(handle, "fasta"):
                    variant_id = int(chimera_seqrecord.id.split(';')[0])
                    filter_output_borderline_df.loc[(filter_output_borderline_df['run_id'] == run_id)
                                                    & (filter_output_borderline_df['marker_id'] == marker_id)
                                                    & (filter_output_borderline_df['sample_id'] == sample_id)
                                                    & (filter_output_borderline_df['variant_id'] == variant_id), 'filter_delete'] = True

        return filter_output_chimera_df, filter_output_borderline_df
