import os

from Bio import SeqIO

from vtam.utils.PathManager import PathManager
from vtam.utils.VariantDFutils import VariantDFutils
from vtam.utils.VariantReadCountDF import VariantReadCountDF
from vtam.utils.VSearch import VSearch


class FilterChimeraRunner(object):

    def __init__(self, variant_df, variant_read_count_df):
        """Carries out a chimera analysis"""

        self.variant_df = variant_df
        self.variant_read_count_df = variant_read_count_df

    def run(self, tmp_dir):

        this_step_tmp_dir = os.path.join(tmp_dir, os.path.basename(__name__))
        PathManager.instance().mkdir_p(this_step_tmp_dir)

        filter_output_chimera_df = self.variant_read_count_df.copy()
        filter_output_chimera_df['filter_delete'] = False
        #
        filter_output_borderline_df = self.variant_read_count_df.copy()
        filter_output_borderline_df['filter_delete'] = False

        run_marker_biosample_df = self.variant_read_count_df[['run_id', 'marker_id', 'biosample_id']].drop_duplicates()
        for row in run_marker_biosample_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id

            variant_read_count_df = self.variant_read_count_df.loc[(self.variant_read_count_df.run_id == run_id)
                                                                   & (self.variant_read_count_df.marker_id == marker_id)
                                                                   & (self.variant_read_count_df.biosample_id == biosample_id)]

            variant_read_count_df_obj = VariantReadCountDF(variant_read_count_df=variant_read_count_df)
            N_i_df = variant_read_count_df_obj.get_N_i_df()

            variant_df = self.variant_df.merge(N_i_df, left_on='id', right_on='variant_id')
            variant_df.rename(columns = {'N_i': 'size'}, inplace=True)

            ###################################################################
            #
            # Sort variants by abundance and write to fasta_path
            #
            ###################################################################

            variant_df.sort_values(by='size', ascending=False, inplace=True)

            variant_df_utils_obj = VariantDFutils(variant_df)

            uchime_fasta_path = os.path.join(tmp_dir, os.path.basename(__name__), 'run_{}_marker_{}_biosample_{}.fasta'
                                      .format(run_id, marker_id, biosample_id))
            variant_df_utils_obj.to_fasta(fasta_path=uchime_fasta_path, add_column="size")

            ###################################################################
            #
            # Run uchime_denovo
            #
            ###################################################################
            uchime_borderline_fasta_path = os.path.join(this_step_tmp_dir,
                                                     'run_{}_marker_{}_biosample_{}_borderline.fasta'
                                                         .format(run_id, marker_id, biosample_id))
            uchime_nonchimeras_fasta_path = os.path.join(this_step_tmp_dir,
                                                      'run_{}_marker_{}_biosample_id_{}_nonchimeras.fasta'
                                                         .format(run_id, marker_id, biosample_id))
            uchime_chimeras_fasta_path = os.path.join(this_step_tmp_dir,
                                                   'run_{}_marker_{}_biosample_{}_chimeras.fasta'
                                                         .format(run_id, marker_id, biosample_id))

            #
            # Create object and run vsearch
            vsearch_parameters = {'--uchime_denovo': uchime_fasta_path,
                                  '--borderline': uchime_borderline_fasta_path,
                                  '--nonchimeras': uchime_nonchimeras_fasta_path,
                                  '--chimeras': uchime_chimeras_fasta_path,
                                  }
            vsearch_cluster = VSearch(parameters=vsearch_parameters)
            vsearch_cluster.run()


            ###################################################################
            #
            # 4. Delete variant from replicate/sample if chimeras
            #
            ###################################################################

            with open(uchime_chimeras_fasta_path, "r") as handle:
                for chimera_seqrecord in SeqIO.parse(handle, "fasta"):
                    variant_id = int(chimera_seqrecord.id.split(';')[0])
                    filter_output_chimera_df.loc[(filter_output_chimera_df['run_id'] == run_id)
                                         & (filter_output_chimera_df['marker_id'] == marker_id)
                                         & (filter_output_chimera_df['biosample_id'] == biosample_id)
                                         & (filter_output_chimera_df['variant_id'] == variant_id), 'filter_delete'] = True

            with open(uchime_borderline_fasta_path, "r") as handle:
                for chimera_seqrecord in SeqIO.parse(handle, "fasta"):
                    variant_id = int(chimera_seqrecord.id.split(';')[0])
                    filter_output_borderline_df.loc[(filter_output_borderline_df['run_id'] == run_id)
                                         & (filter_output_borderline_df['marker_id'] == marker_id)
                                         & (filter_output_borderline_df['biosample_id'] == biosample_id)
                                         & (filter_output_borderline_df['variant_id'] == variant_id), 'filter_delete'] = True

        return filter_output_chimera_df, filter_output_borderline_df
