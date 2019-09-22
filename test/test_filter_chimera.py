
from unittest import TestCase
from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.VSearch import Vsearch3
from Bio import SeqIO
import os
from wopmetabarcoding.utils.utilities import create_step_tmp_dir, tempdir
from wopmetabarcoding.wrapper.FilterChimera import f11_filter_chimera
import pandas

class TestChimera(TestCase):

    def setUp(self):
        # Input from min_replicate_number
        # Variants 1 and 2 are ok but 3-5 are chimeras
        self.variant_df = pandas.DataFrame({
            'id': list(range(1,7)),
            'sequence': [
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'ACTAACCAGGTGATTTGAAGTAAATTAGTTGAGGATTTAGCCGCGCTATCCGGTAATCTCCAAATTAAAACATACCGTTCCATGAGGGCTAGAATTACTTACCGGCCTTCACCATGCCTGCGCTATACGCGCCCACTCTCCCGTTTATCCGTCCAAGCGGATGCAATGCGATCCTCCGCTAAGATATTCTTACGTGTAACGTAGCTATGTATTTTACAGAGCTGGCGTACGCGTTGAACACTTCACAGATGATAGGGATTCGGGTAAAGAGCGTGTTATTGGGGACTTACACAGGCGTAG',
                'CCTGGGTGAGCTCGAGACTCGGGGTGACAGCTCTTCATACATAGAGCGGGGGCGTCGAACGGTCGTGAAAGTCATAGTACCCCGGGTACCAACTTACTGAGGATATTGCTTGAAGCTGTACCGTTTTAGGGGGGGAACGCTGAAGATCTCTTCTTCTCATGACTGAACTCGCGAGGGTCGTGATGTCGGTTCCTTCAAAGGTTAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTAGCGGCTGGCGTCTTGAATCCTCGGTCCCCCTTGTCTTTCCAGATTAATCCATTTCC',
                'AACTATGTACACAAATTTTAGTATATTGGCAGGGATAGTAGGAACTTTACTATCGTTAGTTATCAGAATGGAATTATCAACAGGAAACATGTTAGATGGAGACGGTCAACAATATAACGTAATCGTAACCGCACATGGATTAATAATGATATTCTTCGTGGTTATGCCGGCAATGTTAGGAGGATTTGCAAACTGGTTCATACCAATAATGGTAGGATCACCAGATGTAGCTTTTCCAAGATTAAACAACATTAGCTTATGGTTAATATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTATTCTTCGTGGTTATGCCGGCAATGTTAGGAGGATTTGCAAACTGGTTCATACCAATAATGGTAGGATCACCAGATGTAGCTTTTCCAAGATTAAACAACATTAGCTTATGGTTAATATTATTGCCCCCTAGTTTATTATTATTAGTTGG',

                         ],
        })
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 12,
            'marker_id': [1]*12,
            'variant_id': [6]*2 + [1]*2 + [2]*2 + [3]*2 + [4]*2 + [5]*2,
            'biosample_id': [1] * 12,
            'replicate_id': [1, 2] * 6,
            'read_count':[
                25, 25, 350, 360, 335, 325, 350, 350, 325, 325, 35, 25
                  ],
        })

        self.this_step_tmp_dir = create_step_tmp_dir(__file__)


    def test_02_f11_chimera(self):
        this_filter_name = 'FilterRenkonen'
        filter_output_df = self.variant_read_count_df.copy()
        filter_output_df['filter_id'] = 11
        filter_output_df['filter_delete'] = 0
        #
        df_grouped_variant_read_count = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).sum().reset_index()
        df_grouped_variant_read_count = df_grouped_variant_read_count[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'read_count']] # keep columns
        #
        df_grouped_run_marker_biosample = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'biosample_id']).sum().reset_index()
        df_grouped_run_marker_biosample = df_grouped_run_marker_biosample[['run_id', 'marker_id', 'biosample_id']] # keep columns
        ###################################################################
        #
        # 1. Make a fasta file with all variants for each run, marker, biosample
        #
        ###################################################################
        for row_run_marker_biosample in df_grouped_run_marker_biosample.iterrows():
            run_id = row_run_marker_biosample[1].run_id
            marker_id = row_run_marker_biosample[1].marker_id
            biosample_id = row_run_marker_biosample[1].biosample_id
            #
            df_grouped_biosample_id = df_grouped_variant_read_count.loc[df_grouped_variant_read_count.biosample_id==biosample_id, ['variant_id', 'read_count']]

            ###################################################################
            #
            # 2. Sort variants by abundance and write to fasta
            #
            ###################################################################
            df_grouped_biosample_id.sort_values(by='read_count', ascending=False, inplace=True)
            #
            chimera1_fasta = os.path.join(self.this_step_tmp_dir, 'chimera1_biosample_id_{}.fasta'.format(biosample_id))
            #
            # Prepare 1 fasta file
            # PathFinder.mkdir_p(os.path.join(self.tempdir, this_filter_name))
            with open(chimera1_fasta, 'w') as fasta_fout:
                for variant_id in df_grouped_biosample_id.variant_id.unique():
                    variant_sequence = self.variant_df.loc[self.variant_df.id == variant_id,'sequence'].values[0]
                    read_count = df_grouped_biosample_id.loc[df_grouped_variant_read_count.variant_id == variant_id, 'read_count'].values[0]
                    fasta_fout.write('>{};size={}\n{}\n'.format(variant_id, read_count, variant_sequence))

            ###################################################################
            #
            # 3. Run uchime_denovo
            #
            ###################################################################
            chimear2_borderline_fasta = os.path.join(self.this_step_tmp_dir, 'chimear2_biosample_id_{}_borderline.fasta'.format(biosample_id))
            chimear2_nonchimeras_fasta = os.path.join(self.this_step_tmp_dir, 'chimear2_biosample_id_{}_nonchimeras.fasta'.format(biosample_id))
            chimear2_chimeras_fasta = os.path.join(self.this_step_tmp_dir, 'chimear2_biosample_id_{}_chimeras.fasta'.format(biosample_id))
            #
            vsearch_chimera_args = {
                "uchime_denovo": chimera1_fasta,
                "borderline": chimear2_borderline_fasta,
                "nonchimeras": chimear2_nonchimeras_fasta,
                "chimeras": chimear2_chimeras_fasta
            }
            vsearch_chimera = Vsearch3(**vsearch_chimera_args)
            vsearch_chimera.run()

            ###################################################################
            #
            # 4. Delete variant from replicate/sample if chimeras
            #
            ###################################################################
            with open(chimear2_chimeras_fasta) as fasta_fin:
                chimera_seqio = SeqIO.parse(chimear2_chimeras_fasta, 'fasta')
                for chimera_seqrecord in chimera_seqio:
                    variant_id = int(chimera_seqrecord.id.split(';')[0])
                    filter_output_df.loc[(filter_output_df['run_id'] == run_id)
                                         & (filter_output_df['marker_id'] == marker_id)
                                         & (filter_output_df['biosample_id'] == biosample_id)
                                         & (filter_output_df['variant_id'] == variant_id), 'filter_delete'] = 1
        #
        self.assertTrue(filter_output_df.loc[(filter_output_df.run_id == 1)
                                                                         & (filter_output_df.marker_id == 1)
                                                                         & (filter_output_df.variant_id == 6)
                                                                         & (filter_output_df.biosample_id == 1)
                                                                         & (filter_output_df.replicate_id == 1),
                                                                         # & (filter_output_df.filter_id == 11),
                                                                        'filter_delete'].values[0])
        #
        self.assertTrue(not filter_output_df.loc[(filter_output_df.run_id == 1)
                                                                         & (filter_output_df.marker_id == 1)
                                                                         & (filter_output_df.variant_id == 1)
                                                                         & (filter_output_df.biosample_id == 1)
                                                                         & (filter_output_df.replicate_id == 1),
                                                                         # & (filter_output_df.filter_id == 11),
                                                                        'filter_delete'].values[0])

        #
        df_output,df_output_borderline = f11_filter_chimera(self.variant_read_count_df, self.variant_df, this_step_tmp_dir=self.this_step_tmp_dir)

        self.assertTrue(df_output.loc[(df_output.run_id == 1)
                                             & (df_output.marker_id == 1)
                                             & (df_output.variant_id == 6)
                                             & (df_output.biosample_id == 1)
                                             & (df_output.replicate_id == 1),
                                             # & (df_output.filter_id == 11),
                                             'filter_delete'].values[0])
        #

        self.assertTrue(not df_output.loc[(df_output.run_id == 1)
                                                 & (df_output.marker_id == 1)
                                                 & (df_output.variant_id == 1)
                                                 & (df_output.biosample_id == 1)
                                                 & (df_output.replicate_id == 1),
                                                 # & (df_output.filter_id == 11),
                                                 'filter_delete'].values[0])

        #


