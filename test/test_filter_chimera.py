
import os
import pandas
from unittest import TestCase
from sqlalchemy import select
from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.VSearch import VSearch1, Vsearch2, Vsearch3
import pandas, itertools
from Bio import SeqIO
import re
import os
from wopmetabarcoding.utils.constants import tempdir








class TestChimera(TestCase):

    def setUp(self):
        # Input from min_replicate_number
        self.variant_df = pandas.DataFrame({
            'id': list(range(1,6)),
            'sequence': [
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG',
                'AACTATGTACACAAATTTTAGTATATTGGCAGGGATAGTAGGAACTTTACTATCGTTAGTTATCAGAATGGAATTATCAACAGGAAACATGTTAGATGGAGACGGTCAACAATATAACGTAATCGTAACCGCACATGGATTAATAATGATATTCTTCGTGGTTATGCCGGCAATGTTAGGAGGATTTGCAAACTGGTTCATACCAATAATGGTAGGATCACCAGATGTAGCTTTTCCAAGATTAAACAACATTAGCTTATGGTTAATATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTATTCTTCGTGGTTATGCCGGCAATGTTAGGAGGATTTGCAAACTGGTTCATACCAATAATGGTAGGATCACCAGATGTAGCTTTTCCAAGATTAAACAACATTAGCTTATGGTTAATATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTATTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTATTCTTCGTGGTTATGCCGGCAATGTTAGGAGGATTTGCAAACTGGTTCATACCAATAATGGTAGGATCACCAGATGTAGCTTTTCCAAGATTAAACAACATTAGCTTATGGTTAATATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG'
                         ],
        })
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1]*10,
            'marker_id': [1]*10,
            'variant_id': [1]*2 + [2]*2 + [3]*2 + [4]*2 + [5]*2,
            'biosample_id': [ 1,2,1,2,1,2,1,2,1, 2],
            'replicate_id': [1, 2]*5,
            'read_count':[
                350, 300, 300, 400, 50, 0, 260, 240,  0, 50,
                  ],
        })

        self.tempdir = os.path.join(tempdir, "FilterUtilities", "FilterUtilities",
                                    self.__class__.__name__)
        PathFinder.mkdir_p(self.tempdir)



        #

    def test_02_f11_chimera(self):
        #
        # - Test 1: Given variants and read_count for grouped variant for each biosample_id, prepare fasta repl.fas
        #
        this_filter_name="chimera_filter"
        variant_id_list_of_lists = []
        sample_list = self.variant_read_count_df.biosample_id.unique().tolist()

        for row in self.variant_read_count_df[['biosample_id']].drop_duplicates(['biosample_id']).itertuples():
            biosample_id = row.biosample_id
            subset_name = "{}".format(biosample_id)
            df_subset = self.variant_read_count_df.loc[
                                                     (self.variant_read_count_df.biosample_id == biosample_id)]
            variant_id_list_of_lists.append(df_subset.variant_id.unique().tolist())

        for variant_id_list in variant_id_list_of_lists:

                ###################################################################
                # 1. Make a fasta file with all variants of the sample or replicate
                ###################################################################
                PathFinder.mkdir_p(os.path.join(self.tempdir, this_filter_name))
                subset_fasta = os.path.join(self.tempdir, this_filter_name, '{}.fasta'.format(subset_name))
                subset_sortbysize_fasta = os.path.join(self.tempdir, this_filter_name,
                                                       '{}_sortbysize.fasta'.format(subset_name))

                ###################################################################
                # 2. Sort variants by abundance
                ###################################################################
                with open(subset_fasta, 'w') as fout:
                    variant_sequence_list = []
                    for variant_id in variant_id_list:
                        # with engine.connect() as conn:
                        #     stmt = select([self.variant_df.__table__.c.sequence]).where(
                            df = self.variant_df.loc[(self.variant_df.id == variant_id)]
                            variant_sequence = df.sequence
                            variant_sequence_list.append(variant_sequence)
                            fout.write(">{}\n{}\n".format(variant_id, variant_sequence))
                            del df
                vsearch_sortbysize_args = {"sortbysize": subset_fasta, "output": subset_sortbysize_fasta}
                vsearch_sortbysize = Vsearch2(**vsearch_sortbysize_args)
                vsearch_sortbysize.run()

                ###################################################################
                # 3. Run uchime_denovo
                ###################################################################
                subset_borderline_fasta = os.path.join(self.tempdir, this_filter_name,'{}_borderline.fasta'.format(subset_name))
                subset_nonchimeras_fasta = os.path.join(self.tempdir, this_filter_name,'{}_nonchimeras.fasta'.format(subset_name))
                subset_chimeras_fasta = os.path.join(self.tempdir, this_filter_name,'{}_chimeras.fasta'.format(subset_name))
                vsearch_chimera_args = {
                    "uchime_denovo": subset_sortbysize_fasta,
                    "borderline": subset_borderline_fasta,
                    "nonchimeras": subset_nonchimeras_fasta,
                    "chimeras": subset_chimeras_fasta
                }
                vsearch_chimera = Vsearch3(**vsearch_chimera_args)
                vsearch_chimera.run()

                ###################################################################
                # 4. Delete variant from replicate/sample if chimeras
                ###################################################################
                subset_chimera_seqio = SeqIO.parse(open(subset_chimeras_fasta), 'fasta')
                for subset_chimera in subset_chimera_seqio:
                    self.variant_read_count_df.loc[variant_id, 'passed'] = False
                    self.variant_read_count_df.loc[subset_chimera.id, 'f10_chimera'] = False  # does not pass
        import pdb;
        pdb.set_trace()

