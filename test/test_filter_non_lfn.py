import os
import pandas
from unittest import TestCase

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner
from wopmetabarcoding.wrapper.FilterMinReplicateNumber import f9_delete_min_replicate_number
from wopmetabarcoding.wrapper.FilterNonLFN import FilterNonLFNRunner


class TestFilterNonLFN(TestCase):

    def setUp(self):
        pass

    def test_01_f9_delete_min_repln(self):
        #
        variant_read_count_df = pandas.DataFrame({
            'run_id': [1]*3,
            'marker_id': [1]*3,
            'variant_id': [1]*2 + [2]*1,
            'biosample_id':[1,1,1],
            'replicate_id':[1,2,1],
            'read_count':[
                10,5,249,
                  ],
        })
        min_replicate_number = 2
        #
        df_filter_output = f9_delete_min_replicate_number(variant_read_count_df = variant_read_count_df, min_replicate_number=min_replicate_number)
        #
        self.assertTrue(not df_filter_output.loc[
                            (df_filter_output.run_id == 1)
                            & (df_filter_output.marker_id == 1)
                            & (df_filter_output.variant_id == 1)
                         & (df_filter_output.biosample_id == 1)
                         & (df_filter_output.replicate_id == 1)
                         & (df_filter_output.filter_id == 9),
                                   'filter_delete'].values[0])
        self.assertTrue(df_filter_output.loc[
                            (df_filter_output.run_id == 1)
                            & (df_filter_output.marker_id == 1)
                            & (df_filter_output.variant_id == 2)
                         & (df_filter_output.biosample_id == 1)
                         & (df_filter_output.replicate_id == 1)
                         & (df_filter_output.filter_id == 9),
                                  'filter_delete'].values[0])



    def test_02_f10_pcr_error(self):
        #
        variant_df = pandas.DataFrame({
            'id': list(range(1,5)),
            'sequence': [
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGATGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTCGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTCGCTGTAACTTTATCATTATCAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                         ],
        })
        #
        variant_read_count_df = pandas.DataFrame({
            'run_id': [1]*4,
            'marker_id': [1]*4,
            'variant_id': list(range(1,5)),
            'biosample_id': [1]*4,
            'replicate_id': [1]*4,
            'read_count':[
                650,520,60,2,
                  ],
        })
        #
        from math import floor
        length_min = min(variant_df.sequence.apply(len).tolist())  # length of smallest sequence
        self.assertTrue(length_min, 300)
        identity = floor((length_min - 1) / length_min * 100) / 100  # 0.99
        self.assertTrue(length_min, 0.99)
        #
        ###################################################################
        # 1. Make a fasta file with all variants of the sample or replicate
        ###################################################################
        from wopmetabarcoding.utils.constants import tempdir
        PathFinder.mkdir_p(os.path.join(tempdir, "f10_pcr_error"))
        variant_fasta = os.path.join(tempdir, "f10_pcr_error", '{}.fasta'.format("pcr_error"))
        with open(variant_fasta, 'w') as fout:
            for row in variant_df.itertuples():
                id = row.id
                sequence = row.sequence
                fout.write(">{}\n{}\n".format(id, sequence))
        ###################################################################
        # 3 Detect all pairs of variants with only 1 difference in the sequences and strong difference in abundance (readcounts)
        # 3.1 vsearch
        ###################################################################
        # import pdb; pdb.set_trace()
        sample_tsv = os.path.join(tempdir, '{}.tsv'.format("pcr_error"))
        vsearch_usearch_global_args = {'db': variant_fasta,
                                       'usearch_global': variant_fasta,
                                       'id': str(identity),
                                       'maxrejects': 0,
                                       'maxaccepts': 0,
                                       'userout': sample_tsv,
                                       'userfields': "query+target+alnlen+ids+mism+gaps",
                                       }
        from wopmetabarcoding.utils.VSearch import VSearch1
        vsearch_usearch_global = VSearch1(**vsearch_usearch_global_args)
        vsearch_usearch_global.run()

    def test_03_f10_pcr_error_vsearch_output_processing(self):
        #
        # Input from min_replicate_number
        variant_read_count_df = pandas.DataFrame({
            'run_id': [1]*8,
            'marker_id': [1]*8,
            'variant_id': [1]*2 + [2]*2 + [3]*2 + [4]*2,
            'biosample_id': [1]*8,
            'replicate_id': [1, 2]*4,
            'read_count':[
                350, 300, 300, 220, 60, 0, 2, 0,
                  ],
        })
        #
        # Output of vsearch
        pcr_error_vsearch_output_df = pandas.DataFrame({
            'query' : [3, 3, 3, 3, 2, 2, 2, 2, 4, 4, 4, 4, 1, 1],
            'target' : [3, 1, 4, 2, 2, 1, 3, 4, 4, 3, 1, 2, 1, 2],
            'alnlen' : [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300],
            'ids' : [300, 299, 299, 298, 300, 299, 298, 297, 300, 299, 298, 297, 300, 299],
            'mism' : [0, 1, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1],
            'gaps' : [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        })
        #
        # Aggregate by biosample
        variant_read_count_grouped_df = variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).sum().reset_index()
        #
        pcr_error_vsearch_output_df[
            'sum_mism_gaps'] = pcr_error_vsearch_output_df.mism + pcr_error_vsearch_output_df.gaps
        # mism + gaps = 1
        check_read_count_df = pcr_error_vsearch_output_df.loc[pcr_error_vsearch_output_df.sum_mism_gaps == 1,['query', 'target']]
        import pdb; pdb.set_trace()



