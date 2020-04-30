import pandas
from unittest import TestCase

from vtam.utils.FilterCodonStopRunner2 import FilterCodonStopRunner2
# from vtam.wrapper.FilterCodonStop import f14_filter_codon_stop
from io import StringIO


class TestFilterCodonStop(TestCase):
    """Emese Meglecz

    codeX refers to the genetic code (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    frameX the reading frame with the smallest number of codon STOP
    stopX refers to the number of codon STOP with using codeX in frameX

    code5 (The Invertebrate Mitochondrial Code (transl_table=5))
    is a code that contains only two codon STOPs (TAA,TAG),
    the most frequently used ones among all genetic codes

    >seq_stop0_code5_stop2_code1_frame2
    AATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATATTTAGTTGG
    >seq_stop0_code5_stop2_code1_frame1
    ATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATATTTAGTTGG
    >seq_stop1_code5_stop3_code1_frame2
    AATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG
    >seq_stop2_code5_stop4_code1_frame2_TAA_TAG
    ATAAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG
    >seq_stop2_code5_stop5_code1_frame1_TAA_TAG_TGA
    TAAGTATTTTGACTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG
    """

    def setUp(self):
        #
        variant_str = """id	sequence
seq_stop0_code5_stop2_code1_frame2	AATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATATTTAGTTGG
seq_stop0_code5_stop2_code1_frame1	ATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATATTTAGTTGG
seq_stop1_code5_stop3_code1_frame2	AATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG
seq_stop2_code5_stop4_code1_frame2_TAA_TAG	ATAAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG
seq_stop2_code5_stop5_code1_frame1_TAA_TAG_TGA	TAAGTATTTTGACTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG"""
        self.variant_df = pandas.read_csv(StringIO(variant_str), sep="\t", header=0)


    def test_count_codon_stop_nb_one_seq(self):

        ################################################################################################################
        # next seq

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=1)
        seqi = self.variant_df.loc[self.variant_df.id == "seq_stop0_code5_stop2_code1_frame2", "sequence"].values[0]
        framei = 2

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=1)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(2 == codon_stop_nb)

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=5)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(0 == codon_stop_nb)

        ################################################################################################################
        # next seq
        seqi = self.variant_df.loc[self.variant_df.id == "seq_stop0_code5_stop2_code1_frame1", "sequence"].values[0]
        framei = 1

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=1)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(2 == codon_stop_nb)

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=5)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(0 == codon_stop_nb)

        ################################################################################################################
        # next seq
        seqi = self.variant_df.loc[self.variant_df.id == "seq_stop1_code5_stop3_code1_frame2", "sequence"].values[0]
        framei = 2

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=1)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(3 == codon_stop_nb)

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=5)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(1 == codon_stop_nb)

        ################################################################################################################
        # next seq
        seqi = self.variant_df.loc[self.variant_df.id == "seq_stop2_code5_stop4_code1_frame2_TAA_TAG", "sequence"].values[0]
        framei = 2

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=1)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(4 == codon_stop_nb)

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=5)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(2 == codon_stop_nb)

        ################################################################################################################
        # next seq
        seqi = self.variant_df.loc[self.variant_df.id == "seq_stop2_code5_stop5_code1_frame1_TAA_TAG_TGA", "sequence"].values[0]
        framei = 1

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=1)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(5 == codon_stop_nb)

        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=5)
        codon_stop_nb = filter_codon_stop_runner_obj.count_codon_stop_nb_one_seq(seq=seqi, frame=framei)
        self.assertTrue(2 == codon_stop_nb)

    def test_annotate_stop_codon_count(self):
        filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=1)
        variant_has_stop_codon_df = filter_codon_stop_runner_obj.annotate_stop_codon_count(self.variant_df)

        variant_stop_codon_count_df_bak_str = """                                               id  has_stop_codon
0              seq_stop0_code5_stop2_code1_frame2               1
1              seq_stop0_code5_stop2_code1_frame1               1
2              seq_stop1_code5_stop3_code1_frame2               1
3      seq_stop2_code5_stop4_code1_frame2_TAA_TAG               1
4  seq_stop2_code5_stop5_code1_frame1_TAA_TAG_TGA               1"""
        variant_stop_codon_count_df_str = variant_has_stop_codon_df[['id', 'has_stop_codon']].to_string()
        self.assertTrue(variant_stop_codon_count_df_str == variant_stop_codon_count_df_bak_str)
