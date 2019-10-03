from io import StringIO
import pandas
from unittest import TestCase


class TestSingleton(TestCase):

    def test_1(self):

        variant_known_str = """
marker	run	biosample_name	biosample_type	variant_id	action
MFZR	prerun	Tpos1_prerun	mock	3	keep
MFZR	prerun	Tpos1_prerun	mock	66	tolerate
MFZR	prerun	Tpos2_prerun	mock	3	keep
MFZR	prerun	Tpos2_prerun	mock	66	tolerate
MFZR	prerun	TnegPai1_prerun	negative		delete
MFZR	prerun	TnegPai2_prerun	negative		delete
MFZR	prerun	P1	real	9	delete
MFZR	prerun	P1	real	7	delete"""

        variant_known_str = """run_id	marker_id	biosample_id	biosample_type	variant_id	action
1	1	2	mock	3	keep
1	1	2	mock	66	tolerate
1	1	5	mock	3	keep
1	1	5	mock	66	tolerate
1	1	6	negative		delete
1	1	7	negative		delete
1	1	7	real	9	delete
1	1	7	real	7	delete"""

        variant_read_count_str = """run_id	marker_id	variant_id	biosample_id	replicate_id	read_count
1	1	3	2	1	1310
1	1	3	2	2	1551
1	1	3	2	3	1412
1	1	3	5	1	1774
1	1	3	5	2	1406
1	1	3	5	3	2260
1	1	66	2	1	26
1	1	66	2	2	61
1	1	66	2	3	49
1	1	66	5	1	69
1	1	66	5	2	53
1	1	66	5	3	85
1	1	1257	6	1	1
1	1	943	7	1	3
1	1	7	44	1	3
1	1	7	44	2	6
1	1	7	44	3	5"""

        variant_known_strio = StringIO(variant_known_str)
        variant_known_df = pandas.read_csv(variant_known_strio, sep="\t", header=0)

        variant_read_count_strio = StringIO(variant_read_count_str)
        variant_read_count_df = pandas.read_csv(variant_read_count_strio, sep="\t", header=0)

        ##########################################################
        #
        # Get keep variants, that is variants marked as keep in either mock or real biosamples
        #
        ##########################################################
        variant_read_count_keep_df = variant_read_count_df.merge(variant_known_df,
                                                                 on=['run_id', 'marker_id', 'biosample_id', 'variant_id'])
        import pdb; pdb.set_trace()
