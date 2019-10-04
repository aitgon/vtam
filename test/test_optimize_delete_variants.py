from io import StringIO
import pandas
from unittest import TestCase

from vtam.wrapper.OptimizeUtillities import KnownVariantAnalyzer

class TestSingleton(TestCase):

    def test_1(self):

        variant_known_str = """run_id	marker_id	biosample_id	variant_id	biosample_type	action
1	1	2	3	mock	keep
1	1	2	66	mock	tolerate
1	1	5	3	mock	keep
1	1	5	66	mock	tolerate
1	1	6		negative	delete
1	1	7		negative	delete
1	1	8	7	real	delete"""

        variant_read_count_str = """run_id	marker_id	biosample_id	variant_id	replicate_id	read_count
1	1	2	3	1	1310
1	1	2	3	2	1551
1	1	2	3	3	1412
1	1	5	3	1	1774
1	1	5	3	2	1406
1	1	5	3	3	2260
1	1	2	10	1	2
1	1	2	66	1	26
1	1	2	66	2	61
1	1	2	66	3	49
1	1	5	66	1	69
1	1	5	66	2	53
1	1	5	66	3	85
1	1	6	1257	1	1
1	1	7	943	1	3
1	1	8	7	1	3
1	1	8	7	2	6
1	1	8	7	3	5
1	1	8	8	1	7
1	1	8	8	2	8
1	1	8	8	3	9"""

        result_delete_str = """run_id	marker_id	biosample_id	variant_id
1	1	2	10
1	1	6	1257
1	1	7	943
1	1	8	7"""

        variant_known_strio = StringIO(variant_known_str)
        variant_known_df = pandas.read_csv(variant_known_strio, sep="\t", header=0)

        variant_read_count_strio = StringIO(variant_read_count_str)
        variant_read_count_df = pandas.read_csv(variant_read_count_strio, sep="\t", header=0)

        result_delete_strio = StringIO(result_delete_str)
        result_delete_df = pandas.read_csv(result_delete_strio, sep="\t", header=0)

        variant_biosample_df = variant_read_count_df[['run_id', 'marker_id', 'variant_id', 'biosample_id']]\
            .drop_duplicates(inplace=False)

        ##########################################################
        #
        # Get keep variants, that is variants marked as keep in either mock or real biosamples
        #
        ##########################################################
        known_variant_analyzer = KnownVariantAnalyzer(variant_known_df=variant_known_df, variant_read_count_df=variant_read_count_df)
        variant_keep_df = known_variant_analyzer.get_variant_keep_df()

        ##########################################################
        #
        # Get delete variants, that are not keep in mock samples
        #
        ##########################################################

        #
        variant_delete_df = known_variant_analyzer.get_variant_delete_df()
        #
        self.assertTrue(result_delete_df.equals(variant_delete_df))
