import tempfile

import os

tempdir = tempfile.mkdtemp()
print(tempdir)

wopmetabarcoding_filter_test_data = os.path.join(os.path.dirname(__file__), "../../../wopmetabarcoding_filter_test_data")