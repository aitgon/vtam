import tempfile

import os

tempdir = tempfile.mkdtemp()
print(tempdir)

order = [100.0, 97.0, 95.0, 90.0]
# order = [100.0, 97.0, 95.0, 90.0, 85.0, 80.0]

wopmetabarcoding_filter_test_data = os.path.join(os.path.dirname(__file__), "../../../wopmetabarcoding_filter_test_data")
