import tempfile

import os

tempdir = tempfile.mkdtemp()

order = [100.0, 97.0, 95.0, 90.0]
# order = [100.0, 97.0, 95.0, 90.0, 85.0, 80.0]
taxonomic_levels = {"family": 5, "order": 4, "genus": 3, "species": 2, "subspecies": 1}


wopmetabarcoding_filter_test_data = os.path.join(os.path.dirname(__file__), "../../../wopmetabarcoding_filter_test_data")
