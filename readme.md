# Installation and use

~~~
conda create --name wopmetabarcoding python=3.6
source activate wopmetabarcoding
pip install -e . # development mode
~~~

# Run

~~~
wopmars -w Wopfile.yml -D sqlite:///db.sqlite -p -v  -F
~~~

# Tests

~~~
python -m unittest discover
python -m unittest test.test_wopmetabarcoding.TestWopMetabarcoding.test_02sample_information
python -m unittest test.test_wopmetabarcoding.TestWopMetabarcoding.test_02sample_information_error
python -m unittest test.test_wopmetabarcoding.TestWopMetabarcoding.test_03sort_reads
python -m unittest test.test_wopmetabarcoding.TestWopMetabarcoding.test_04filter_store_index_below_lfn1_per_replicate
python -m unittest test.test_wopmetabarcoding.TestWopMetabarcoding.test_04filter_store_index_below_lfn2_per_variant
~~~

