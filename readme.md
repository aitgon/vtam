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
python -m unittest test.test_wopmetabarcoding.TestWopMetabarcoding.test_01_sample_information
python -m unittest test.test_wopmetabarcoding.TestWopMetabarcoding.test_03_sort_reads
~~~

