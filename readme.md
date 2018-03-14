# Installation and use

~~~
conda create --name wopmetabarcoding python=3.6
source activate wopmetabarcoding
pip install -e . # development mode
~~~

# Run

~~~
wopmars -w Wopfile.yml -D sqlite:///output/db.sqlite -p -v  -F 
~~~

