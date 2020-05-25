% Developer Notes

# Test tutorial section "Run VTAM with snakemake

## One marker analysis

This snakefile will download these data to the specified working directory

- The fastq test dataset
- The non-rendundant coi_blast_db data
- The taxonomy TSV file

Marker MFZR: data

~~~
mkdir -p out && snakemake -p --cores $(grep -c ^processor /proc/cpuinfo) -d out -s tools/snake.tuto.data.yml --config MARKER=mfzr PROJECT=asper1 --until all_one_marker
~~~

Marker MFZR: pipeline

~~~
snakemake -p --resources db=1 --snakefile doc/data/snakefile.yml --cores 4 -d out --configfile out/asper1/user_input/snakeconfig_mfzr.yml --until asvtable_taxa
~~~

~~~
snakemake -p --resources db=1 --snakefile doc/data/snakefile.yml --cores 4 -d out --configfile out/asper1/user_input/snakeconfig_mfzr.yml --until optimize -n
~~~

~~~
snakemake -p --resources db=1 --snakefile doc/data/snakefile.yml --cores 4 -d out --configfile out/asper1/user_input/snakeconfig_mfzr.yml --until asvtable_optimized_taxa -n
~~~

Marker ZFZR: data

~~~
mkdir -p out && snakemake -p --cores $(grep -c ^processor /proc/cpuinfo) -d out -s tools/snake.tuto.data.yml --config MARKER=zfzr PROJECT=asper1 --until all_one_marker
~~~

Marker ZFZR: pipeline

~~~
snakemake -p --resources db=1 --snakefile doc/data/snakefile.yml --cores 4 -d out --configfile out/asper1/user_input/snakeconfig_zfzr.yml --until asvtable_taxa
~~~

~~~
snakemake -p --resources db=1 --snakefile doc/data/snakefile.yml --cores 4 -d out --configfile out/asper1/user_input/snakeconfig_zfzr.yml --until optimize
~~~

~~~
snakemake -p --resources db=1 --snakefile doc/data/snakefile.yml --cores 4 -d out --configfile out/asper1/user_input/snakeconfig_zfzr.yml --until asvtable_optimized_taxa
~~~

Pooled both markers

~~~
vtam pool --db asper1/db.sqlite --runmarker asper1/user_input/pool_run_marker.tsv --output asper1/asvtable_pooled_mfzr_zfzr.tsv --log asper1/vtam.log -v
~~~

## Two marker analysis

Get data

~~~
mkdir -p out && snakemake -p --cores $(grep -c ^processor /proc/cpuinfo) -d out -s tools/snake.tuto.data.yml --config PROJECT=asper2 --until all_two_markers
~~~

Then we can run the commands

~~~
snakemake -p --resources db=1 --snakefile doc/data/snakefile.yml --cores 4 -d out --configfile out/asper2/user_input/snakeconfig.yml --until asvtable_taxa
~~~

~~~
snakemake -p --resources db=1 --snakefile doc/data/snakefile.yml --cores 4 -d out --configfile out/asper2/user_input/snakeconfig.yml --until optimize
~~~

