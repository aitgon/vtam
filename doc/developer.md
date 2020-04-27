% Developer Notes

# Test tutorial section "Run VTAM with snakemake"

## Download data for VTAM

This snakefile will download these data to the specified working directory

- The fastq test dataset
- The non-rendundant coi_blast_db data
- The taxonomy TSV file


Marker MFZR

~~~
mkdir -p out && snakemake -p --cores $(grep -c ^processor /proc/cpuinfo) -d out -s tools/snake.tuto.data.yml --config MARKER=mfzr PROJECT=aspersnake --until all_one_marker
~~~

Marker ZFZR

~~~
mkdir -p out && snakemake -p --cores $(grep -c ^processor /proc/cpuinfo) -d out -s tools/snake.tuto.data.yml --config MARKER=mfzr PROJECT=aspersnake --until one_marker
~~~

Two markers MFZR and ZFZR simultaneously

~~~
mkdir -p out && snakemake -p --cores $(grep -c ^processor /proc/cpuinfo) -d out -s tools/snake.tuto.data.yml PROJECT=aspersnake2 --until all_two_markers
~~~



