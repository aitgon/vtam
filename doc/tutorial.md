% Tutorial
% Aitor Gonzalez, Emese Meglecz

# Installation

## Build Conda environment

Create a conda environment with vsearch and blast

~~~
conda create --name vtam python=3.7 -y
~~~

Activate conda

~~~
conda activate vtam
~~~

## Install VTAM and its dependencies

Download latest release of VTAM: https://www.dropbox.com/sh/g1mshtlet0ymyud/AAA4BNPGpJlewTvJ172B3TFHa?dl=0

~~~
tar zxvf vtam-0.0.1.tar.gz
cd vtam-0.0.1
~~~

_TODO_

VTAM depends on these List of dependencies:

* vsearch=2.7.0
* blast=2.9.0
* wopmars=0.0.8 ( https://github.com/aitgon/wopmars )

These dependencies and VTAM can be installed automatically using within the Conda environment with this command:

~~~
make
~~~

# Merge FASTQ files

The first step is to merge the FASTQ files.

## Data

We can use the dataset from our previous publication: [PMID 28776936](https://pubmed.ncbi.nlm.nih.gov/28776936). 
You can download these FASTQ files from the Dryad website [doi:10.5061/dryad.f40v5](https://datadryad.org/stash/dataset/doi:10.5061/dryad.f40v5).
Create a *fastq* directory:

~~~
mkdir -p fastq
~~~

Copy the FASTQ file inside and verify the directory:

~~~
$ ls fastq/*.fastq
fastq/MFZR1_S4_L001_R1_001.fastq  fastq/MFZR2_S5_L001_R2_001.fastq  fastq/ZFZR1_S1_L001_R1_001.fastq  fastq/ZFZR2_S2_L001_R2_001.fastq
fastq/MFZR1_S4_L001_R2_001.fastq  fastq/MFZR3_S6_L001_R1_001.fastq  fastq/ZFZR1_S1_L001_R2_001.fastq  fastq/ZFZR3_S3_L001_R1_001.fastq
fastq/MFZR2_S5_L001_R1_001.fastq  fastq/MFZR3_S6_L001_R2_001.fastq  fastq/ZFZR2_S2_L001_R1_001.fastq  fastq/ZFZR3_S3_L001_R2_001.fastq
~~~

## Define FASTQ file sample information

Create a TSV (tab-separated file), with a header and 10 columns with all the information per FASTQ file pair.

These columns are needed

- Tag_fwd
- Primer_fwd
- Tag_rev
- Primer_rev
- Marker
- Biosample
- Replicate
- Run
- Fastq_fw
- Fastq_rv


There is an example in the "tutorial" folder with these first lines

~~~
$ head -n3 fastqinfo.tsv 
TagPair Forward	Primer Forward	TagPair Reverse	Primer Reverse	Marker name	 Biosample	Replicate	Run	Fastq_fw	Fastq_rv
tcgatcacgatgt	TCCACTAATCACAARGATATTGGTAC	tgtcgatctacagc	WACTAATCAATTWCCAAATCCTCC	MFZR	Tpos1_prerun	1	prerun	MFZR1_S4_L001_R1_001.fastq	MFZR1_S4_L001_R2_001.fastq
tgatcgatgatcag	TCCACTAATCACAARGATATTGGTAC	tgtcgatctacagc	WACTAATCAATTWCCAAATCCTCC	MFZR	Tpos2_prerun	1	prerun	MFZR1_S4_L001_R1_001.fastq	MFZR1_S4_L001_R2_001.fastq
~~~

## Run the VTAM merge command

In addition to *fastq* and *fastqinfo.tsv*, we need

- Output WopMars DB file
- Output TSV file with FASTA file sample information
- Output directory to write the merged FASTA files

~~~
mkdir -p out/fasta
vtam merge --fastqinfo fastqinfo.tsv --fastqdir fastq --fastainfo out/fastainfo.tsv --fastadir out/fasta --log out/vtam.log -v
~~~

Open the *fastainfo.tsv* file and verify its content. A new column should be written with the names of the merged FASTA files.

Verify also the content of the *out/fasta* with the merged FASTA files.

# Demultiplex and trim the reads

There is a single command *sortreads* to trim and demultiplex the reads. This command takes quite long but its progress can be seen in the log file.

~~~
vtam sortreads --fastainfo out/fastainfo.tsv --fastadir out/fasta --outdir out/fasta_trimmed --log out/vtam.log -vv
~~~

# Filter variants and create the ASV tables

This command filter variants and create the ASV tables. 

~~~
vtam filter --fastainfo out/fastainfo.tsv --fastadir out/fasta --db out/db.sqlite --outdir out --log out/vtam.log -v
~~~

The variants that passed all the filters together with read count in the different biosamples are found in the *out/asvtable.tsv*. The variants that were removed by the different filters can be found in the *out/db.sqlite* database that can be opened with the *sqlitebrowser* program.

# Pool Markers

When variants were amplified with different markers, these variants can be pooled around a variant centroid with the following commands.

An input TSV file must be given with the run and marker combinations that must be pooled. Eg, this is the *pool_run_marker.tsv* file:

~~~
run_name	marker_name
prerun	MFZR
prerun	ZFZR
~~~

Then the *pool_markers* subcommand can be used:


~~~
vtam pool_markers --db ${DB} --runmarker pool_run_marker.tsv --output out/pooled_markers.tsv
~~~

# Taxon Assignation

There is the 'taxassign' subcommand that can assign taxa. 

To assign variants to taxa, we need the COI blast DB and the taxonomy information.

Precomputed versions of these files can be downloaded in the following way:

~~~
vtam taxonomy -o out/taxonomy.tsv --precomputed
vtam coi_blast_db --coi_blast_db out/coi_blast_db
~~~

The input file is a TSV file, where the last column are the sequence of the variants. Both the *out/asvtable.tsv* and *pool_run_marker.tsv* can be used for the assignation.

The command to carry out the taxon assignation is:

~~~
vtam taxassign --variants out/pooled_markers.tsv --output out/pooled_markers_taxa.tsv --db out/db.sqlite --taxonomy out/taxonomy.tsv --blastdbdir out/coi_blast_db --blastdbname coi_blast_db --log out/vtam.log
~~~

# Parameter Optimization

To help the user select the parameters, VTAM has an *optimize* subcommand that will compute different values based on positive and negative variants present in the mock, negative and real biosamples. The set of known variants are defined in a TSV file like this: :download:`variant_known.tsv <variant_known.tsv>`

~~~
vtam optimize --fastainfo out/fastainfo.tsv --fastadir out/fasta --variant_known variant_known.tsv --db out/db.sqlite --outdir out --log out/vtam.log -v
~~~

