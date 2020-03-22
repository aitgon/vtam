% Tutorial
% Aitor Gonzalez, Emese Meglecz

# Install

Create a conda environment and activate conda

~~~
conda create --name vtam python=3.7 -y
conda activate vtam
~~~

- Download the latest release of VTAM
- Untar the file
- Change to the VTAM package

VTAM depends on these dependencies

* vsearch=2.7.0
* blast=2.9.0
* wopmars=0.0.8 ( https://github.com/aitgon/wopmars )

These dependencies can be installed within the Conda environment with this command:

~~~
make
~~~

# Data

In this tutorial, we can use the dataset from our previous publication: [PMID 28776936](https://pubmed.ncbi.nlm.nih.gov/28776936). 
You can download these FASTQ files from the Dryad website [doi:10.5061/dryad.f40v5](https://datadryad.org/stash/dataset/doi:10.5061/dryad.f40v5).
Create a *FASTQ* directory, copy the FASTQ files inside and verify the directory:

~~~
$ mkdir -p fastq
$ ls fastq/*.fastq``
fastq/MFZR1_S4_L001_R1_001.fastq  fastq/MFZR2_S5_L001_R2_001.fastq  fastq/ZFZR1_S1_L001_R1_001.fastq  fastq/ZFZR2_S2_L001_R2_001.fastq
fastq/MFZR1_S4_L001_R2_001.fastq  fastq/MFZR3_S6_L001_R1_001.fastq  fastq/ZFZR1_S1_L001_R2_001.fastq  fastq/ZFZR3_S3_L001_R1_001.fastq
fastq/MFZR2_S5_L001_R1_001.fastq  fastq/MFZR3_S6_L001_R2_001.fastq  fastq/ZFZR2_S2_L001_R1_001.fastq  fastq/ZFZR3_S3_L001_R2_001.fastq
~~~

# Merge FASTQ files

The first step is to merge the FASTQ files.
Create a TSV (tab-separated file), with a header and 10 columns with all the information per FASTQ file pair.

These columns are needed

- TagFwd
- PrimerFwd
- TagRev
- PrimerRev
- Marker
- Biosample
- Replicate
- Run
- FastqFwd
- FastqRev

There is an example in the "tutorial" folder with these first lines

~~~
TagFwd	PrimerFwd	TagRev	PrimerRev	Marker	Biosample	Replicate	Run	FastqFwd	FastqRev
tcgatcacgatgt	TCCACTAATCACAARGATATTGGTAC	tgtcgatctacagc	WACTAATCAATTWCCAAATCCTCC	MFZR	Tpos1_prerun	1	prerun	MFZR1_S4_L001_R1_001.fastq	MFZR1_S4_L001_R2_001.fastq
tgatcgatgatcag	TCCACTAATCACAARGATATTGGTAC	tgtcgatctacagc	WACTAATCAATTWCCAAATCCTCC	MFZR	Tpos2_prerun	1	prerun	MFZR1_S4_L001_R1_001.fastq	MFZR1_S4_L001_R2_001.fastq
~~~

We define an output directory and then we run the *merge* command:

~~~
mkdir -p out/fasta
vtam merge --fastqinfo fastqinfo.tsv --fastqdir fastq --fastainfo out/fastainfo.tsv --fastadir out/fasta --log out/vtam.log -v
~~~

The result is a directory with the merged files in FASTA format and a *fastainfo.tsv* with the FASTA information.
The first lines of the *fastainfo.tsv* look like this:

~~~
TagFwd	PrimerFwd	TagRev	PrimerRev	Marker	Biosample	Replicate	Run	FastqFwd	FastqRev	Fasta
tcgatcacgatgt	TCCACTAATCACAARGATATTGGTAC	tgtcgatctacagc	WACTAATCAATTWCCAAATCCTCC	MFZR	Tpos1_prerun	1	prerun	MFZR1_S4_L001_R1_001.fastq	MFZR1_S4_L001_R2_001.fastq	MFZR1_S4_L001_R1_001_merged.fasta
tgatcgatgatcag	TCCACTAATCACAARGATATTGGTAC	tgtcgatctacagc	WACTAATCAATTWCCAAATCCTCC	MFZR	Tpos2_prerun	1	prerun	MFZR1_S4_L001_R1_001.fastq	MFZR1_S4_L001_R2_001.fastq	MFZR1_S4_L001_R1_001_merged.fasta
~~~

# Demultiplex, sort and trim the reads

There is a single command *sortreads* to trim and demultiplex the reads.
This command takes quite long.

~~~
vtam sortreads --fastainfo out/fastainfo.tsv --fastadir out/fasta --outdir out/sorted --log out/vtam.log
~~~

The sorted reads will be written in individual text files in the *out/sorted* directory.
In the *out/sorted* directory, there is also a *readinfo.tsv* file with the information about each sorted read file that looks like this:

~~~
Run	Marker	Biosample	Replicate	SortedReadFile
prerun	MFZR	Tpos1_prerun	1	MFZR1_S4_L001_R1_001_merged_000.txt
prerun	MFZR	Tpos2_prerun	1	MFZR1_S4_L001_R1_001_merged_001.txt
~~~

# Filter variants and create the ASV table

This command filters variants and create the ASV table. 

~~~
vtam filter --db out/db.sqlite --readinfo out/sorted/readinfo.tsv --readdir out/sorted --outdir out --log out/vtam.log -v
~~~

The variants that passed all the filters together with read count in the different biosamples are found in the *out/asvtable.tsv*. 
The variants that were removed by the different filters can be found in the *out/db.sqlite* database that can be opened with the *sqlitebrowser* program.
These are the first lines of the *out/asvtable.tsv*

~~~
variant_id	marker	run	sequence_length	read_count	Tpos1_prerun	Tpos2_prerun	chimera_borderline	sequence
27	MFZR	prerun	181	564	237	327	False	ACTATACCTTATCTTCGCAGTATTCTCAGGAATGCTAGGAACTGCTTTTAGTGTTCTTATTCGAATGGAACTAACATCTCCAGGTGTACAATACCTACAGGGAAACCACCAACTTTACAATGTAATCATTACAGCTCACGCATTCCTAATGATCTTTTTCATGGTTATGCCAGGACTTGTT
75	MFZR	prerun	175	587	318	269	False	ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC
~~~

# Asign variants of ASV table to taxa

The 'taxassign' command assigns variant sequences in the last column of a TSV file to taxa. The 'taxassign' command need a blast DB and the taxonomy information. Precomputed versions of the blast DB for the COI marker and the taxonomy information can be generated with these commands:

~~~
vtam taxonomy -o out/taxonomy.tsv --precomputed
vtam coi_blast_db --coi_blast_db out/coi_blast_db
~~~

The input file of the 'taxassign' command is a TSV file, where the last column are the sequence of the variants. Both the *out/asvtable.tsv* and *pool_run_marker.tsv* can be used for the assignation.

The command to carry out the taxon assignation with the *asvtable.tsv* is:

~~~
vtam taxassign --db out/db.sqlite --variants out/asvtable.tsv --output out/asvtable_taxa.tsv --taxonomy out/taxonomy.tsv --blastdbdir out/coi_blast_db --blastdbname coi_blast_db --log out/vtam.log -v
~~~

# Pool markers

When variants were amplified with different markers, these variants can be pooled around a centroid variant.

An input TSV file must be given with the run and marker combinations that must be pooled. Eg, this is the *pool_run_marker.tsv* file:

~~~
run_name	marker_name
prerun	MFZR
prerun	ZFZR
~~~

Then the *pool_markers* subcommand can be used:

~~~
vtam pool --db out/db.sqlite --runmarker pool_run_marker.tsv --output out/pooled_markers.tsv
~~~

# Parameter Optimization

To help the user select the parameters, VTAM has an *optimize* subcommand that will compute different values based on positive and negative variants present in the mock, negative and real biosamples. The set of known variants are defined in a TSV file like this: :download:`variant_known.tsv <variant_known.tsv>`

~~~
vtam optimize --db out/db.sqlite --readinfo out/sorted/readinfo.tsv --readdir out/sorted --variant_known  variant_known.tsv --outdir out --log out/vtam.log -v
~~~

