Tutorial
============

To setup a new Wopmetabarcodign workflow, follow these steps:

1. Download Wopfile
2. Check FASTQ location
3. Specify FASTQ pair input and FASTA output information (*sample2fastq.csv*)
4. Specify sample information with FASTA paths (*sample2fasta.csv*)
5. Update *sample2fasta.csv* paths in Wopfile
6. Update workflow parameters in Wopfile (Optional)
7. Run wopmars

Download Wopfile or Workflow definition file
---------------------------------------------------

This is a full :download:`Wopfile example </data/Wopfile.yml>`. Here there are the first lines:

.. code-block:: bash

    rule SampleInformation:
        tool: wopmetabarcoding.wrapper.SampleInformation
        input:
            file:
                csv: data/1_sample_info.csv
    ...

Specify FASTQ pair input and FASTA output information (TO DISCUSS)
---------------------------------------------------------------------

In my case, FASTQ files look this and pairs are given by R1/R2:

.. code-block:: bash

    $ ls /home/gonzalez/Data/2017_meglecz_metabarcoding/data/nuxeo/input/input_fastq
    MFZR1_S4_L001_R1_001.fastq  MFZR2_S5_L001_R2_001.fastq  ZFZR1_S1_L001_R1_001.fastq  ZFZR2_S2_L001_R2_001.fastq
    MFZR1_S4_L001_R2_001.fastq  MFZR3_S6_L001_R1_001.fastq  ZFZR1_S1_L001_R2_001.fastq  ZFZR3_S3_L001_R1_001.fastq
    MFZR2_S5_L001_R1_001.fastq  MFZR3_S6_L001_R2_001.fastq  ZFZR2_S2_L001_R1_001.fastq  ZFZR3_S3_L001_R2_001.fastq

Therefore, my *fastq_pairs.csv* file is composed of three columns with forward and reversed FASTQ and output FASTA:

.. code-block:: bash

    fastq_fw,fastq_rev,fasta
    MFZR1_S4_L001_R1_001.fastq,MFZR1_S4_L001_R2_001.fastq,prerun_MFZR_repl1.fasta
    MFZR2_S5_L001_R1_001.fastq,MFZR2_S5_L001_R2_001.fastq,prerun_MFZR_repl2.fasta
    MFZR3_S6_L001_R1_001.fastq,MFZR3_S6_L001_R2_001.fastq,prerun_MFZR_repl3.fasta
    ZFZR1_S1_L001_R1_001.fastq,ZFZR1_S1_L001_R2_001.fastq,prerun_ZFZR_repl1.fasta
    ZFZR2_S2_L001_R1_001.fastq,ZFZR2_S2_L001_R2_001.fastq,prerun_ZFZR_repl2.fasta
    ZFZR3_S3_L001_R1_001.fastq,ZFZR3_S3_L001_R2_001.fastq,prerun_ZFZR_repl3.fasta


Specify sample information with FASTA paths
---------------------------------------------------------

The sample information file (*sample_information.csv*) is a CSV file with header and 9 columns:

- Forward tag seq
- Forward primer seq
- Reverse tag seq
- Reverse primer seq
- Marker name
- Biosample name
- Replicate
- FASTA file path
- Run name

The absolute FASTA path or relative FASTA path relative to the *sample2fasta.csv* must be correct. For instance, you can have this file hierarchy:

.. code-block:: bash

    .
    ├── fasta
    │   ├── prerun_MFZR_repl2.fasta
    │   ├── prerun_MFZR_repl3.fasta
    │   ├── prerun_ZFZR_repl1.fasta
    │   ├── prerun_ZFZR_repl2.fasta
    │   └── prerun_ZFZR_repl3.fasta
    └── sample2fasta.csv


These are the first lines of the *sample2fasta.csv* file

.. code-block:: bash

    Forward tag seq,Forward primer seq,Reverse tag seq,Reverse primer seq,Marker name,Biosample name,Replicate,FASTA file path,Run name
    cgatcgtcatcacg,TCCACTAATCACAARGATATTGGTAC,cgcgatctgtagag,WACTAATCAATTWCCAAATCCTCC,MFZR,14Mon01,repl2,fasta/prerun_MFZR_repl2.fasta,prerun


Update the *sample2fasta.csv* file path parameter in Wopfile
---------------------------------------------------------------

The simplest is to run Wopmars with the working directory in the same folder as *sample2fasta.csv*. Then you only need to change the *input/file/sample2fasta.csv* in the *SampleInformation* rule as shown here:

.. code-block:: bash

    rule SampleInformation:
        tool: wopmetabarcoding.wrapper.SampleInformation
        input:
            file:
                csv: sample2fasta.csv
        output:
            table:
                File: wopmetabarcoding.model.File
                SampleInformation: wopmetabarcoding.model.SampleInformation
                Marker: wopmetabarcoding.model.Marker
                PrimerPair: wopmetabarcoding.model.PrimerPair
                Biosample: wopmetabarcoding.model.Biosample
                TagPair: wopmetabarcoding.model.TagPair
                Replicate: wopmetabarcoding.model.Replicate


Update workflow parameters in Wopfile (Optional)
--------------------------------------------------------------------------------

In the simplest case, just keep the default parameters.


Run wopmars
--------------------------------------------------------------------------------

A simple wopmars command to test the *SampleInformation* rule is here

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -F -t SampleInformation

This command will store the primers, tags and other informations in the database *db.sqlite*, which can be opened with eg, sqlitebrowser.

In order to run the whole workflow, the files *cutoff_variant.tsv* and *genetic_codes.csv* are necessary and a valid path relative to Wopmars working directory. This path must be defined in the Wopfile:

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -F -t SampleInformation
    
    

