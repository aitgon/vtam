**Changes In Version In Version 0.1.18 (Sep 23, 2020)**

- BUG bugs fixed

**Changes In Version In Version 0.1.17 (Sep 19, 2020)**

- DOC Renamed and improved command-line interface help
- ENH Added a 'vtam example' command to generate a file tree for a quick start

**Changes In Version In Version 0.1.16 (Sep 12, 2020)**

- DOC updated docs and doc files

**Changes In Version In Version 0.1.15 (Sep 10, 2020)**

- BUG Compatible biopython >= 1.78

**Changes In Version In Version 0.1.14 (Sep 8, 2020)**

- ENH Cluster sequences in ASV table
- ENH Label keep occurrences in ASV table
- RFR Renamed biosample to sample

**Changes In Version In Version 0.1.13 (Jul 12, 2020)**

- BUG fixed read fasta from sorted reads
- ENH partially windows compatible

**Changes In Version In Version 0.1.12 (Jul 3, 2020)**

- BUG fixed FilterPCRerror runned with all variants

**Changes In Version In Version 0.1.11 (Jun 22, 2020)**

- ENH verification of --cutoff_specific and --lfn_variant_replicate arguments

**Changes In Version In Version 0.1.10 (Jun 18, 2020)**

- RFR Refactored optimize lfn variant or variant replicate
- TST New tests optimize lfn variant or variant replicate

**Changes In Version In Version 0.1.9 (Jun 12, 2020)**

- BLD requirements.txt missing in build

**Changes In Version In Version 0.1.8 (Jun 12, 2020)**

- BUG Fixed that vtam must run again if --params file is updated (Mantis issues 0002620, 0002621) 
- ENH vtam coi_blast_db progressbar

**Changes In Version In Version 0.1.7 (Jun 8, 2020)**

- RFR vtam coi_blast_db refactor/added --blastdbdir BLASTDBDIR and --blastdbname BLASTDBNAME arguments

**Changes In Version In Version 0.1.6 (Jun 3, 2020)**

- BUG 0002609 Fixed lfn_variant(_replicate) cutoff specific
- BUG Various bugs

**Changes In Version In Version 0.1.5 (Mai 31, 2020)**

- ENH: issue 0002607. check paramater names in params.yml
- BUG: issue 0002606. When using --params in filter lfn_read_count_cutoff takes the value of lfn_sample_replicate_cutoff

**Changes In Version In Version 0.1.4 (Mai 29, 2020)**

- BUG: Fixed filter

**Changes In Version In Version 0.1.3 (Mai 27, 2020)**

- BUG: Fixed OptimizePCRerror

**Changes In Version In Version 0.1.2 (Mai 25, 2020)**

- Bug FilterMinReplicateNumber fixed

**Changes In Version In Version 0.1.1 (Mai 24, 2020)**

- Added tests
- Refactored code
- Change optimization based on explicit 'keep' and 'delete' variants

**Changes In Version In Version 0.1.0 (April 30, 2020)**

- Bugs fixed
- Added tests
- Created test dataset
- Made faster some parts of the code

**Changes In Version In Version 0.0.1.4 (April 13, 2020)**

- Bugs fixed
- Updated to use uchime3_denovo
- Updated to use wopmars 11

**Changes In Version In Version 0.0.1.3 (April 7, 2020)**

- VariantReadCount: Fixed "insert variants"

**Changes In Version In Version 0.0.1.2 (April 5, 2020)**

- 'sortreads' based on cutadapt
- 'filter' commands output to asvtable file instead to output directory
- new 'global_read_count_threshold' that will stop variants below this parameter to entering the database

**Changes In Version In Version 0.0.1.1 (March 22, 2020)**

- Change subcommand "poolmarkers" to "pool"
- Reorder optimize columns and other minor output improvements
- Fixed FilterLFNreference
- Renkonen filter does not run if only one replicate

**Changes In Version In Version 0.0.1 (March 18, 2020)**

-  First version running until the end without apparent bugs affecting results


