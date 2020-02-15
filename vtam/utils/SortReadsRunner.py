import multiprocessing
import os
import pathlib

from Bio import SeqIO

from vtam.utils.PathManager import PathManager
from vtam.utils.ReadTrimmer import ReadTrimmer


class SortReadsRunner(object):
    """Returns sort_reads.tsv file for a given FASTA"""

    def __init__(self, fasta_information_df, fasta_path, alignement_parameters, outdir, num_threads=multiprocessing.cpu_count()):

        self.fasta_path = fasta_path
        self.alignement_parameters = alignement_parameters
        self.outdir = outdir

        self.fasta_information_df = fasta_information_df

        fasta_information_df.TagFwd = fasta_information_df.TagFwd.str.lower()
        fasta_information_df.TagRev = fasta_information_df.TagRev.str.lower()
        fasta_information_df.PrimerFwd = fasta_information_df.PrimerFwd.str.lower()
        fasta_information_df.PrimerRev = fasta_information_df.PrimerRev.str.lower()

        self.tag_fwd_sequence_list = self.fasta_information_df.TagFwd.tolist()
        self.tag_rev_sequence_list = self.fasta_information_df.TagRev.tolist()
        self.primer_fwd_sequence_list = self.fasta_information_df.PrimerFwd.tolist()
        self.primer_rev_sequence_list = self.fasta_information_df.PrimerRev.tolist()

        self.trimmed_fasta_df = self.fasta_information_df.copy()
        self.trimmed_fasta_df['FastaTrimmed'] = None
        for i, row in enumerate(self.trimmed_fasta_df.itertuples()):
            fasta_trimmed_filename = (row.Fasta).replace('.fasta', '_%03d.fasta' % i)
            self.trimmed_fasta_df.loc[row.Index, 'FastaTrimmed'] = fasta_trimmed_filename
            fasta_trimmed_path = os.path.join(outdir, fasta_trimmed_filename)
            if os.path.isfile(fasta_trimmed_path):
                pathlib.Path(fasta_trimmed_path).unlink()

        self.this_tempdir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(self.this_tempdir).mkdir(exist_ok=True)

        self.num_threads = num_threads

    def run(self):

        #######################################################################
        #
        # Trim forward
        #
        #######################################################################

        # tag_fwd_sequence_list = self.fasta_information_df.tag_fwd_sequence.tolist()
        # primer_fwd_sequence_list = self.fasta_information_df.primer_fwd_sequence.tolist()

        this_tempdir_fwd = os.path.join(self.this_tempdir, 'fwd')
        pathlib.Path(this_tempdir_fwd).mkdir(exist_ok=True)

        # Create ReadTrimmer
        reads_trimmer_fwd = ReadTrimmer(reads_fasta_path=self.fasta_path, tag_sequence_list=self.tag_fwd_sequence_list,
                                            primer_sequence_list=self.primer_fwd_sequence_list,
                                            align_parameters=self.alignement_parameters, tempdir=this_tempdir_fwd)
        reads_trimmer_fwd.trim_reads()

        #######################################################################
        #
        # Reverse-complement fasta_path read with trimmed 5 prime
        #
        #######################################################################

        reads_reversed_5prime_trimmed_fasta_path = os.path.join(this_tempdir_fwd, "reads_5prime_trimmed_reversed.fasta")
        ReadTrimmer.fasta_file_to_reverse_complement(reads_trimmer_fwd.reads_trimmed_fasta_path,
                                                     reads_reversed_5prime_trimmed_fasta_path)

        #######################################################################
        #
        # Trim reverse
        #
        #######################################################################

        # tag_rev_sequence_list = self.fasta_information_df.tag_rev_sequence.tolist()
        # primer_rev_sequence_list = self.fasta_information_df.primer_rev_sequence.tolist()

        this_tempdir_rev = os.path.join(self.this_tempdir, 'rev')

        # Create ReadTrimmer with 5 prime trimmed reversed reads
        reads_trimmer_rev = ReadTrimmer(reads_fasta_path=reads_reversed_5prime_trimmed_fasta_path,
                                        tag_sequence_list=self.tag_rev_sequence_list,
                                            primer_sequence_list=self.primer_rev_sequence_list,
                                            align_parameters=self.alignement_parameters, tempdir=this_tempdir_rev,
                                        num_threads=self.num_threads)
        reads_trimmer_rev.trim_reads()

        #######################################################################
        #
        # Reverse-complement fasta_path read with trimmed 5 and 3 primes
        #
        #######################################################################

        reads_trimmed_fasta_path = os.path.join(this_tempdir_rev, "reads_trimmed.fasta")
        ReadTrimmer.fasta_file_to_reverse_complement(reads_trimmer_rev.reads_trimmed_fasta_path,
                                                     reads_trimmed_fasta_path)

        ################################################################################################################
        #
        # Write to trimmed fasta
        #
        ################################################################################################################

        for record in SeqIO.parse(reads_trimmed_fasta_path, "fasta"):  # Loop over each read
            read_id = record.id.split(';')[0]
            tag_fwd = record.id.split(';')[1].split('=')[1]
            primer_fwd = record.id.split(';')[2].split('=')[1]
            tag_rev = record.id.split(';')[3].split('=')[1]
            primer_rev = record.id.split(';')[4].split('=')[1]
            read_sequence = record.seq
            fasta_file_name = os.path.basename( self.fasta_path)
            fasta_trimmed_filename = self.trimmed_fasta_df.loc[
                (self.trimmed_fasta_df.TagFwd == tag_fwd)
                & (self.trimmed_fasta_df.TagRev == tag_rev)
                & (self.trimmed_fasta_df.PrimerFwd == primer_fwd)
                & (self.trimmed_fasta_df.PrimerRev == primer_rev)
                & (self.trimmed_fasta_df.Fasta == fasta_file_name), 'FastaTrimmed'].item()

            fasta_trimmed_path = os.path.join(self.outdir, fasta_trimmed_filename)

            with open(fasta_trimmed_path, 'a') as fout:
                out_line = ">{}\n{}\n".format(read_id, read_sequence)
                fout.write(out_line)

        ################################################################################################################
        #
        # Return trimmed fasta index
        #
        ################################################################################################################

        self.trimmed_fasta_df = self.trimmed_fasta_df[['Run', 'Marker', 'Biosample', 'Replicate', 'FastaTrimmed']]
        self.trimmed_fasta_df.rename({'FastaTrimmed': 'Fasta'}, inplace=True, axis=1)
        # fasta_trimmed_info_tsv = os.path.join(self.outdir, 'fasta_info.tsv')
        # self.trimmed_fasta_df.to_csv(fasta_trimmed_info_tsv, sep="\t", header=True, index=False)
        return self.trimmed_fasta_df

    # def annotate_reads(self, reads_trimmed_fasta_path):
    #
    #     """Will annotate reads based on the FASTA annotation
    #
    #     :param reads_trimmed_final_fasta_path: Fasta file with header: >read_id;tag_sequence=...;primer_sequence=...;
    #         tag_seq...
    #     :param fasta_information: DataFrame with FASTA information in these columns:
    #         tag_sequence_fwd, tag_sequence_rev, primer_sequence_fwd, primer_sequence_rev, run_id, marker_id, biosample_id
    #         replicate, fasta_path
    #     :param sort_reads_tsv_path: Path to TSV with final read annotations
    #     """
    #
    #     # with open(self.sort_reads_tsv, 'w') as fout:
    #     for record in SeqIO.parse(reads_trimmed_fasta_path, "fasta"):
    #         read_id = record.id.split(';')[0]
    #         tag_fwd = record.id.split(';')[1].split('=')[1]
    #         primer_fwd = record.id.split(';')[2].split('=')[1]
    #         tag_rev = record.id.split(';')[3].split('=')[1]
    #         primer_rev = record.id.split(';')[4].split('=')[1]
    #         read_sequence = record.seq
    #         fasta_file_name = os.path.basename( self.fasta_path)
    #         annotation_df = self.fasta_information_df.loc[(self.fasta_information_df.tag_fwd_sequence == tag_fwd) & (self.fasta_information_df.tag_rev_sequence == tag_rev) & (self.fasta_information_df.primer_fwd_sequence == primer_fwd) & (self.fasta_information_df.primer_rev_sequence == primer_rev) & (self.fasta_information_df.fasta_file_name == fasta_file_name)]
    #         import pdb; pdb.set_trace()
    #         if annotation_df.shape[0] > 0: # there is a match
    #             annotation_series = annotation_df.iloc[0] # get row as series
    #             run_id = annotation_series.run_id
    #             marker_id = annotation_series.marker_id
    #             biosample_id = annotation_series.biosample_id
    #             replicate = annotation_series.replicate
    #             out_line = "\t".join([str(i) for i in [read_id, run_id, fasta_id, marker_id, biosample_id, replicate, read_sequence]]) + "\n"
    #             fout.write(out_line)


