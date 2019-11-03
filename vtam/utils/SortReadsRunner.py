import os
import sys

from Bio import SeqIO

from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.ReadTrimmer import ReadTrimmer


class SortReadsRunner(object):
    """Returns sort_reads.tsv file for a given FASTA"""

    def __init__(self, fasta_path, fasta_id, alignement_parameters, fasta_information_df, sort_reads_tsv):

        self.fasta_path = fasta_path
        self.fasta_id = fasta_id
        self.alignement_parameters = alignement_parameters



        # This is the only order allowed of variant_read_count_df columns
        column_list = ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'tag_fwd_sequence',
                            'primer_fwd_sequence', 'tag_rev_sequence', 'primer_rev_sequence', 'fasta_file_name']
        try:
            fasta_information_df[column_list] = fasta_information_df
            assert fasta_information_df.columns.tolist() == \
                   ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'tag_fwd_sequence', 'primer_fwd_sequence',
                    'tag_rev_sequence', 'primer_rev_sequence', 'fasta_file_name']
        except:
            Logger.instance().error(VTAMexception("This fasta_information_df DataFrame must be composed of columns: "
                                                  "'run_id', 'marker_id', 'biosample_id', 'replicate_id', "
                                                  "'tag_fwd_sequence, primer_fwd_sequence', 'tag_rev_sequence', "
                                                  "'primer_rev_sequence', 'fasta_file_name'. The workflow will exit"))
            sys.exit(1)

        fasta_information_df[['tag_fwd_sequence', 'primer_fwd_sequence', 'tag_rev_sequence', 'primer_rev_sequence']] \
            = fasta_information_df[['tag_fwd_sequence', 'primer_fwd_sequence', 'tag_rev_sequence', 'primer_rev_sequence']]\
            .apply(lambda x: x.str.lower(), axis=0)
        self.fasta_information_df = fasta_information_df

        self.this_tempdir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))

        #######################################################################
        #
        # Output
        #
        #######################################################################

        self.sort_reads_tsv = sort_reads_tsv

    def run(self):

        #######################################################################
        #
        # Trim forward
        #
        #######################################################################

        tag_fwd_sequence_list = self.fasta_information_df.tag_fwd_sequence.tolist()
        primer_fwd_sequence_list = self.fasta_information_df.primer_fwd_sequence.tolist()

        this_tempdir_fwd = os.path.join(self.this_tempdir, 'fwd')
        PathManager.instance().mkdir_p(this_tempdir_fwd)

        # Create ReadTrimmer
        reads_trimmer_fwd = ReadTrimmer(reads_fasta_path=self.fasta_path, tag_sequence_list=tag_fwd_sequence_list,
                                            primer_sequence_list=primer_fwd_sequence_list,
                                            align_parameters=self.alignement_parameters, tempdir=this_tempdir_fwd)
        reads_trimmer_fwd.trim_reads()

        #######################################################################
        #
        # Reverse-complement fasta_path read with trimmed 5 prime
        #
        #######################################################################

        reads_reversed_5prime_trimmed_fasta_path = os.path.join(this_tempdir_fwd, "reads_5prime_trimmed_reversed.fasta_path")
        ReadTrimmer.fasta_file_to_reverse_complement(reads_trimmer_fwd.reads_trimmed_fasta_path,
                                                     reads_reversed_5prime_trimmed_fasta_path)

        #######################################################################
        #
        # Trim reverse
        #
        #######################################################################

        tag_rev_sequence_list = self.fasta_information_df.tag_rev_sequence.tolist()
        primer_rev_sequence_list = self.fasta_information_df.primer_rev_sequence.tolist()

        this_tempdir_rev = os.path.join(self.this_tempdir, 'rev')

        # Create ReadTrimmer with 5 prime trimmed reversed reads
        reads_trimmer_rev = ReadTrimmer(reads_fasta_path=reads_reversed_5prime_trimmed_fasta_path, tag_sequence_list=tag_rev_sequence_list,
                                            primer_sequence_list=primer_rev_sequence_list,
                                            align_parameters=self.alignement_parameters, tempdir=this_tempdir_rev)

        reads_trimmer_rev.trim_reads()

        #######################################################################
        #
        # Reverse-complement fasta_path read with trimmed 5 and 3 primes
        #
        #######################################################################

        reads_trimmed_fasta_path = os.path.join(this_tempdir_rev, "reads_trimmed.fasta")
        ReadTrimmer.fasta_file_to_reverse_complement(reads_trimmer_rev.reads_trimmed_fasta_path,
                                                     reads_trimmed_fasta_path)

        self.annotate_reads(reads_trimmed_fasta_path=reads_trimmed_fasta_path, fasta_id=self.fasta_id)


    def annotate_reads(self, reads_trimmed_fasta_path, fasta_id):
        """Will annotate reads based on the FASTA annotation

        :param reads_trimmed_final_fasta_path: Fasta file with header: >read_id;tag_sequence=...;primer_sequence=...;
            tag_seq...
        :param fasta_information: DataFrame with FASTA information in these columns:
            tag_sequence_fwd, tag_sequence_rev, primer_sequence_fwd, primer_sequence_rev, run_id, marker_id, biosample_id
            replicate_id, fasta_path
        :param sort_reads_tsv_path: Path to TSV with final read annotations
        """

        with open(self.sort_reads_tsv, 'w') as fout:
            for record in SeqIO.parse(reads_trimmed_fasta_path, "fasta"):
                read_id = record.id.split(';')[0]
                tag_fwd = record.id.split(';')[1].split('=')[1]
                primer_fwd = record.id.split(';')[2].split('=')[1]
                tag_rev = record.id.split(';')[3].split('=')[1]
                primer_rev = record.id.split(';')[4].split('=')[1]
                read_sequence = record.seq
                fasta_file_name = os.path.basename( self.fasta_path)
                # import pdb; pdb.set_trace()
                annotation = self.fasta_information_df.loc[(self.fasta_information_df.tag_fwd_sequence == tag_fwd)
                                                           & (self.fasta_information_df.tag_rev_sequence == tag_rev)
                                                           & (self.fasta_information_df.primer_fwd_sequence == primer_fwd)
                                                           & (self.fasta_information_df.primer_rev_sequence == primer_rev)
                                                           & (self.fasta_information_df.fasta_file_name == fasta_file_name)].iloc[0]
                run_id = annotation.run_id
                marker_id = annotation.marker_id
                biosample_id = annotation.biosample_id
                replicate_id = annotation.replicate_id
                out_line = "\t".join([str(i) for i in [read_id, run_id, fasta_id, marker_id, biosample_id, replicate_id, read_sequence]]) + "\n"
                fout.write(out_line)

