import os
import re
import sqlite3
from Bio import SeqIO

from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.VSearch import VSearch


class ReadTrimmer(object):
    """Sort reads according to informations in sample_information file
    """

    def __init__(self, reads_fasta_path, tag_sequence_list, primer_sequence_list, align_parameters, tempdir):
        """
        Inits this class

        :param reads_fasta_path: Path to the FASTA file
        :param fasta_information_df: DataFrame with sample information about FASTA file. Example
                fasta_information_df = pandas.DataFrame(
            {
                'tag_fwd': ['tcgatcacgatgt', 'tgatcgatgatcag'],
                'primer_fwd': ['TCCACTAATCACAARGATATTGGTAC', 'TCCACTAATCACAARGATATTGGTAC'],
                'tag_rev': ['tgtcgatctacagc', 'tgtcgatctacagc'],
                'primer_rev': ['WACTAATCAATTWCCAAATCCTCC', 'WACTAATCAATTWCCAAATCCTCC'],
                'run_id': [1, 1],
                'marker_id': [1, 1],
                'biosample_id': [2, 5],
                'replicate_id': [2, 2],
                'reads_fasta_path': [reads_fasta_path, reade_fasta_path],
             }
        """
        self.reads_fasta_path = reads_fasta_path
        self.align_parameters = align_parameters

        self.tag_sequence_list = [x.lower() for x in tag_sequence_list] # LOWER
        self.primer_sequence_list = [x.upper() for x in primer_sequence_list] # UPPER

        self.tempdir = tempdir
        PathManager.instance().mkdir_p(tempdir)
        self.tag_primer_fasta_path = None

        # TSV file with alignement from vsearch
        self.alignements_tsv_path = None

        # TSV file that was filered for high alignements
        self.alignements_with_high_quality_tsv_path = None

        # Fasta file with trimmed reads that were trimmed in 5'
        self.reads_trimmed_fasta_path = None


    def write_tag_primer_fasta(self):
        """
        Creates a primer-tag fasta_path file to run vsearch

        :param reads_fasta_path: Path to read FASTA file
        :param is_forward_strand: Boolean stating whether forward or reverse strand
        :param tag_primer_fasta_path: Output path for the primer tag FASTA
        :return:
        """
        if self.tag_primer_fasta_path is None:
            self.tag_primer_fasta_path = os.path.join(self.tempdir, "tag_primer.fasta_path")

        tag_primer_sequence_list = \
            ["{}{}".format(tag_sequence.lower(), primer_sequence.upper())
             for tag_sequence, primer_sequence in zip(self.tag_sequence_list, self.primer_sequence_list)]

        with open(self.tag_primer_fasta_path, 'w') as fout:
            # tag sequence MUST be lower and primer sequence MUST be upper, because it is important for the
            # alignements_with_high_quality_tsv_path function below
            for tag_primer_sequence in tag_primer_sequence_list:
                fout.write(">{}\n{}\n".format(tag_primer_sequence,tag_primer_sequence))

    def align_tag_primer_to_reads(self):
        """Aligns tag_primer fasta_path to reads

        :param align_parameters: Dictionary with alignement parameters
        """
        if self.tag_primer_fasta_path is None:
            self.write_tag_primer_fasta()

        self.alignements_tsv_path = os.path.join(self.tempdir, "vsearch_align.tsv")

        # Create object and run vsearch
        vsearch_parameters = {'--db': self.tag_primer_fasta_path,
                              '--usearch_global': self.reads_fasta_path,
                              '--id': self.align_parameters["min_id"],
                              '--maxhits': 1,
                              '--maxrejects': 0,
                              '--maxaccepts': 0,
                              '--minseqlength': self.align_parameters["minseqlength"],
                              '--userfields': "query+target+tl+qilo+qihi+tilo+tihi+qrow",
                              '--userout': self.alignements_tsv_path,
                              }
        vsearch_cluster = VSearch(parameters=vsearch_parameters)
        vsearch_cluster.run()


    def keep_alignements_with_high_quality(self):
        """
        Function to select lines of vsearch_align_tsv according to these SRS criteria:
        1. first position of target in the alignment (tilo) is 1
        2. last position of target in the alignment (tihi) equals length of the target (tl)
        3. lower case part of the target (tag) has an exact match (non-case sensitive) to the aligned part of query (qrow)
        4. this match is within -overhang bases in the 5’ end of qrow

        Example of vsearch_output_csv:
        M00842:118:000000000-ABGKE:1:1101:19978:2032|tgatcgatgatcag|AGATATTGGAACWTTATATTTTATTTTTGG	atatcgacagtgWACTAATCAATTWCCAAATCCTCC	36	1	36	1	36	ATATCGACAGTGTACTAATCAATTTCCAAATCCTCC
        M00842:118:000000000-ABGKE:1:1101:19826:1953|acgatccacagtg|AGATATTGGAACWTTATATTTTATTTTTGG	tctcgatgatcagWACTAATCAATTWCCAAATCCTCC	37	1	37	1	37	TCTCGATGATCAGAACTAATCAATTTCCAAATCCTCC

        :param session: Current session of the database
        :param model: Model of the table
        :param vsearch_output_tsv: csv containing the results of the alignment
        :return: void
        """
        # TODO Emese: ask where overhang parameter is used.

        if self.alignements_tsv_path is None:
            self.align_tag_primer_to_reads()

        self.alignements_with_high_quality_tsv_path = os.path.join(self.tempdir, "vsearch_align_high_quality.tsv")

        vsearch_align_tsv_columns = ['fasta_entry_id', 'tag_primer_sequence', 'tl', 'qilo', 'qihi', 'tilo', 'tihi', 'tag_primer_to_read_alignment']
        with open(self.alignements_with_high_quality_tsv_path, 'w') as fout:
            with open(self.alignements_tsv_path, 'r') as fin:
                nb_discarded_reads = 0
                for line in fin:
                    line_values = line.strip().split()
                    line_dic = dict(zip(vsearch_align_tsv_columns, line_values))
                    if line_dic['tilo'] == "1":
                        tag_primer_sequence = line_dic['tag_primer_sequence']
                        tag_sequence = re.sub('[A-Z]', '', tag_primer_sequence)
                        tag_sequence = tag_sequence.upper()
                        vsearch_tl = line_dic['tl']
                        vsearch_tilo = line_dic['tilo']
                        vsearch_tihi = line_dic['tihi']
                        alignement_of_tag_primer_and_read_sequences = line_dic['tag_primer_to_read_alignment']
                        if vsearch_tilo == "1" and vsearch_tihi == vsearch_tl and tag_sequence in alignement_of_tag_primer_and_read_sequences:
                            fout.write(line)
                        else:
                            nb_discarded_reads += 1
        Logger.instance().info("Number of discarded reads: {}".format(nb_discarded_reads))


    def trim_reads(self):

        """
        Takes an alignement file from vsearch and the read fasta_path and trim reads

        :return:
        """


        if self.alignements_with_high_quality_tsv_path is None:
            self.keep_alignements_with_high_quality()

        temp_db_sqlite = os.path.join(self.tempdir, "reads.db")
        self.reads_trimmed_fasta_path = os.path.join(self.tempdir, "read_trimmed.fasta_path")

        try:
            with sqlite3.connect(temp_db_sqlite) as conn:
                conn.execute("DROP TABLE IF EXISTS read_fasta")
                conn.execute("CREATE TABLE  read_fasta (id VARCHAR PRIMARY KEY , seq VARCHAR)")
                for record in SeqIO.parse(self.reads_fasta_path, 'fasta'):
                    conn.execute("INSERT INTO read_fasta (id, seq) VALUES (?, ?)",
                                 (str(record.description.split()[0].split(';')[0]), str(record.seq)))
                conn.commit()
        except UnicodeDecodeError:
            pass
        with open(self.reads_trimmed_fasta_path, 'w') as fout:

            vsearch_align_tsv_columns = ['fasta_entry_id', 'tag_primer_sequence', 'tl', 'qilo', 'qihi', 'tilo', 'tihi',
                                         'tag_primer_to_read_alignment']
            # Go through each alignement
            with open(self.alignements_with_high_quality_tsv_path, 'r') as fin:
                for line in fin:
                    line_values = line.strip().split()
                    line_dic = dict(zip(vsearch_align_tsv_columns, line_values))
                    fasta_entry_id = line_dic['fasta_entry_id']
                    read_id = fasta_entry_id.split(";")[0]
                    vsearch_qihi = line_dic['qihi']
                    tag_primer_sequence = line_dic['tag_primer_sequence']
                    tag_sequence = re.sub('[A-Z]', '', tag_primer_sequence)
                    primer_sequence = re.sub('[a-z]', '', tag_primer_sequence)
                    # line_info = line.strip().split('\t')

                    # Get read sequence
                    with sqlite3.connect(temp_db_sqlite) as conn:
                        c = conn.cursor()
                        c.execute('SELECT seq FROM read_fasta WHERE id=?', (read_id,))
                        read_sequence = c.fetchone()[0]
                        c.close()

                    # Get end position of superposition between tag_primer sequence and read
                    # vsearch_qihi = line_dic['qihi']
                    tag_primer_sequence_superposed_with_read = read_sequence[0:int(vsearch_qihi)]

                    # Trim forward
                    read_sequence_trimmed = read_sequence.replace(tag_primer_sequence_superposed_with_read, "")
                    # line = line.strip() + '\t' + str(reverse_read) + '\n'
                    fasta_entry = ">{};tag_sequence={};primer_sequence={}\n{}\n"\
                        .format(fasta_entry_id, tag_sequence.lower(), primer_sequence.lower(), read_sequence_trimmed)
                    fout.write(fasta_entry)


    @staticmethod
    def fasta_file_to_reverse_complement(in_fasta, out_fasta):
        with open(out_fasta, 'w') as fout:
            for record in SeqIO.parse(in_fasta, "fasta"):
                fasta_entry = ">{}\n{}\n".format(record.id, record.seq.reverse_complement())
                fout.write(fasta_entry)

