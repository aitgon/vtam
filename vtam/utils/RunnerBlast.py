import inspect
import os
import pathlib
import sys

import pandas
from vtam.utils.VTAMexception import VTAMexception

from vtam.utils.PathManager import PathManager

from vtam.utils.Logger import Logger
from Bio.Blast.Applications import NcbiblastnCommandline


class RunnerBlast(object):
    """Runs Blast. Used by Taxassign"""

    def __init__(self, variant_fasta, blast_db_dir, blast_db_name, num_threads,
            qcov_hsp_perc):

        self.variant_fasta = variant_fasta
        self.blast_db_dir = blast_db_dir
        self.blast_db_name = blast_db_name
        # self.ltg_rule_threshold = ltg_rule_threshold
        # self.include_prop = include_prop
        # self.min_number_of_taxa = min_number_of_taxa
        self.num_threads = num_threads
        self.qcov_hsp_perc = qcov_hsp_perc

        self.this_temp_dir = os.path.join(PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(self.this_temp_dir).mkdir(exist_ok=True, parents=True)

    def run_local_blast(self):

        """Runs a local blast and returns the path to the output TSV file"""

        #######################################################################
        #
        # 3 Run local blast
        #
        #######################################################################

        Logger.instance().debug("file: {}; line: {}; Running local blast with FASTA input {}".format(__file__, inspect.currentframe().f_lineno, self.variant_fasta))

        # Run and read local blast result
        blast_output_tsv = os.path.join(self.this_temp_dir, 'blast_output.tsv')
        # blast_output_tsv = "/home/gonzalez/tmp/blast/blast_output.tsv" # uncomment for testing
        # get blast db dir and filename prefix from NHR file
        os.environ['BLASTDB'] = self.blast_db_dir

        blastn_cline = NcbiblastnCommandline(
            query=self.variant_fasta,
            db=self.blast_db_name,
            evalue=1e-5,
            outfmt='"6 qseqid sacc pident evalue qcovhsp staxids"',
            dust='yes',
            qcov_hsp_perc=self.qcov_hsp_perc,
            num_threads=self.num_threads,
            out=blast_output_tsv)
        Logger.instance().debug(
            "file: {}; line: {}; {}".format(
                __file__,
                inspect.currentframe().f_lineno,
                str(blastn_cline)))
        #
        # Run blast
        stdout, stderr = blastn_cline()
        return blast_output_tsv

    @staticmethod
    def process_blast_result(blast_output_tsv):
        """Reads blast_output_tsv and creates a DF that is compatible to the following taxassign. If this DF is empty, vtam will exit with a warning

        """

        Logger.instance().debug(
            "file: {}; line: {}; Reading Blast output from: {}".format(
                __file__, inspect.currentframe().f_lineno, blast_output_tsv))
        blast_output_df = pandas.read_csv(blast_output_tsv, sep='\t',
                                          header=None, names=[
                'variant_id',
                'target_id',
                'identity',
                'evalue',
                'coverage',
                'target_tax_id'])
        # Remove null target tax ids
        blast_output_df = blast_output_df.loc[~blast_output_df.target_tax_id.isnull()]

        # expand multiple target_tax_ids
        # first convert as string
        blast_output_df.target_tax_id = blast_output_df.target_tax_id.astype('str')
        # split by ';' to keep just one target_tax_id and reassign in DF
        blast_output_df.target_tax_id = blast_output_df.target_tax_id.str.split(pat=';', n=1,
                                                                                expand=True)
        # Convert back to numeric/int
        blast_output_df.target_tax_id = blast_output_df.target_tax_id.astype('float').astype('int')
        # blast_output_df = (pandas.concat([
        #     blast_output_df, blast_output_df.target_tax_id.str.split(pat=';', n=1, expand=True)],
        #         axis=1))
        # blast_output_df.drop(['target_tax_id'], axis=1, inplace=true)
        # Blast output extract
        """   variant_id  target_id  identity        evalue  coverage  target_tax_id
0           2  MF7836761    99.429  1.620000e-86       100        1469487
1           2  MF7836761    99.429  1.620000e-86       100         189839
2           2  KY2618191    98.857  7.520000e-85       100         189839
3           2  MF7834791    98.857  7.520000e-85       100         189839
4           2  KU9559321    98.857  7.520000e-85       100         189839
"""

        if blast_output_df.shape[0] == 0:
            Logger.instance().warning(
                VTAMexception("Blast did not find any target. "
                              "VTAM will stop here."))
            sys.exit(0)
        return blast_output_df