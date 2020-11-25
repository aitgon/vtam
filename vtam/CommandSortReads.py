import multiprocessing
import os
import sys

import pandas
import pathlib
import shlex
import shutil
import subprocess

# Compatible with both pre- and post Biopython 1.78:
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None

from Bio.Seq import Seq
from vtam.utils.Logger import Logger
from vtam.utils.FileParams import FileParams
from vtam.utils.PathManager import PathManager
from vtam.utils.FileSampleInformation import FileSampleInformation


class CommandSortReads(object):
    """Class for the Merge command"""

    @staticmethod
    def main(fastainfo, fastadir, sorteddir, params=None, num_threads=multiprocessing.cpu_count()):

        if sys.platform.startswith('win'):
            num_threads = 1

        ############################################################################################
        #
        # params.yml parameters
        #
        ############################################################################################

        params_dic = FileParams(params).get_params_dic()

        cutadapt_error_rate = params_dic['cutadapt_error_rate']
        cutadapt_minimum_length = params_dic['cutadapt_minimum_length']
        cutadapt_maximum_length = params_dic['cutadapt_maximum_length']

        ############################################################################################
        #
        # Loop over tag and primer pairs to demultiplex and trim reads
        #
        ############################################################################################

        merged_fastainfo_df = FileSampleInformation(fastainfo).read_tsv_into_df()

        pathlib.Path(sorteddir).mkdir(parents=True, exist_ok=True)
        tempdir = PathManager.instance().get_tempdir()

        sorted_read_info_df = pandas.DataFrame()

        for i in range(0, merged_fastainfo_df.shape[0]):
            fasta_info_series = merged_fastainfo_df.iloc[i]

            tag_fwd = fasta_info_series.tagfwd
            tag_rev = fasta_info_series.tagrev
            primer_fwd = fasta_info_series.primerfwd
            primer_rev = fasta_info_series.primerrev
            in_fasta_basename = fasta_info_series.mergedfasta

            Logger.instance().debug("Analysing FASTA file: {}".format(in_fasta_basename))

            fasta_info_df_i = fasta_info_series.to_frame().T
            in_raw_fasta_path = os.path.join(fastadir, in_fasta_basename)

            ########################################################################################
            #
            # Cut adapt tag of forward reads
            # cutadapt --cores=8 --no-indels --error-rate 0 --trimmed-only
            # --front 'tcgatcacgatgt;min_overlap=13...gctgtagatcgaca;min_overlap=14'
            # --output /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_sorted_000.fasta
            # out/control_mfzr/merged/MFZR1_S4_L001_R1_001_merged.fasta
            #
            ########################################################################################

            if generic_dna:  # Biopython <1.78
                tag_rev_rc = str(Seq(tag_rev, generic_dna).reverse_complement())
            else:  # Biopython =>1.78
                tag_rev_rc = str(Seq(tag_rev).reverse_complement())

            out_fasta_basename = os.path.basename(in_raw_fasta_path).replace(
                '.fasta', '_sorted_%03d.fasta' % i)
            out_fasta_path = os.path.join(tempdir, out_fasta_basename)

            cmd_cutadapt_tag_dic = {
                'tag_fwd': tag_fwd,
                'tag_fwd_len': len(tag_fwd),
                'tag_rev_rc': tag_rev_rc,
                'tag_rev_rc_len': len(tag_rev_rc),
                'in_fasta_path': in_raw_fasta_path,
                'out_fasta': out_fasta_path,
                'num_threads': num_threads,
            }

            cmd_cutadapt_tag_str = 'cutadapt --cores={num_threads} --no-indels --error-rate 0 --trimmed-only ' \
                '--front "{tag_fwd};min_overlap={tag_fwd_len}...{tag_rev_rc};min_overlap={tag_rev_rc_len}" ' \
                '--output {out_fasta} {in_fasta_path}'.format(**cmd_cutadapt_tag_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_tag_str))

            if sys.platform.startswith("win"):
                args = cmd_cutadapt_tag_str
            else:
                args = shlex.split(cmd_cutadapt_tag_str)
            run_result = subprocess.run(args=args, capture_output=True, check=True)

            Logger.instance().info(run_result.stdout.decode())
            Logger.instance().info(run_result.stderr.decode())

            ########################################################################################
            #
            # Trim primers from output
            # cutadapt --cores=8 --no-indels --error-rate 0.1 --minimum-length 50 --maximum-length 500 --trimmed-only
            # --front 'TCCACTAATCACAARGATATTGGTAC;min_overlap=26...GGAGGATTTGGWAATTGATTAGTW;min_overlap=24'
            # --output /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_sorted_trimmed_000.fasta
            # /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_sorted_000.fasta
            #
            ########################################################################################

            if generic_dna:  # Biopython <1.78
                primer_rev_rc = str(Seq(primer_rev, generic_dna).reverse_complement())
            else:  # Biopython =>1.78
                primer_rev_rc = str(Seq(primer_rev).reverse_complement())

            in_fasta_path = out_fasta_path
            out_fasta_basename = os.path.basename(in_fasta_path).replace(
                '_sorted_%03d.fasta' % i, '_sorted_trimmed_%03d.fasta' % i)
            out_fasta_path = os.path.join(tempdir, out_fasta_basename)

            cmd_cutadapt_primer_dic = {
                'primer_fwd': primer_fwd,
                'primer_fwd_len': len(primer_fwd),
                'primer_rev_rc': primer_rev_rc,
                'primer_rev_rc_len': len(primer_rev_rc),
                'in_fasta_path': in_fasta_path,
                'out_fasta': out_fasta_path,
                'error_rate': cutadapt_error_rate,
                'read_min_length': cutadapt_minimum_length,
                'read_max_length': cutadapt_maximum_length,
                'num_threads': num_threads,
            }

            cmd_cutadapt_primer_str = 'cutadapt --cores={num_threads} --no-indels --error-rate {error_rate} ' \
                                      '--minimum-length {read_min_length} ' \
                                      '--maximum-length {read_max_length} --trimmed-only  ' \
                                      '--front "{primer_fwd};min_overlap={primer_fwd_len}...{primer_rev_rc};min_overlap={primer_rev_rc_len}" '  \
                '--output {out_fasta} {in_fasta_path}'.format(**cmd_cutadapt_primer_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_primer_str))

            if sys.platform.startswith("win"):
                args = cmd_cutadapt_primer_str
            else:
                args = shlex.split(cmd_cutadapt_primer_str)
            run_result = subprocess.run(args=args, capture_output=True)

            Logger.instance().info(run_result.stdout.decode())
            Logger.instance().info(run_result.stderr.decode())

            ########################################################################################
            #
            # Cut adapt tag of reverse-complement reads
            # cutadapt --cores=8 --no-indels --error-rate 0 --trimmed-only
            # --front 'tgtcgatctacagc;min_overlap=14...acatcgtgatcga;min_overlap=13'
            # --output /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_rc_sorted_000.fasta
            # out/control_mfzr/merged/MFZR1_S4_L001_R1_001_merged.fasta
            #
            ########################################################################################

            if generic_dna:  # Biopython <1.78
                tag_fwd_rc = str(Seq(tag_fwd, generic_dna).reverse_complement())
            else:  # Biopython =>1.78
                tag_fwd_rc = str(Seq(tag_fwd).reverse_complement())

            out_rc_fasta_basename = os.path.basename(in_raw_fasta_path).replace(
                '.fasta', '_rc_sorted_%03d.fasta' % i)
            out_rc_fasta_path = os.path.join(tempdir, out_rc_fasta_basename)

            cmd_cutadapt_tag_dic = {
                'tag_fwd': tag_rev,
                'tag_fwd_len': len(tag_rev),
                'tag_rev_rc': tag_fwd_rc,
                'tag_rev_rc_len': len(tag_fwd_rc),
                'in_fasta_path': in_raw_fasta_path,
                'out_fasta': out_rc_fasta_path,
                'num_threads': num_threads,
            }

            cmd_cutadapt_tag_str = 'cutadapt --cores={num_threads} --no-indels --error-rate 0 --trimmed-only ' \
                '--front "{tag_fwd};min_overlap={tag_fwd_len}...{tag_rev_rc};min_overlap={tag_rev_rc_len}" ' \
                '--output {out_fasta} {in_fasta_path}'.format(**cmd_cutadapt_tag_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_tag_str))

            if sys.platform.startswith("win"):
                args = cmd_cutadapt_tag_str
            else:
                args = shlex.split(cmd_cutadapt_tag_str)
            run_result = subprocess.run(args=args, capture_output=True)

            Logger.instance().info(run_result.stdout.decode())
            Logger.instance().info(run_result.stderr.decode())

            ###################################################################
            #
            # Trim primers from output
            # cutadapt --cores=8 --no-indels --error-rate 0.1 --minimum-length 50 --maximum-length 500 --trimmed-only
            # --front 'WACTAATCAATTWCCAAATCCTCC;min_overlap=24...GTACCAATATCYTTGTGATTAGTGGA;min_overlap=26'
            # --output /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_rc_sorted_trimmed_000.fasta
            # /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_rc_sorted_000.fasta
            #
            ###################################################################

            if generic_dna:  # Biopython <1.78
                primer_fwd_rc = str(Seq(primer_fwd, generic_dna).reverse_complement())
            else:  # Biopython =>1.78
                primer_fwd_rc = str(Seq(primer_fwd).reverse_complement())

            in_fasta_path = out_rc_fasta_path
            out_rc_fasta_basename = os.path.basename(in_fasta_path).replace(
                '_rc_sorted_%03d.fasta' % i, '_rc_sorted_trimmed_%03d.fasta' % i)
            out_rc_fasta_path = os.path.join(tempdir, out_rc_fasta_basename)

            cmd_cutadapt_primer_dic = {
                'primer_fwd': primer_rev,
                'primer_fwd_len': len(primer_rev),
                'primer_rev_rc': primer_fwd_rc,
                'primer_rev_rc_len': len(primer_fwd_rc),
                'in_fasta_path': in_fasta_path,
                'out_fasta': out_rc_fasta_path,
                'error_rate': cutadapt_error_rate,
                'read_min_length': cutadapt_minimum_length,
                'read_max_length': cutadapt_maximum_length,
                'num_threads': num_threads,
            }
            cmd_cutadapt_primer_str = 'cutadapt --cores={num_threads} --no-indels --error-rate {error_rate} ' \
                '--minimum-length {read_min_length} ' \
                '--maximum-length {read_max_length} --trimmed-only  ' \
                '--front "{primer_fwd};min_overlap={primer_fwd_len}...{primer_rev_rc};min_overlap={primer_rev_rc_len}" ' \
                '--output {out_fasta} {in_fasta_path}'.format(**cmd_cutadapt_primer_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_primer_str))

            if sys.platform.startswith("win"):
                args = cmd_cutadapt_primer_str
            else:
                args = shlex.split(cmd_cutadapt_primer_str)
            run_result = subprocess.run(args=args, capture_output=True)

            Logger.instance().info(run_result.stdout.decode())
            Logger.instance().info(run_result.stderr.decode())

            ###################################################################
            #
            # Reverse complement back rc fasta and pool
            #
            ###################################################################

            out_final_fasta_basename = os.path.basename(
                in_raw_fasta_path).replace('.fasta', '_%03d.fasta' % i)
            out_final_fasta_path = os.path.join(
                sorteddir, out_final_fasta_basename)
            shutil.copy(out_fasta_path, out_final_fasta_path)

            Logger.instance().debug("Pooling fwd and rc reads...")
            with open(out_final_fasta_path, 'a') as fout:
                with open(out_rc_fasta_path, 'r') as fin:
                    for line in fin:
                        if not line.startswith('>'):

                            if generic_dna:  # Biopython <1.78
                                fout.write("%s\n" % str(
                                    Seq(line.strip(), generic_dna).reverse_complement()))
                            else:  # Biopython =>1.78
                                fout.write("%s\n" % str(
                                    Seq(line.strip()).reverse_complement()))

                        else:
                            fout.write(line)

            fasta_info_df_i = fasta_info_df_i[[
                'run', 'marker', 'sample', 'replicate']]
            fasta_info_df_i['sortedfasta'] = out_final_fasta_basename
            sorted_read_info_df = pandas.concat(
                [sorted_read_info_df, fasta_info_df_i], axis=0)

        fasta_trimmed_info_tsv = os.path.join(sorteddir, 'sortedinfo.tsv')
        sorted_read_info_df.to_csv(fasta_trimmed_info_tsv, sep="\t", header=True, index=False)
