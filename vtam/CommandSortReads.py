from audioop import reverse
import multiprocessing
import os
import sys

import pandas
import pathlib
import shlex
import shutil
import subprocess
import gzip 
import bz2
from functools import partial
import shlex

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
from vtam.utils.FilesInputCutadapt import FilesInputCutadapt

class CommandSortReads(object):
    """Class for the Merge command"""

    @staticmethod
    def main(fastainfo, fastadir, sorteddir, params=None, num_threads=multiprocessing.cpu_count(), 
        no_reverse=False, tag_to_end=False, primer_to_end=False):

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

        merged_fastainfo_df = pandas.read_csv(fastainfo, sep='\t', header=0)
        pathlib.Path(sorteddir).mkdir(parents=True, exist_ok=True)
        tempdir = PathManager.instance().get_tempdir()

        sorted_read_info_df = pandas.DataFrame()
        print(merged_fastainfo_df.loc[:,'mergedfasta'].to_list())
        for mergedfasta in merged_fastainfo_df.loc[:,'mergedfasta'].to_list():

            inputFiles = FilesInputCutadapt(fastainfo, mergedfasta, tag_to_end, primer_to_end)
            
            tagFile = inputFiles.tags_file()
            primerFile = inputFiles.primers_file()

            Logger.instance().debug("Analysing FASTA file: {}".format(mergedfasta))

            #fasta_info_df_i = fasta_info_series.to_frame().T
            in_raw_fasta_path = os.path.join(fastadir, mergedfasta)

            ########################################################################################
            #
            # Cut adapt tag of forward reads
            # cutadapt --cores=8 --no-indels --error-rate 0 --trimmed-only
            # --front 'tcgatcacgatgt;min_overlap=13...gctgtagatcgaca;min_overlap=14'
            # --output /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_sorted_000.fasta
            # out/control_mfzr/merged/MFZR1_S4_L001_R1_001_merged.fasta
            #
            #   cutadapt --cores=0 -e 0 --no-indels --trimmed-only -g tagFile:$tagfile --overlap length -o "tagtrimmed.{name}.fasta" in_raw_fasta_path
            ########################################################################################


            out_fasta_basename = os.path.basename(in_raw_fasta_path).replace(
                '.fasta', f'_sorted_{mergedfasta}')

            out_fasta_path = os.path.join(tempdir, out_fasta_basename)

            cmd_cutadapt_tag_dic = {
                'in_fasta_path': in_raw_fasta_path,
                'out_fasta': out_fasta_path,
                'num_threads': num_threads,
                'tagFile': tagFile,
            }

            cmd_cutadapt_tag_str = 'cutadapt --cores={num_threads} --no-indels --error-rate 0 --trimmed-only ' \
                '-g file:{tagFile} --overlap --output {out_fasta}.{{name}}.fasta {in_fasta_path}' \
                .format(**cmd_cutadapt_tag_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_tag_str))

            if sys.platform.startswith("win"):
                args = cmd_cutadapt_tag_str
            else:
                args = shlex.split(cmd_cutadapt_tag_str)
            run_result = subprocess.run(args=args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            Logger.instance().info(run_result.stdout.decode())

            ########################################################################################
            #
            # Trim primers from output
            # cutadapt --cores=8 --no-indels --error-rate 0.1 --minimum-length 50 --maximum-length 500 --trimmed-only
            # --front 'TCCACTAATCACAARGATATTGGTAC;min_overlap=26...GGAGGATTTGGWAATTGATTAGTW;min_overlap=24'
            # --output /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_sorted_trimmed_000.fasta
            # /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_sorted_000.fasta
            #
            ########################################################################################

            in_fasta_path = out_fasta_path
            out_fasta_basename = os.path.basename(in_fasta_path).replace(
                f'_sorted_{mergedfasta}', f'_sorted_trimmed_{mergedfasta}')
            out_fasta_path = os.path.join(tempdir, out_fasta_basename)

            cmd_cutadapt_primer_dic = {
                'in_fasta_path': in_fasta_path,
                'out_fasta': out_fasta_path,
                'error_rate': cutadapt_error_rate,
                'num_threads': num_threads,
                'primerFile': primerFile,
            }

            cmd_cutadapt_primer_str = 'cutadapt --cores={num_threads} --no-indels --error-rate {error_rate} ' \
                '--trimmed-only --g file:{primerFile} --output {out_fasta}.{{name}}.fasta {in_fasta_path}'\
                .format(**cmd_cutadapt_primer_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_primer_str))

            if sys.platform.startswith("win"):
                args = cmd_cutadapt_primer_str
            else:
                args = shlex.split(cmd_cutadapt_primer_str)
            run_result = subprocess.run(args=args, stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)

            Logger.instance().info(run_result.stdout.decode())

            ########################################################################################
            #
            # Cut adapt tag of reverse-complement reads
            # cutadapt --cores=8 --no-indels --error-rate 0 --trimmed-only
            # --front 'tgtcgatctacagc;min_overlap=14...acatcgtgatcga;min_overlap=13'
            # --output /tmp/tmpcqlhktae/MFZR1_S4_L001_R1_001_merged_rc_sorted_000.fasta
            # out/control_mfzr/merged/MFZR1_S4_L001_R1_001_merged.fasta
            #
            #######################################################################################
            if no_reverse: #no_reverse stores False, if the option is selected no_reverse == False
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
                run_result = subprocess.run(args=args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

                Logger.instance().info(run_result.stdout.decode())

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
                run_result = subprocess.run(args=args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

                Logger.instance().info(run_result.stdout.decode())

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

            if no_reverse:

                Logger.instance().debug("Pooling fwd and rc reads...")

                if out_final_fasta_path.endswith(".gz"):      
                    _open = partial(gzip.open) 
                elif out_final_fasta_path.endswith(".bz2"):
                    _open = partial(bz2.open)
                else:
                    _open = open

                if out_rc_fasta_path.endswith(".gz"):
                    _open2 = partial(gzip.open) 
                elif out_rc_fasta_path.endswith(".bz2"):
                    _open2 = partial(bz2.open) 
                else: 
                    _open2 = open
                
                with _open(out_final_fasta_path, 'at') as fout:
                    with _open2(out_rc_fasta_path, 'rt') as fin:
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
