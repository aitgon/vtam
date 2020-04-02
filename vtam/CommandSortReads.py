import multiprocessing
import os
import shlex
import shutil
import subprocess

import pandas
import pathlib

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from vtam.utils.PathManager import PathManager
from vtam.utils.Logger import Logger


class CommandSortReads(object):
    """Class for the Merge command"""

    @classmethod
    def main(cls, fastainfo, fastadir, outdir, params=None, num_threads=multiprocessing.cpu_count()):

        ################################################################################################################
        #
        # Loop over fasta files to sort reads per fasta
        #
        ################################################################################################################

        fastainfo_df = pandas.read_csv(fastainfo, sep='\t', header=0)
        fastainfo_df.columns = fastainfo_df.columns.str.lower()

        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

        tempdir = PathManager.instance().get_tempdir()
        sorted_read_info_df = pandas.DataFrame()

        read_min_length = 150
        read_max_length = 200

        for i, fasta_filename in enumerate(sorted(fastainfo_df.fasta.unique().tolist())):

            Logger.instance().debug("Analysing FASTA file: {}".format(fasta_filename))

            fasta_info_df_i = fastainfo_df.loc[fastainfo_df.fasta == fasta_filename]

            tag_fwd = fasta_info_df_i.tagfwd[0]
            tag_rev = fasta_info_df_i.tagrev[0]
            primer_fwd = fasta_info_df_i.primerfwd[0]
            primer_rev = fasta_info_df_i.primerrev[0]
            in_fasta_basename = fasta_info_df_i.fasta[0]
            in_raw_fasta_path = os.path.join(fastadir, in_fasta_basename)

            ############################################################################################################
            #
            # Cut adapt tag of forward reads
            #
            ############################################################################################################

            tag_rev_rc = str(Seq(tag_rev, generic_dna).reverse_complement())
            out_fasta_basename = os.path.basename(in_raw_fasta_path).replace('.fasta', '_sorted_%03d.fasta' % i)
            out_fasta_path = os.path.join(tempdir, out_fasta_basename)

            cmd_cutadapt_tag_dic = {'tag_fwd': tag_fwd, 'tag_fwd_len': len(tag_fwd), 'tag_rev_rc': tag_rev_rc,
                        'tag_rev_rc_len': len(tag_rev_rc), 'num_threads': num_threads, 'in_fasta_path': in_raw_fasta_path,
                                'out_fasta': out_fasta_path}
            cmd_cutadapt_tag_str = "cutadapt --cores={num_threads} --no-indels -e 0 --trimmed-only " \
                               "-g '{tag_fwd};min_overlap={tag_fwd_len}...{tag_rev_rc};min_overlap={tag_rev_rc_len}' " \
                               "-o {out_fasta} {in_fasta_path}".format(**cmd_cutadapt_tag_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_tag_str))
            run_result = subprocess.run(shlex.split(cmd_cutadapt_tag_str), stdout=subprocess.PIPE)
            Logger.instance().info(str(run_result.stdout))

            ############################################################################################################
            #
            # Trim primers from output
            # (Error allowed, 0.1 by default; minimum overlap: length of the primer; no indel; multithreading)
            #
            ############################################################################################################
            #   cutadapt --cores=THREADS --no-indels -e ERROR_RATE -m MINIMUM_LENGTH -M MAXIMUM_LENGTH --trimmed-only
            #   -g "PRIMER_FW_SEQ;min_overlap=LENGTH_OF_PRIMER_FW...PRIMER_REV_REVCOMP_SEQ;min_overlap=LENGTH_OF_PRIMER_REV"
            #   -o TRIMMED_MARKER_RUN_SAMPLE_REPLICATE.fasta DEMULT_MARKER_RUN_SAMPLE_REPLICATE.fasta

            primer_rev_rc = str(Seq(primer_rev, generic_dna).reverse_complement())
            in_fasta_path = out_fasta_path
            out_fasta_basename = os.path.basename(in_fasta_path).replace('_sorted_%03d.fasta'%i, '_sorted_trimmed_%03d.fasta'%i)
            out_fasta_path = os.path.join(tempdir, out_fasta_basename)

            cmd_cutadapt_primer_dic = {'primer_fwd': primer_fwd, 'primer_fwd_len': len(primer_fwd), 'primer_rev_rc': primer_rev_rc,
                        'primer_rev_rc_len': len(primer_rev_rc), 'num_threads': num_threads, 'in_fasta_path': in_fasta_path,
                                'out_fasta': out_fasta_path, 'error_rate': 0.1, 'read_min_length': read_min_length,
                                       'read_max_length': read_max_length}
            cmd_cutadapt_primer_str = "cutadapt --cores={num_threads} --no-indels -e {error_rate} -m {read_min_length} " \
                                      "-M {read_max_length} --trimmed-only  " \
                                      "-g '{primer_fwd};min_overlap={primer_fwd_len}...{primer_rev_rc};min_overlap={primer_rev_rc_len}' " \
                               "-o {out_fasta} {in_fasta_path}".format(**cmd_cutadapt_primer_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_primer_str))
            run_result = subprocess.run(shlex.split(cmd_cutadapt_primer_str), stdout=subprocess.PIPE)
            Logger.instance().info(str(run_result.stdout))

            ############################################################################################################
            #
            # Cut adapt tag of reverse-complement reads
            #
            ############################################################################################################

            tag_fwd_rc = str(Seq(tag_fwd, generic_dna).reverse_complement())
            out_rc_fasta_basename= os.path.basename(in_raw_fasta_path).replace('.fasta', '_rc_sorted_%03d.fasta' % i)
            out_rc_fasta_path = os.path.join(tempdir, out_rc_fasta_basename)

            cmd_cutadapt_tag_dic = {'tag_fwd': tag_rev, 'tag_fwd_len': len(tag_rev), 'tag_rev_rc': tag_fwd_rc,
                        'tag_rev_rc_len': len(tag_fwd_rc), 'num_threads': num_threads, 'in_fasta_path': in_raw_fasta_path,
                                'out_fasta': out_rc_fasta_path}
            cmd_cutadapt_tag_str = "cutadapt --cores={num_threads} --no-indels -e 0 --trimmed-only " \
                               "-g '{tag_fwd};min_overlap={tag_fwd_len}...{tag_rev_rc};min_overlap={tag_rev_rc_len}' " \
                               "-o {out_fasta} {in_fasta_path}".format(**cmd_cutadapt_tag_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_tag_str))
            run_result = subprocess.run(shlex.split(cmd_cutadapt_tag_str), stdout=subprocess.PIPE)
            Logger.instance().info(str(run_result.stdout))

            ############################################################################################################
            #
            # Trim primers from output
            # (Error allowed, 0.1 by default; minimum overlap: length of the primer; no indel; multithreading)
            #
            ############################################################################################################
            #   cutadapt --cores=THREADS --no-indels -e ERROR_RATE -m MINIMUM_LENGTH -M MAXIMUM_LENGTH --trimmed-only
            #   -g "PRIMER_FW_SEQ;min_overlap=LENGTH_OF_PRIMER_FW...PRIMER_REV_REVCOMP_SEQ;min_overlap=LENGTH_OF_PRIMER_REV"
            #   -o TRIMMED_MARKER_RUN_SAMPLE_REPLICATE.fasta DEMULT_MARKER_RUN_SAMPLE_REPLICATE.fasta

            primer_fwd_rc = str(Seq(primer_fwd, generic_dna).reverse_complement())
            in_fasta_path = out_rc_fasta_path
            out_rc_fasta_basename = os.path.basename(in_fasta_path).replace('_rc_sorted_%03d.fasta'%i, '_rc_sorted_trimmed_%03d.fasta'%i)
            out_rc_fasta_path = os.path.join(tempdir, out_rc_fasta_basename)

            cmd_cutadapt_primer_dic = {'primer_fwd': primer_rev, 'primer_fwd_len': len(primer_rev), 'primer_rev_rc': primer_fwd_rc,
                        'primer_rev_rc_len': len(primer_fwd_rc), 'num_threads': num_threads, 'in_fasta_path': in_fasta_path,
                                'out_fasta': out_rc_fasta_path, 'error_rate': 0.1, 'read_min_length': read_min_length,
                                       'read_max_length': read_max_length}
            cmd_cutadapt_primer_str = "cutadapt --cores={num_threads} --no-indels -e {error_rate} -m {read_min_length} " \
                                      "-M {read_max_length} --trimmed-only  " \
                                      "-g '{primer_fwd};min_overlap={primer_fwd_len}...{primer_rev_rc};min_overlap={primer_rev_rc_len}' " \
                               "-o {out_fasta} {in_fasta_path}".format(**cmd_cutadapt_primer_dic)

            Logger.instance().debug("Running: {}".format(cmd_cutadapt_primer_str))
            run_result = subprocess.run(shlex.split(cmd_cutadapt_primer_str), stdout=subprocess.PIPE)
            Logger.instance().info(str(run_result.stdout))

            ############################################################################################################
            #
            # Reverse complement back rc fasta and pool
            #
            ############################################################################################################

            out_final_fasta_basename = os.path.basename(in_raw_fasta_path)
            out_final_fasta_path = os.path.join(outdir, out_final_fasta_basename)
            shutil.copy(out_fasta_path, out_final_fasta_path)

            Logger.instance().debug("Pooling fwd and rc reads...")
            with open(out_final_fasta_path, 'a') as fout:
                with open(out_rc_fasta_path, 'r') as fin:
                    for line in fin:
                        if not line.startswith('>'):
                            fout.write("%s\n"%str(Seq(line.strip(), generic_dna).reverse_complement()))
                        else:
                            fout.write(line)

            fasta_info_df_i.drop(['tagfwd', 'primerfwd', 'tagrev', 'primerrev', 'fastqfwd', 'fastqrev'], axis=1, inplace=True)
            fasta_info_df_i['sorted'] = out_final_fasta_basename

        fasta_trimmed_info_tsv = os.path.join(outdir, 'readinfo.tsv')
        sorted_read_info_df.to_csv(fasta_trimmed_info_tsv, sep="\t", header=True, index=False)
