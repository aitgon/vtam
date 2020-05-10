# TODO fix it Mai 9, 2020
# -*- coding: utf-8 -*-
import shlex
import sys
import filecmp
import os
import pathlib
import shutil
import subprocess
import tarfile
import urllib
import unittest

from vtam.utils.constants import fastq_tar_gz_url
from vtam import CommandMerge, CommandSortReads, WopmarsRunner
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from urllib import request


@unittest.skipIf(request.urlopen(fastq_tar_gz_url).getcode() != 200,
                 "Test dataset not available online!")
class TestTutorialCommands(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        # vtam needs to be in the path
        subprocess.run([sys.executable, '-m', 'pip', 'install', '{}/.'.format(PathManager.get_package_path()),
                        '--upgrade'])

        cls.test_outdir_path = os.path.join(PathManager.get_test_path(), 'outdir')
        shutil.rmtree(cls.test_outdir_path, ignore_errors=True)  # during development of the test, this prevents errors
        pathlib.Path(cls.test_outdir_path).mkdir(parents=True, exist_ok=True)

        ################################################################################################################
        #
        # Download fastq test dataset
        #
        ################################################################################################################

        fastq_tar_path = os.path.join(cls.test_outdir_path, "fastq.tar.gz")
        if not os.path.isfile(fastq_tar_path):
            urllib.request.urlretrieve(fastq_tar_gz_url, fastq_tar_path)
        tar = tarfile.open(fastq_tar_path, "r:gz")
        tar.extractall(path=cls.test_outdir_path)
        tar.close()

        # Set test paths
        cls.fastqinfo_path = os.path.join(PathManager.get_package_path(), "doc/data/dryad.f40v5_small/fastqinfo.tsv")
        cls.fastqdir_path = os.path.join(cls.test_outdir_path, "fastq")
        cls.fastainfo_path = os.path.join(cls.test_outdir_path, "fastainfo.tsv")
        cls.fastadir_path = os.path.join(cls.test_outdir_path, "merged")

        cls.sorted_dir_path = os.path.join(cls.test_outdir_path, "sorted")
        cls.sortedreadinfo_path = os.path.join(cls.sorted_dir_path, "readinfo.tsv")

        cls.log_path = os.path.join(cls.test_outdir_path, "vtam.log")

        cls.asvtable_path = os.path.join(cls.test_outdir_path, "asvtable_default.tsv")

        cls.args = {}
        cls.args['fastqinfo'] = cls.fastqinfo_path
        cls.args['fastqdir'] = cls.fastqdir_path
        cls.args['fastainfo'] = cls.fastainfo_path
        cls.args['fastadir'] = cls.fastadir_path
        cls.args['sorted'] = cls.sorted_dir_path
        cls.args['db'] = os.path.join(cls.test_outdir_path, "db.sqlite")
        cls.args['readinfo'] = cls.sortedreadinfo_path
        cls.args['readdir'] = cls.sorted_dir_path
        cls.args['asvtable'] = cls.asvtable_path
        cls.args['log'] = cls.log_path

    def test_step01_merge(self):

        ################################################################################################################
        #
        # Command Merge
        #
        ################################################################################################################

        cmd_merge = "vtam merge --fastqinfo {fastqinfo} --fastqdir {fastqdir} --fastainfo {fastainfo} --fastadir {fastadir} " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd_merge))

        self.fastainfo_path_bak = os.path.join(os.path.dirname(__file__), "fastainfo.tsv")
        self.fastadir_path_bak = os.path.join(os.path.dirname(__file__), "merge")
        self.assertTrue(filecmp.cmp(self.fastainfo_path, self.fastainfo_path_bak, shallow=True))
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'mfzr_1_fw.fasta')) >= 11608260)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'mfzr_1_fw.fasta')) <= 11608270)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'zfzr_3_fw.fasta')) >= 11658700)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'zfzr_3_fw.fasta')) <= 11658710)

    def test_step02_sort_reads(self):

        ################################################################################################################
        #
        # Command SortReads
        #
        ################################################################################################################

        cmd = "vtam sortreads --fastainfo {fastainfo} --fastadir {fastadir} --outdir {sorted} " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

        self.sortedreadinfo_path_bak = os.path.join(os.path.dirname(__file__), "sortedreadinfo.tsv")
        self.assertTrue(filecmp.cmp(self.sortedreadinfo_path, self.sortedreadinfo_path_bak, shallow=True))
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'mfzr_1_fw_000.fasta')) >= 5131890)  # 5131896
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'mfzr_1_fw_000.fasta')) <= 5131900)
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'zfzr_3_fw_023.fasta')) >= 909500)  # 909507
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'zfzr_3_fw_023.fasta')) <= 909510)

    def test_step03_filter(self):

        ################################################################################################################
        #
        # Command Filter
        #
        ################################################################################################################

        cmd = "vtam filter --db {db} --readinfo {readinfo} --readdir {readdir} --asvtable {asvtable} " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

        asvtable_bak_path = os.path.join(os.path.dirname(__file__), "asvtable_default.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    # def test_step04_optimize(self):
    #
    #     ################################################################################################################
    #     #
    #     # Command Optimize
    #     #
    #     ################################################################################################################
    #
    #     # asvtable_path = os.tsv_path.join(self.outdir_path, "asvtable.tsv")
    #     db_path = os.path.join(self.test_outdir_path, "db.sqlite")
    #     known_occurrences_path = os.path.join(PathManager.get_package_path(), "doc/data/known_occurrences.tsv")
    #
    #     cli_args_dic = {'readinfo': self.sortedreadinfo_path,
    #                     'readdir': self.sorteddir_path,
    #                     'outdir': self.test_outdir_path,
    #                     'known_occurrences': known_occurrences_path,
    #                     'db': db_path,
    #                     'params': None,
    #                     'log_file': self.log_path}
    #
    #     wopmars_runner = WopmarsRunner(command='optimize', cli_args_dic=cli_args_dic)
    #     wopmars_command = wopmars_runner.get_wopmars_command()
    #
    #     # Some arguments will be passed through environmental variables
    #     Logger.instance().info(wopmars_command)
    #
    #     run_result = subprocess.run(wopmars_command, shell=True)
    #
    #     optimize_lfn_biosample_replicate = os.path.join(self.test_outdir_path, "optimize_lfn_biosample_replicate.tsv")
    #     optimize_lfn_read_count_and_lfn_variant = os.path.join(self.test_outdir_path, "optimize_lfn_read_count_and_lfn_variant.tsv")
    #     optimize_lfn_variant_specific = os.path.join(self.test_outdir_path, "optimize_lfn_variant_specific.tsv")
    #     optimize_pcr_error = os.path.join(self.test_outdir_path, "optimize_pcr_error.tsv")
    #
    #     # Remove if exists
    #     if os.path.exists(optimize_lfn_biosample_replicate):
    #         pathlib.Path(optimize_lfn_biosample_replicate).unlink()
    #     if os.path.exists(optimize_lfn_read_count_and_lfn_variant):
    #         pathlib.Path(optimize_lfn_read_count_and_lfn_variant).unlink()
    #     if os.path.exists(optimize_lfn_variant_specific):
    #         pathlib.Path(optimize_lfn_variant_specific).unlink()
    #     if os.path.exists(optimize_pcr_error):
    #         pathlib.Path(optimize_pcr_error).unlink()
    #
    #     # This second run is necessary but I do not know why
    #     run_result = subprocess.run(wopmars_command, shell=True)
    #
    #     optimize_lfn_biosample_replicate_bak = os.path.join(self.test_outdir_bak_path, "optimize_lfn_biosample_replicate.tsv")
    #     # TODO must be fixed
    #     # self.assertTrue(filecmp.cmp(optimize_lfn_biosample_replicate, optimize_lfn_biosample_replicate_bak, shallow=False))
    #
    #     optimize_lfn_read_count_and_lfn_variant_bak = os.path.join(self.test_outdir_bak_path, "optimize_lfn_read_count_and_lfn_variant.tsv")
    #     # TODO must be fixed
    #     # self.assertTrue(filecmp.cmp(optimize_lfn_read_count_and_lfn_variant, optimize_lfn_read_count_and_lfn_variant_bak, shallow=False))
    #
    #     optimize_lfn_variant_specific_bak = os.path.join(self.test_outdir_bak_path, "optimize_lfn_variant_specific.tsv")
    #     # TODO must be fixed
    #     # self.assertTrue(filecmp.cmp(optimize_lfn_variant_specific, optimize_lfn_variant_specific_bak, shallow=False))
    #
    #     optimize_pcr_error_bak = os.path.join(self.test_outdir_bak_path, "optimize_pcr_error.tsv")
    #     # TODO must be fixed
    #     # self.assertTrue(filecmp.cmp(optimize_pcr_error, optimize_pcr_error_bak, shallow=False))

    @classmethod
    def tearDownClass(cls):

        subprocess.run([sys.executable, '-m', 'pip', 'uninstall', 'vtam', '-y'])
        shutil.rmtree(cls.test_outdir_path, ignore_errors=True)

