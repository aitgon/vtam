# -*- coding: utf-8 -*-
from unittest import TestCase

from vtam import CommandMerge, CommandSortReads, WopmarsRunner
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager

import filecmp
import os
import pathlib
import shutil
import subprocess
import tarfile
import urllib


class TestMain(TestCase):
    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        url_taxonomy_fastq_tar = "http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/fastq.tar.gz"
        cls.test_outdir_path = os.path.join(PathManager.get_test_path(), 'outdir', cls.__qualname__)
        cls.test_outdir_bak_path = os.path.join(PathManager.get_test_path(), 'test_files', cls.__qualname__)

        cls.fastqdir_path = os.path.join(cls.test_outdir_path, "fastq")
        pathlib.Path(cls.test_outdir_path).mkdir(parents=True, exist_ok=True)

        fastq_tar_path = os.path.join(cls.test_outdir_path, "fastq.tar.gz")

        if not os.path.isfile(fastq_tar_path):
            urllib.request.urlretrieve(url_taxonomy_fastq_tar, fastq_tar_path)

        tar = tarfile.open(fastq_tar_path, "r:gz")
        tar.extractall(path=cls.test_outdir_path)
        tar.close()

    def test_01(self):

        ################################################################################################################
        #
        # Command Merge
        #
        ################################################################################################################

        fastqinfo_path = os.path.join(PathManager.get_package_path(), "doc/data/fastqinfo.tsv")
        # fastqdir =
        fastainfo_path = os.path.join(self.test_outdir_path, "fastainfo.tsv")
        fastadir_path = os.path.join(self.test_outdir_path, "merged")
        CommandMerge.main(fastqinfo=fastqinfo_path, fastqdir=self.fastqdir_path, fastainfo=fastainfo_path,
                          fastadir=fastadir_path)

        fastainfo_bak_path = os.path.join(self.test_outdir_bak_path, "fastainfo.tsv")
        self.assertTrue(filecmp.cmp(fastainfo_path, fastainfo_path, shallow=False))

        ################################################################################################################
        #
        # Command SortReads
        #
        ################################################################################################################

        sorteddir_path = os.path.join(self.test_outdir_path, "sorted")
        sortedreadinfo_path = os.path.join(sorteddir_path, "readinfo.tsv")

        CommandSortReads.main(fastainfo=fastainfo_path, fastadir=fastadir_path, outdir=sorteddir_path)

        sortedreadinfo_bak_path = os.path.join(self.test_outdir_bak_path, "readinfo.tsv")
        self.assertTrue(filecmp.cmp(sortedreadinfo_path, sortedreadinfo_bak_path, shallow=False))

        ################################################################################################################
        #
        # Command Filter
        #
        ################################################################################################################

        asvtable_path = os.path.join(self.test_outdir_path, "asvtable.tsv")
        db_path = os.path.join(self.test_outdir_path, "db.sqlite")

        # Remove if exists, equivalent force
        if os.path.exists(db_path):
            pathlib.Path(db_path).unlink()

        cli_args_dic = {'readinfo': sortedreadinfo_path,
                        'readdir': sorteddir_path,
                        'asvtable': asvtable_path,
                        'db': db_path,
                        'params': None,
                        'threshold_specific': None}

        wopmars_runner = WopmarsRunner(command='filter', cli_args_dic=cli_args_dic)
        wopmars_command = wopmars_runner.get_wopmars_command()

        # Some arguments will be passed through environmental variables
        Logger.instance().info(wopmars_command)
        run_result = subprocess.run(wopmars_command, shell=True)

        asvtable_bak_path = os.path.join(self.test_outdir_bak_path, "asvtable.tsv")
        self.assertTrue(filecmp.cmp(asvtable_path, asvtable_bak_path, shallow=False))

        ################################################################################################################
        #
        # Command Optimize
        #
        ################################################################################################################

        # asvtable_path = os.path.join(self.test_outdir_path, "asvtable.tsv")
        db_path = os.path.join(self.test_outdir_path, "db.sqlite")
        known_occurrences_path = os.path.join(PathManager.get_package_path(), "doc/data/known_occurrences.tsv")

        optimize_lfn_biosample_replicate = os.path.join(self.test_outdir_path, "optimize_lfn_biosample_replicate.tsv")
        optimize_lfn_read_count_and_lfn_variant = os.path.join(self.test_outdir_path, "optimize_lfn_read_count_and_lfn_variant.tsv")
        optimize_lfn_variant_specific = os.path.join(self.test_outdir_path, "optimize_lfn_variant_specific.tsv")
        optimize_pcr_error = os.path.join(self.test_outdir_path, "optimize_pcr_error.tsv")

        # Remove if exists
        if os.path.exists(optimize_lfn_biosample_replicate):
            pathlib.Path(optimize_lfn_biosample_replicate).unlink()
        if os.path.exists(optimize_lfn_read_count_and_lfn_variant):
            pathlib.Path(optimize_lfn_read_count_and_lfn_variant).unlink()
        if os.path.exists(optimize_lfn_variant_specific):
            pathlib.Path(optimize_lfn_variant_specific).unlink()
        if os.path.exists(optimize_pcr_error):
            pathlib.Path(optimize_pcr_error).unlink()

        cli_args_dic = {'readinfo': sortedreadinfo_path,
                        'readdir': sorteddir_path,
                        'outdir': self.test_outdir_path,
                        'known_occurrences': known_occurrences_path,
                        'db': db_path,
                        'params': None}

        wopmars_runner = WopmarsRunner(command='optimize', cli_args_dic=cli_args_dic)
        wopmars_command = wopmars_runner.get_wopmars_command()

        # Some arguments will be passed through environmental variables
        run_result = subprocess.run(wopmars_command, shell=True)

        optimize_lfn_biosample_replicate_bak = os.path.join(self.test_outdir_bak_path, "optimize_lfn_biosample_replicate.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_biosample_replicate, optimize_lfn_biosample_replicate_bak, shallow=False))

        optimize_lfn_read_count_and_lfn_variant_bak = os.path.join(self.test_outdir_bak_path, "optimize_lfn_read_count_and_lfn_variant.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_read_count_and_lfn_variant, optimize_lfn_read_count_and_lfn_variant_bak, shallow=False))

        optimize_lfn_variant_specific_bak = os.path.join(self.test_outdir_bak_path, "optimize_lfn_variant_specific.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_variant_specific, optimize_lfn_variant_specific_bak, shallow=False))

        optimize_pcr_error_bak = os.path.join(self.test_outdir_bak_path, "optimize_pcr_error.tsv")
        self.assertTrue(filecmp.cmp(optimize_pcr_error, optimize_pcr_error_bak, shallow=False))

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.test_outdir_path)
