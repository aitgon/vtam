from bdb import set_trace
from vtam.utils import pip_install_vtam_for_tests
from vtam.utils.PathManager import PathManager
import filecmp
import os
import pathlib
import shlex
import shutil
import subprocess
import sys
import tarfile
import unittest
import urllib.request

from vtam.utils.constants import fastq_tar_gz_url1, fastq_tar_gz_url2, fastq_tar_gz_url3
# from vtam.utils.MyProgressBar import MyProgressBar
from tqdm import tqdm
from vtam.utils import tqdm_hook


class TestTutorialCommands(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        pip_install_vtam_for_tests()  # vtam needs to be in the path

        cls.package_path = PathManager.get_package_path()
        cls.test_path = os.path.join(PathManager.get_test_path())

        cls.outdir_path = os.path.join(cls.test_path, 'outdir')
        cls.outdir_data_path = os.path.join(cls.outdir_path, 'data')
        shutil.rmtree(cls.outdir_path, ignore_errors=True)  # during development of the test, this prevents errors
        pathlib.Path(cls.outdir_data_path).mkdir(parents=True, exist_ok=True)

        ############################################################################################
        #
        # Download fastq test dataset
        #
        ############################################################################################

        fastq_tar_path = os.path.join(cls.outdir_data_path, "fastq.tar.gz")
        # Test first in local dir, otherwise in the remote URLs
        if not os.path.isfile(fastq_tar_path) or pathlib.Path(fastq_tar_path).stat().st_size < 1000000:
            try:
                # urllib.request.urlretrieve(fastq_tar_gz_url1, fastq_tar_path, MyProgressBar())
                with tqdm(...) as t:
                    t.set_description(os.path.basename(fastq_tar_path))
                    urllib.request.urlretrieve(fastq_tar_gz_url1, fastq_tar_path, reporthook=tqdm_hook(t))
            except Exception:
                try:
                    # urllib.request.urlretrieve(fastq_tar_gz_url2, fastq_tar_path, MyProgressBar())
                    with tqdm(...) as t:
                        t.set_description(os.path.basename(fastq_tar_path))
                        urllib.request.urlretrieve(fastq_tar_gz_url2, fastq_tar_path, reporthook=tqdm_hook(t))
                except Exception:
                    # urllib.request.urlretrieve(fastq_tar_gz_url3, fastq_tar_path, MyProgressBar())
                    with tqdm(...) as t:
                        t.set_description(os.path.basename(fastq_tar_path))
                        urllib.request.urlretrieve(fastq_tar_gz_url3, fastq_tar_path, reporthook=tqdm_hook(t))

        tar = tarfile.open(fastq_tar_path, "r:gz")
        tar.extractall(path=cls.outdir_path)
        tar.close()

        # Set test paths
        cls.fastqinfo_path = os.path.join(PathManager.get_package_path(), "data/example/fastqinfo.tsv")
        cls.fastqdir_path = os.path.join(cls.outdir_path, "fastq")
        cls.fastainfo_path = os.path.join(cls.outdir_path, "fastainfo.tsv")
        cls.fastadir_path = os.path.join(cls.outdir_path, "merged")

        cls.sorted_dir_path = os.path.join(cls.outdir_path, "sorted")
        cls.sortedinfo_path = os.path.join(cls.sorted_dir_path, "sortedinfo.tsv")

        cls.log_path = os.path.join(cls.outdir_path, "vtam.log")

        cls.asvtable_path = os.path.join(cls.outdir_path, "asvtable_default.tsv")

        cls.args = {}
        cls.args['fastqinfo'] = cls.fastqinfo_path
        cls.args['fastqdir'] = cls.fastqdir_path
        cls.args['fastainfo'] = cls.fastainfo_path
        cls.args['fastadir'] = cls.fastadir_path
        cls.args['sorted'] = cls.sorted_dir_path
        cls.args['db'] = os.path.join(cls.outdir_path, "db.sqlite")
        cls.args['sortedinfo'] = cls.sortedinfo_path
        cls.args['sorteddir'] = cls.sorted_dir_path
        cls.args['asvtable'] = cls.asvtable_path
        cls.args['log'] = cls.log_path

        ################################################################################################################
        #
        # Command Merge
        #
        ################################################################################################################

        cmd = "vtam merge --fastqinfo {fastqinfo} --fastqdir {fastqdir} --fastainfo {fastainfo} --fastadir {fastadir} " \
              "-v --log {log}".format(**cls.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args)

    def test_step01_merge(self):

        self.fastainfo_path_bak = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "fastainfo.tsv")
        self.fastadir_path_bak = os.path.join(os.path.dirname(__file__), "merge")

        self.assertTrue(filecmp.cmp(self.fastainfo_path, self.fastainfo_path_bak, shallow=True))
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'mfzr_1_fw.fasta')) >= 11608260)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'mfzr_1_fw.fasta')) <= 11795030)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'zfzr_3_fw.fasta')) >= 11658700)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'zfzr_3_fw.fasta')) <= 11838710)

    def test_step02_sortreads(self):

        ################################################################################################################
        #
        # Command SortReads
        #
        ################################################################################################################

        cmd = "vtam sortreads --fastainfo {fastainfo} --fastadir {fastadir} --sorteddir {sorted} " \
              "-v --log {log}".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args)

        self.sortedinfo_path_bak = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "sortedinfo.tsv")
        

        self.assertTrue(filecmp.cmp(self.sortedinfo_path, self.sortedinfo_path_bak, shallow=True))

        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'mfzr_1_fw_trimmed_0_tpos1_run1.fasta')) >= 5131890)  # 5131896 linux, 5155350 windows
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'mfzr_1_fw_trimmed_0_tpos1_run1.fasta')) <= 5155360)
        
        import pdb; pdb.set_trace()
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'zfzr_3_fw_trimmed_0_14ben02.fasta')) >= 909500)  # 909507 linux, 913883 windows
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'zfzr_3_fw_trimmed_0_14ben02.fasta')) <= 913890)
        # 'zfzr_3_fw_023.fasta'

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)
