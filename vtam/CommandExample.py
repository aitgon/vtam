import os
import pathlib
import shlex
import subprocess
import sys
import tarfile
import urllib

from vtam.utils.PathManager import PathManager
from urllib import request
from vtam.utils.constants import fastq_tar_gz_url1, fastq_tar_gz_url2, fastq_tar_gz_url3
# from vtam.utils.MyProgressBar import MyProgressBar
from tqdm import tqdm
from vtam.utils import tqdm_hook


class CommandExample(object):
    """Generates data in the outdir folder that can be used for a quick example"""

    @staticmethod
    def main(outdir):

        package_path = PathManager.get_package_path()
        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

        #######################################################################
        #
        # Download fastq
        #
        #######################################################################

        fastq_tar_path = os.path.join(outdir, "fastq.tar.gz")
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
        tar.extractall(path=outdir)
        tar.close()

        os.remove(fastq_tar_path)

        #######################################################################
        #
        # Set command args
        #
        #######################################################################

        args = {}
        args['package_path'] = package_path
        args['snake_tuto_data'] = os.path.join(package_path, "data/snake.tuto.data.yml")

        #######################################################################
        #
        # Copy data to directory tree
        #
        #######################################################################

        cmd = "snakemake --cores 1 -s {snake_tuto_data} --config MARKER=mfzr " \
              "PROJECT=asper1 PACKAGE_PATH={package_path} --until all_one_marker".format(**args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=outdir)

