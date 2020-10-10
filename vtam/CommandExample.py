import os
import shlex
import subprocess
import sys
import tarfile
import urllib

from vtam.utils.PathManager import PathManager
from urllib import request
from vtam.utils.constants import fastq_tar_gz_url2, fastq_tar_gz_url1
from vtam.utils.MyProgressBar import MyProgressBar


class CommandExample(object):
    """Generates data in the outdir folder that can be used for a quick example"""

    @staticmethod
    def main(outdir):

        package_path = PathManager.get_package_path()

        ############################################################################################
        #
        # Download fastq
        #
        ############################################################################################

        # fastq_tar_path = os.path.join(outdir, "fastq.tar.gz")
        # if not os.path.isfile(fastq_tar_path):
        #     Logger.instance().info(fastq_tar_path)
        #     urllib.request.urlretrieve(fastq_tar_gz_url1, fastq_tar_path, MyProgressBar())
        # tar = tarfile.open(fastq_tar_path, "r:gz")
        # tar.extractall(path=outdir)
        # tar.close()

        fastq_tar_path = os.path.join(package_path, "..", "data", "fastq.tar.gz")  # tests
        if not os.path.isfile(fastq_tar_path):  # in the distribution
            fastq_tar_path = os.path.join(outdir, "fastq.tar.gz")
            try:
                urllib.request.urlretrieve(fastq_tar_gz_url1, fastq_tar_path, MyProgressBar())
            except:
                urllib.request.urlretrieve(fastq_tar_gz_url2, fastq_tar_path, MyProgressBar())
        tar = tarfile.open(fastq_tar_path, "r:gz")
        tar.extractall(path=outdir)
        tar.close()

        os.remove(fastq_tar_path)

        ############################################################################################
        #
        # Set command args
        #
        ############################################################################################

        args = {}
        args['package_path'] = package_path
        args['snake_tuto_data'] = os.path.join(package_path, "data/snake.tuto.data.yml")

        ############################################################################################
        #
        # Copy data to directory tree
        #
        ############################################################################################

        cmd = "snakemake --cores 1 -s {snake_tuto_data} --config MARKER=mfzr " \
              "PROJECT=asper1 PACKAGE_PATH={package_path} --until all_one_marker".format(**args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=outdir)

