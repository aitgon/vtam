Installation
=================================================

Linux
-------------------------------------------------

The easiest way to install and run vtam is via conda (`<https://docs.conda.io/projects/conda/en/latest/index.html>`_) environment. Download and install Miniconda `<https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_. Use 3.X python version.

Create a conda environment.

.. code-block:: bash

    conda create --name vtam python=3.7 -y

Activate the conda environment.

.. code-block:: bash

    conda activate vtam

Install Cutadapt, VSEARCH and BLAST in the environment

.. code-block:: bash

    conda install -c bioconda blast
    conda install -c bioconda vsearch

Install VTAM via pip

.. code-block:: bash

    pip install vtam
    
You can also verify the BLAST (>=  v2.2.26), CutAdapt (>= 2.10) and VSEARCH (>= 2.15.0) versions:

.. code-block:: bash

    blastn -version
    cutadapt --version
    vsearch --version

Singularity
-------------------------------------------------

We provide a singularity container in the VTAM github: https://github.com/aitgon/vtam .
First you build the image with this command as root:

.. code-block:: bash

    sudo singularity build vtam.sif vtam.singularity

Then you can use VTAM directly from the singularity image:

.. code-block:: bash

    singularity run --app vtam vtam.sif --help
    singularity run --app vtam vtam.sif merge --help

Windows
-------------------------------------------------

You can run VTAM on a windows machine using the Windows Subsystem for Linux (WSL)

Install Windows Subsystem for Linux and Ubuntu following these instructions: `<https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_

Open WSL Ubuntu in a terminal and install conda using Linux instructions: `<https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html>`_

If it has not been done, you need to install make and gcc:

.. code-block:: bash

    sudo apt-get update
    sudo apt install make
    sudo apt-get install gcc
    conda update -n base -c defaults conda

Go on with the VTAM installation as described here `Linux`_

.. note::
    You can access your files Windows system from the */mnt* directory of your WSL. 
    For example, execute the following command from the Ubuntu terminal to copy the *file.txt* from *c:\temp\file.txt* to your current directory in WSL:

.. code-block:: bash

    cp /mnt/c/temp/file.txt ./file.txt

.. _example_installation:

Test your VTAM installation
-------------------------------------------------

You can verify the installation of vtam by looking at the VTAM version

.. code-block:: bash

    cd vtam
    conda activate vtam
    vtam --version
    
The **vtam example** command downloads the example files and create a file structure.
A snakemake command will run the whole pipeline. If you get through without an error message your VTAM installation is fully functional.

.. code-block:: bash

    vtam example
    cd example
    snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile asper1/user_input/snakeconfig_mfzr.yml --until asvtable_optimized_taxa

