SHELL=/bin/bash

CONDABIN=$$(dirname `which python`)

install: install_wopmars install_deps install_vtam

help:
	@echo "Usage: make"

install_wopmars:
	wget https://github.com/aitgon/wopmars/archive/0.0.9.tar.gz -O wopmars-0.0.9.tar.gz
	tar zxvf wopmars-0.0.9.tar.gz
	${CONDABIN}/pip install wopmars-0.0.9/. --upgrade
	rm -f wopmars-0.0.9.tar.gz
	rm -rf wopmars-0.0.9

install_deps:
	conda install -c bioconda vsearch=2.7.0 -y
	conda install -c bioconda blast=2.9.0 -y

install_vtam:
	${CONDABIN}/pip install -e . --upgrade

