SHELL=/bin/bash

CONDABIN=$$(dirname `which python`)

install: install_wopmars install_deps install_vtam

help:
	@echo "Usage: make"

WOPMARSVERSION=0.0.11

install_wopmars:
	wget https://github.com/aitgon/wopmars/archive/${WOPMARSVERSION}.tar.gz -O wopmars-${WOPMARSVERSION}.tar.gz
	tar zxvf wopmars-${WOPMARSVERSION}.tar.gz
	${CONDABIN}/pip install wopmars-${WOPMARSVERSION}/. --upgrade
	rm -f wopmars-${WOPMARSVERSION}.tar.gz
	rm -rf wopmars-${WOPMARSVERSION}

install_deps:
	conda install -c bioconda vsearch=2.7.0 -y
	conda install -c bioconda blast=2.9.0 -y
	pip install cutadapt --upgrade

install_vtam:
	${CONDABIN}/pip install -e . --upgrade

