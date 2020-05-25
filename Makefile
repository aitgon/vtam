SHELL=/bin/bash

CONDABIN=$$(dirname `which python`)

install: install_deps install_vtam

help:
	@echo "Usage: make"

install_deps:
	conda env update -f environment.yml

install_vtam:
	${CONDABIN}/pip install -e . --upgrade

