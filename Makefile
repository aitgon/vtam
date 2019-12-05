SHELL=/bin/bash

REPOSITORIES=${HOME}/Software/repositories

CONDABIN=$$(dirname `which python`)
CONDAENV=$$(conda env list |grep "*" |cut -d ' ' -f 1,1) # current conda env

install: install_wopmars

help:
	@echo "Usage: make install REPOSITORIES=${REPOSITORIES}"

wopmars_branch=develop
wopmars_sha=116273dc3e50c72b262698c7a95d1f9120199283
wopmars_path=${REPOSITORIES}/wopmars

install_wopmars:
	#${CONDABIN}/pip install git+https://github.com/aitgon/wopmars.git@${wopmars_version}
	cd ${wopmars_path}; git fetch origin; git checkout ${wopmars_branch}; git pull origin ${wopmars_branch}; git reset --hard $(wopmars_sha)
	cd ${wopmars_path}; ${CONDABIN}/pip install -e . --upgrade


