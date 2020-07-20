FROM continuumio/miniconda3

# docker build  -t vtam .
# docker run --rm -v $(pwd):/apps  -it vtam vtam merge --help
# docker run --rm -v $(pwd):/apps  -it vtam vtam sortreads --help

MAINTAINER Aitor Gonzalez <aitor.gonzalez@univ-amu.fr>

WORKDIR /app
RUN echo "Europe/Paris" > /etc/timezone && \
    apt-get update && \
    apt-get install -y build-essential python3-dev

RUN conda update -n base -c defaults conda
RUN conda create -n vtam python=3.7
RUN echo "source activate vtam" > ~/.bashrc

# Copy VTAM project
RUN mkdir -p /app/deps/vtam
COPY ./ /app/deps/vtam/
RUN conda env update -n vtam -f /app/deps/vtam/environment-dev.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "vtam", "/bin/bash", "-c"]
RUN python3.7 -m pip install /app/deps/vtam/
ENTRYPOINT ["conda", "run", "-n", "vtam"]

