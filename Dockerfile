FROM condaforge/mambaforge:4.14.0-0

COPY ["environment.yml", "./"]

RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends pigz tzdata procps nano bc && \
    apt-get clean 

RUN mamba env update --name base --file environment.yml && conda clean -a -q -y
