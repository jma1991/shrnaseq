FROM snakemake/snakemake:stable as base

RUN useradd --create-home mamba-user
WORKDIR /home/mamba-user
USER mamba-user

COPY ./workflow ./workflow
COPY ./config ./config

CMD ["snakemake", "--use-conda", "--cores", "all"]
