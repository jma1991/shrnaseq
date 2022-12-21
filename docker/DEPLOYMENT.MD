# shrnaseq deployment

## Overview

Deployment of the shranseq workflow is controlled via an orchestrated series of docker containers. 

Input data is pulled from an AWS S3 bucket within a Python container, transferred to a Snakemake container for processing and then output data transferred to a Python container to be uploaded to the S3 bucket.

## Setting up

AWS credentials are required to pull data from an S3 bucket. You will need to generate the following for a role with S3 read and write permissions:

- AWS_ACCESS_KEY
- AWS_SECRET_KEY

These should be stored in the [.env](./.env) file as outlined [here](https://docs.docker.com/compose/env-file/).

docker can be installed by following the instructions appropriate for your OS [here](https://docs.docker.com/get-docker/).

If you're using Windows, it's recommended that you install [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) and configure an [Ubuntu](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#4-configure-ubuntu) environment. 

You can then install the [docker engine](https://docs.docker.com/engine/install/ubuntu/). Make sure that if you're running docker as a non-root user, which is recommended, that you add the user to the docker group as outlined [here](https://docs.docker.com/engine/install/linux-postinstall/).

## Usage

To run the pipeline, you should execute the [run-docker.sh](../run-docker.sh) file from the root directory.