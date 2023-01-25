# Deployment

This workflow can be deployed to AWS ECS using Terraform. The Terraform scripts are located in the `deployment` directory.

## Prerequisites

- [Terraform](https://www.terraform.io/downloads.html) installed
- [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html) installed and configured with an AWS account
- [Docker](https://docs.docker.com/get-docker/) installed

## Usage

NOTE: The following instructions assume that you are in the `deployment` directory.

You will need a `.env` file with variables required by Terraform. The `.env` file should look like this:

```bash
AWS_ACCOUNT_ID=123456789012
AWS_REGION=eu-west-2
ECR_REGISTRY=123456789012.dkr.ecr.eu-west-2.amazonaws.com
ECR_REPOSITORY=shrnaseq
ECR_IMAGE_TAG=latest
S3_BUCKET_NAME=shrnaseq
S3_DATA_PREFIX=data
S3_EXECUTION_FOLDER=execute
ACCESS_KEY_PARAMETER_NAME=shrnaseq-access-key
SECRET_KEY_PARAMETER_NAME=shrnaseq-secret-key
```

- `AWS_ACCOUNT_ID` and `AWS_REGION` variables can be found in the AWS console.
- `ECR_REGISTRY` variable can be found in the AWS ECR console.
- `ECR_REPOSITORY` variable is the name of the repository in ECR.
- `ECR_IMAGE_TAG` variable is the tag of the image to be deployed.
- `S3_BUCKET_NAME` variable is the name of the S3 bucket to be used for storing data and execution logs.
- `S3_DATA_PREFIX` variable is the prefix for the data folder in the S3 bucket.
- `S3_EXECUTION_FOLDER` variable is the name of the folder that is created to trigger the workflow.
- `ACCESS_KEY_PARAMETER_NAME` and `SECRET_KEY_PARAMETER_NAME` variables are the names of the AWS Systems Manager Parameter Store parameters to store the AWS access key and secret key.

You will need to push a built docker image to the ECR repository to use it in the ECS task.

Once this is done, run the following commands:

```bash
terraform init
terraform apply
```

Upload the data to the `config` folder in the S3 bucket data folder, as described in the Usage section of the root [README](../README.md).

Create the execution folder referenced above in `config` to trigger the workflow.

You can set the number of CPUs and memory allocated to the workflow by changing the appropriate variables in `variables.tf`.

NOTE: The task currently fails due to an issue with snakemake remote S3 because the organism db is not found. This is a known issue and will be fixed in a future release.