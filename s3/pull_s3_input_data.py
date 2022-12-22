import boto3
import botocore
from pathlib import Path
import sys
import os

ROOT_DIR = (Path(__file__).parent / "../").resolve()
sys.path.append(str(ROOT_DIR))

from s3.utils import EnviroVars, AuthError


def main():
    pull_data_from_s3()


def pull_data_from_s3() -> None:
    """
    Download all data files from the shrnaseq folder. Production or test environment specified
    by user pulls files from appropriate sub-folder.

    NOTE: Boto3 will search for the following environment variables:
    - AWS_ACCESS_KEY
    - AWS_SECRET_KEY
    - AWS_BUCKET_NAME
    - AWS_BUCKET_DIR

    Raises
    ------
    AuthError
        Invalid credentials for AWS, check env variables are correct.
    ValueError
        Could not resolve S3 bucket. Check that the target directory on AWS or the AWS bucket name are correct.
    ClientError
        Some other error with AWS client
    """
    data_folder = EnviroVars.AWS_PROD_DIR.value
    if EnviroVars.TEST_ENV.value == "true":
        data_folder = EnviroVars.AWS_TEST_DIR.value

    try:
        session = boto3.Session(
            aws_access_key_id=EnviroVars.AWS_ACCESS_KEY.value,
            aws_secret_access_key=EnviroVars.AWS_SECRET_KEY.value,
        )
        s3 = session.resource("s3")
        bucket = s3.Bucket(EnviroVars.AWS_BUCKET_NAME.value)

        data_folder_prefix = f"{EnviroVars.AWS_BUCKET_DIR.value}/{data_folder}/{EnviroVars.AWS_INPUT_DIR.value}/"

        files = list(bucket.objects.filter(Prefix=data_folder_prefix))[1:]

        os.makedirs(EnviroVars.LOCAL_INPUT_DIR.value, exist_ok=True)

        print("Downloading files from S3...")
        for file in files:
            filename = file.key.split("/")[-1]
            output_path = (
                f"{EnviroVars.LOCAL_INPUT_DIR.value}/{file.key.split('/')[-1]}"
            )
            bucket.download_file(file.key, output_path)
            print(f"Downloaded {filename}")

    except botocore.exceptions.ClientError as e:
        if e.response["Error"]["Code"] == "404":
            raise ValueError(
                "Could not resolve S3 bucket. Check that the target directory on AWS or the AWS bucket name are correct."
            )
        if e.response["Error"]["Code"] == "401":
            raise AuthError(
                "Invalid credentials for AWS, check env variables are correct."
            )
        else:
            raise


if __name__ == "__main__":
    main()
