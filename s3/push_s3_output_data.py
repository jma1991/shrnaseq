import boto3
import botocore
from pathlib import Path
import sys
import os

ROOT_DIR = (Path(__file__).parent / "../").resolve()
sys.path.append(str(ROOT_DIR))

from s3.utils import EnviroVars, AuthError


def main():
    push_data_to_s3()


def push_data_to_s3() -> None:
    """
    Push all data files from output folder to S3.

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
        s3 = session.client("s3")

        data_folder_prefix = f"{EnviroVars.AWS_BUCKET_DIR.value}/{data_folder}/"

        output_roots = ["results", "plots", "logs"]

        print("Uploading files to S3...")
        for root, _, files in os.walk("output"):
            is_output_root = any(
                [root.startswith(f"output/{folder}") for folder in output_roots]
            )
            if is_output_root:
                for file in files:
                    local_path = os.path.join(root, file)
                    aws_path = os.path.join(
                        data_folder_prefix, os.path.join(root, file)
                    )
                    s3.upload_file(
                        local_path, EnviroVars.AWS_BUCKET_NAME.value, aws_path
                    )
                    print(f"Uploaded {local_path}")

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
