import os
from enum import Enum
from dotenv import load_dotenv

load_dotenv("../.env")


class AuthError(Exception):
    pass


class MissingEnvVarError(Exception):
    pass


class EnviroVars(Enum):
    """ """

    try:
        AWS_ACCESS_KEY = os.environ["AWS_ACCESS_KEY"]
        AWS_SECRET_KEY = os.environ["AWS_SECRET_KEY"]
    except KeyError:
        raise MissingEnvVarError(
            "Missing AWS access key and secret in environment variables"
        )
    TEST_ENV = os.getenv("TEST_ENV", "false")
    AWS_BUCKET_NAME = os.getenv("AWS_BUCKET_NAME", "zifo-ds-eu").strip("/")
    AWS_BUCKET_DIR = os.getenv("AWS_BUCKET_DIR", "shrnaseq")
    AWS_PROD_DIR = os.getenv("AWS_PROD_DIR", "production_data")
    AWS_TEST_DIR = os.getenv("AWS_TEST_DIR", "test_data")
    AWS_INPUT_DIR = os.getenv("AWS_PROD_DIR", "input")
    LOCAL_INPUT_DIR = os.getenv("LOCAL_INPUT_DIR", "config")
