FROM python:3.10-slim as base

# Setup env
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

FROM base AS python-deps

COPY requirements-py.txt .
RUN pip install -r requirements-py.txt

# Create and switch to a new user
RUN useradd --create-home python-s3-user
WORKDIR /home/python-s3-user

# Copy scripts into container
COPY ./s3 ./s3