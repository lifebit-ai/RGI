FROM python:3.8.0-slim-buster

LABEL authors="ines@lifebit.ai" \
      description="Docker image containing python dependencies for RGI-nf"

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && \
    apt-get install -y \
    procps build-essential && \
    apt-get clean -y

# Install python dependencies
RUN pip install --upgrade pip setuptools wheel

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt