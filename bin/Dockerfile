FROM python:3.12-slim

RUN apt-get update \
    && apt-get install -y procps

COPY ./requirements.txt /tmp/requirements.txt

RUN pip install --no-cache-dir --prefer-binary -r /tmp/requirements.txt
