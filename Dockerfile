FROM ubuntu:16.04

WORKDIR /app
ADD . /app

RUN \
  apt-get update \
  && apt-get install -y build-essential python python-pip wget vim \
  && apt-get install -y parallel

RUN pip install --upgrade pip
RUN pip install matplotlib

ENTRYPOINT /bin/bash
