FROM python:3.7-buster

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -qqy && apt-get install -qqy \
  make \
  clang \
  curl \
  gcc \
  g++ \
  libgsl23 \
  libgsl-dev
 
RUN python3.7 -m pip install numpy pandas pytest

COPY tests/requirements.txt /requirements.txt
RUN python3 -m pip install -r /requirements.txt

COPY . /IBM_simul

WORKDIR /IBM_simul
RUN cd src && make clean && make all

WORKDIR /IBM_simul
ENTRYPOINT ["/usr/local/bin/pytest", "-s"]
