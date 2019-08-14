FROM ubuntu:18.04

MAINTAINER Eugene Chow "eugenechow823@gmail.com'

RUN apt-get update -y && apt-get install -y python3-pip python3-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

COPY requirements.txt /requirements.txt
COPY setup.py /setup.py
COPY rg4seeker/ /rg4seeker/
RUN python3 setup.py install

ENTRYPOINT [ "rG4-seeker" ]
