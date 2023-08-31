FROM ubuntu:18.04

MAINTAINER Eugene Chow "eugene.chow@link.cuhk.edu.hk'

RUN apt-get update -y && apt-get install -y python3-pip python3-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

RUN apt-get install -y curl
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python3 get-pip.py
RUN mkdir /opt/rG4-seeker-manuscript/
COPY requirements.txt /opt/rG4-seeker-manuscript/requirements.txt
COPY setup.py /opt/rG4-seeker-manuscript/setup.py
COPY rg4seeker/ /opt/rG4-seeker-manuscript/rg4seeker/
WORKDIR /opt/rG4-seeker-manuscript/
RUN pip install .

ENTRYPOINT [ "rG4-seeker" ]
