FROM python:3.9-buster AS fitseq

LABEL version=1.5.0-rc

ENV DEBIAN_FRONTEND="noninteractive"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV LANGUAGE=C.UTF-8
RUN python3 -m pip install git+https://github.com/darachm/fitseq.git
