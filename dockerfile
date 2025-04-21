FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends r-base python3 python3-pip
RUN python3 -m pip install matplotlib==3.7.1
RUN python3 -m pip install numpy===1.24.3
RUN R -e "install.packages('here')"
WORKDIR /comola
COPY . .
CMD python3 __init__.py