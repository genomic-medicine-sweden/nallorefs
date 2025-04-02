# syntax=docker/dockerfile:1
FROM ubuntu:22.04

USER root

RUN apt-get update && apt-get install -y procps sed && rm -rf /var/lib/apt/lists/*
# install app
COPY vcfexpress /usr/bin/vcfexpress

# final configuration
ENV FLASK_APP=vcfexpress
