#!/usr/bin/env bash

set -ex

time sudo apt-get update

time sudo apt-get -q -f install -y \
    libceres-dev \
    libopencv-dev \
    liblensfun-dev \
    liblensfun-data-v1 \
    openimageio-tools libopenimageio-dev \
    nlohmann-json3-dev \
    exiftool
