#!/usr/bin/env bash

set -ex

time sudo apt-get update

time sudo apt-get -q -f install -y \
    libboost-dev libboost-filesystem-dev \
    libboost-test-dev \
    libraw-dev libceres-dev
