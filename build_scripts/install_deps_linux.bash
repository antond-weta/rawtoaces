#!/usr/bin/env bash

set -ex

time sudo apt-get update

apt-cache search openexr

time sudo apt-get -q -f install -y \
    libunwind-dev libopenexr-dev \
    libboost-dev libboost-filesystem-dev \
    libboost-test-dev \
    libraw-dev libceres-dev

sudo apt-get -q -f remove -y libopenexr-dev

git clone https://github.com/AcademySoftwareFoundation/openexr.git
cd openexr
git checkout v3.2.1

cmake -S . -B build
cmake --build build
sudo cmake --install build
