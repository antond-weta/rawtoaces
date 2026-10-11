#!/usr/bin/env bash

set -ex

brew update

brew install \
    nlohmann-json \
    openimageio \
    nanobind \
    robin-map \
    exiftool \
    pytest

if [ "$RTA_ENABLE_EIGEN" != "OFF" ];
then
    brew install eigen
fi

if [ "$RTA_ENABLE_CERES" != "OFF" ];
then
    brew install ceres-solver
fi

if [ "$RTA_ENABLE_LENSFUN" != "OFF" ];
then
    brew install lensfun
fi
