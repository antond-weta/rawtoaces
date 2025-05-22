#!/usr/bin/env bash

# use with
# cmake -D RTA_ENABLE_EIGEN=0 \
#       -D RTA_ENABLE_CERES=0 \
#       -D RTA_BUILD_V1=0

set -ex

brew install openimageio nlohmann-json
