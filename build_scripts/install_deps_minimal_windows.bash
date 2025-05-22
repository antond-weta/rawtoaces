#!/usr/bin/env bash

# use with
# cmake -D RTA_ENABLE_EIGEN=0 \
#       -D RTA_ENABLE_CERES=0 \
#       -D RTA_BUILD_V1=0

set -ex

vcpkg install \
    openimageio:x64-windows \
    nlohmann-json:x64-windows
