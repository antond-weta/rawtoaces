#!/usr/bin/env bash

set -ex

vcpkg install \
    pkgconf:x64-windows

vcpkg install \
    ceres:x64-windows \
    openimageio:x64-windows \
    nlohmann-json:x64-windows \
    lensfun:x64-windows \
    exiftool
    
which exiftool
