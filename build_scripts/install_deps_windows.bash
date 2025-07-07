#!/usr/bin/env bash

set -ex

vcpkg install \
    ceres:x64-windows \
    openimageio:x64-windows \
    nlohmann-json:x64-windows
    
curl --silent --output ./exiftool.zip https://exiftool.org/exiftool-13.32_64.zip
unzip ./exiftool.zip
mv ./exiftool-13.32_64/exiftool\(-k\).exe ./exiftool-13.32_64/exiftool.exe
