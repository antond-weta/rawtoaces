#!/usr/bin/env bash

set -ex

git clone https://github.com/ampas/aces_container.git v1/aces_container

if [[ "$OSTYPE" == "linux-gnu"* ]] || [[ "$OSTYPE" == "darwin"* ]]; then
    cmake \
        -S v1/aces_container \
        -B v1/aces_container/build \
        -DCMAKE_CXX_FLAGS="-Wno-c++11-narrowing"
    cmake --build v1/aces_container/build
    sudo cmake --install v1/aces_container/build
else
    cmake \
        -S v1/aces_container \
        -B v1/aces_container/build \
        -DCMAKE_INSTALL_PREFIX="." \
        -DBUILD_SHARED_LIBS=OFF
    cmake --build v1/aces_container/build --config Release
    cmake --install v1/aces_container/build --config Release
fi

cd ../..
