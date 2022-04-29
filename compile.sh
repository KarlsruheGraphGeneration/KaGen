#!/bin/bash

# Init submodules
git submodule update --init --recursive

# Compile code
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make kagen kagen_library -j4

