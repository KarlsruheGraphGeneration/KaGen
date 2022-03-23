#!/bin/bash

# Init submodules
git submodule update --init --recursive

# Compile code
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make -j4

