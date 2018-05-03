#!/bin/bash

# Init submodules
git submodule update --init --recursive

# Compile code
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=gcc-7 -DCMAKE_CXX_COMPILER=g++-7 ../
make -j4

