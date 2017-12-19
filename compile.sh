#!/bin/bash

mkdir build
cd build
# cmake -DCMAKE_C_COMPILER=gcc-7 -DCMAKE_CXX_COMPILER=g++-7 -DCMAKE_PREFIX_PATH=/home/hpc/pr87si/di36mek/usr ../
cmake -DCMAKE_C_COMPILER=gcc-7 -DCMAKE_CXX_COMPILER=g++-7 ../
make -j1
