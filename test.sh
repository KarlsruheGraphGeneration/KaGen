#!/bin/bash

echo "test generators"
mkdir test

echo "gnm directed"
mpirun -n 4 --oversubscribe ./build/app/generate_kagen -gen gnm_directed -n 16 -m 20 -k 4 -i 1 -seed 26 -output test/test_gnm_directed_even
mpirun -n 5 --oversubscribe ./build/app/generate_kagen -gen gnm_directed -n 16 -m 20 -k 5 -i 1 -seed 26 -output test/test_gnm_directed_odd

echo "gnm undirected"
mpirun -n 4 --oversubscribe ./build/app/generate_kagen -gen gnm_undirected -n 16 -m 20 -k 4 -i 1 -seed 26 -output test/test_gnm_undirected_even
mpirun -n 5 --oversubscribe ./build/app/generate_kagen -gen gnm_undirected -n 16 -m 20 -k 5 -i 1 -seed 26 -output test/test_gnm_undirected_odd

echo "rgg 2d"
mpirun -n 4 --oversubscribe ./build/app/generate_kagen -gen rgg_2d -n 16 -r 0.0072 -k 4 -i 1 -seed 26 -output test/test_rgg_2d_even
mpirun -n 5 --oversubscribe ./build/app/generate_kagen -gen rgg_2d -n 16 -r 0.0072 -k 5 -i 1 -seed 26 -output test/test_rgg_2d_odd

echo "rgg 3d"
mpirun -n 4 --oversubscribe ./build/app/generate_kagen -gen rgg_3d -n 16 -r 0.01 -k 4 -i 1 -seed 26 -output test/test_rgg_3d_even
mpirun -n 5 --oversubscribe ./build/app/generate_kagen -gen rgg_3d -n 16 -r 0.01 -k 5 -i 1 -seed 26 -output test/test_rgg_3d_odd

echo "rdg 2d"
mpirun -n 4 --oversubscribe ./build/app/generate_kagen -gen rdg_2d -n 12 -k 4 -i 1 -seed 26 -output test/test_rdg_2d_even
mpirun -n 5 --oversubscribe ./build/app/generate_kagen -gen rdg_2d -n 12 -k 5 -i 1 -seed 26 -output test/test_rdg_2d_odd

echo "rdg 3d"
mpirun -n 4 --oversubscribe ./build/app/generate_kagen -gen rdg_3d -n 12 -k 4 -i 1 -seed 26 -output test/test_rdg_3d_even
mpirun -n 5 --oversubscribe ./build/app/generate_kagen -gen rdg_3d -n 12 -k 5 -i 1 -seed 26 -output test/test_rdg_3d_odd

echo "rhg"
mpirun -n 4 --oversubscribe ./build/app/generate_kagen -gen rhg -n 16 -d 16 -gamma 3 -k 4 -i 1 -seed 26 -output test/test_rhg_even
mpirun -n 5 --oversubscribe ./build/app/generate_kagen -gen rhg -n 16 -d 16 -gamma 3 -k 5 -i 1 -seed 26 -output test/test_rhg_odd

echo "ba"
mpirun -n 4 --oversubscribe ./build/app/generate_kagen -gen ba -n 16 -md 16 -k 4 -i 1 -seed 26 -output test/test_ba_even
mpirun -n 5 --oversubscribe ./build/app/generate_kagen -gen ba -n 16 -md 16 -k 5 -i 1 -seed 26 -output test/test_ba_odd
