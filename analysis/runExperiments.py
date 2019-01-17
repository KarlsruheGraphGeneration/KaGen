#!/usr/bin/env python3

import glob
import os
import subprocess
import itertools
from multiprocessing import Pool

NODES = [10, 12, 14, 16, 18, 20, 22, 24]
CHUNKS_2D = [1, 4, 9, 16, 36, 64, 121, 256, 529, 1024, 2025, 3025, 4225, 5184, 6084, 7056, 8100, 9025, 10000]
CHUNKS_3D = [1, 8, 27, 64, 125, 216, 512, 1000, 2197, 4096, 5832, 6859, 8000, 9261, 10648]

def runExperiment(params):

    d = params[0]
    n = params[1]
    k = params[2]

    myPID = os.getpid()
    filename = "delaunay_%id_%i_%i" % (d, n, k)

    call = "mpirun -n 1 ../build/app/kagen -gen rdg_" + str(d) + "d -n " + str(n) + " -seed 1337 -k " + str(k) + " -i 10 -output " + filename

    print(myPID, call)
    if subprocess.call(call, shell=True):
        print("ERROR:", call)


if __name__ == '__main__':
    myPID = os.getpid()

    with Pool(os.cpu_count()) as p:
        p.map(runExperiment, itertools.product([2], NODES, CHUNKS_2D))
        p.map(runExperiment, itertools.product([3], NODES, CHUNKS_3D))
