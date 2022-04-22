#!/bin/bash

BINARY=build/app/kagen
PROCS="1 4"
NS="5 10 15"
KAGEN_FLAGS="-V -f none" 
MPI_FLAGS="-oversubscribe"

for PROC in $PROCS; do
	for N in $NS; do 
		echo "##############################"
		echo "### PROC=$PROC N=$N        ###"
		echo "##############################"
		mpirun $MPI_FLAGS -n $PROC $BINARY $KAGEN_FLAGS gnm_undirected -N $N -M $((N + 2)) || exit 1
		mpirun $MPI_FLAGS -n $PROC $BINARY $KAGEN_FLAGS gnp_undirected -N $N -p 0.01 || exit 1
		mpirun $MPI_FLAGS -n $PROC $BINARY $KAGEN_FLAGS rgg_2d -N $N -r 0.125 || exit 1 
		mpirun $MPI_FLAGS -n $PROC $BINARY $KAGEN_FLAGS grid_2d -X $((N / 2)) -Y $((N / 2)) || exit 1 
		mpirun $MPI_FLAGS -n $PROC $BINARY $KAGEN_FLAGS rhg -N $N -g 3.0 -d 8 || exit 1
	done
done
