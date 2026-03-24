/*
 * Pseudo-MPI header for single-process (NOMPI) builds of KaGen.
 *
 * This header provides ONLY type definitions and constants — enough to let
 * kagen.h compile without a real MPI installation.  All actual collective
 * semantics are handled by kagen::SeqComm; no function stubs are needed here.
 *
 * Inspired by: https://github.com/KaHIP/VieClus/blob/master/lib/tools/pseudo_mpi.h
 */
#ifndef PSEUDO_MPI_H
#define PSEUDO_MPI_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================== */
/* Types                                                                       */
/* ========================================================================== */

typedef int MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Op;
typedef int MPI_Request;

typedef struct {
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
} MPI_Status;

/* ========================================================================== */
/* Constants                                                                   */
/* ========================================================================== */

#define MPI_COMM_WORLD  0
#define MPI_COMM_SELF   0
#define MPI_SUCCESS     0

/* MPI_IN_PLACE sentinel -- must be a unique pointer value */
#define MPI_IN_PLACE    ((void*)1)

#define MPI_DATATYPE_NULL   0

/* Datatype constants -- values match pseudo_mpi_sizeof() dispatch in the old
 * stub for any downstream code that still compares against them. */
#define MPI_C_BOOL              1
#define MPI_INT                 2
#define MPI_UNSIGNED            3
#define MPI_LONG_LONG           4
#define MPI_UNSIGNED_LONG_LONG  5
#define MPI_DOUBLE              6
#define MPI_LONG_DOUBLE         7
#define MPI_BYTE                8
#define MPI_UINT64_T            9
#define MPI_INT64_T             10

/* Op constants */
#define MPI_SUM  1
#define MPI_MAX  2
#define MPI_MIN  3
#define MPI_LOR  4

/* Misc */
#define MPI_STATUS_IGNORE   ((MPI_Status*)0)
#define MPI_STATUSES_IGNORE ((MPI_Status*)0)

#ifdef __cplusplus
}
#endif

#endif /* PSEUDO_MPI_H */
