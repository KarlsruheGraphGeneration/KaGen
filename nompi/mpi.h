/*
 * Pseudo-MPI header for single-process (NOMPI) builds of KaGen.
 *
 * This header provides stub implementations of all MPI functions used by KaGen,
 * allowing the library to compile and run without a real MPI installation.
 * All collective operations are either no-ops or simple memcpy for the single-PE case.
 *
 * Inspired by: https://github.com/KaHIP/VieClus/blob/master/lib/tools/pseudo_mpi.h
 */
#ifndef PSEUDO_MPI_H
#define PSEUDO_MPI_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

/* Datatype constants -- unique IDs for pseudo_mpi_sizeof() dispatch */
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

/* ========================================================================== */
/* Helper: datatype -> byte size                                               */
/* ========================================================================== */

/*
 * For derived types created via MPI_Type_contiguous, we encode the total byte
 * size as (size + 1024). Built-in types use IDs 0-10. This lets the helper
 * decode both built-in and derived types from a single integer.
 */
static inline size_t pseudo_mpi_sizeof(MPI_Datatype type) {
    if (type >= 1024) {
        return (size_t)(type - 1024);
    }
    switch (type) {
        case MPI_C_BOOL:              return sizeof(_Bool);
        case MPI_INT:                 return sizeof(int);
        case MPI_UNSIGNED:            return sizeof(unsigned int);
        case MPI_LONG_LONG:           return sizeof(long long);
        case MPI_UNSIGNED_LONG_LONG:  return sizeof(unsigned long long);
        case MPI_DOUBLE:              return sizeof(double);
        case MPI_LONG_DOUBLE:         return sizeof(long double);
        case MPI_BYTE:                return 1;
        case MPI_UINT64_T:            return sizeof(uint64_t);
        case MPI_INT64_T:             return sizeof(int64_t);
        case MPI_DATATYPE_NULL:       return 0;
        default:                      return 0;
    }
}

/* ========================================================================== */
/* Initialization / Finalization                                               */
/* ========================================================================== */

static inline int MPI_Init(int* argc, char*** argv) {
    (void)argc; (void)argv;
    return MPI_SUCCESS;
}

static inline int MPI_Finalize(void) {
    return MPI_SUCCESS;
}

/* ========================================================================== */
/* Communicator queries                                                        */
/* ========================================================================== */

static inline int MPI_Comm_rank(MPI_Comm comm, int* rank) {
    (void)comm;
    *rank = 0;
    return MPI_SUCCESS;
}

static inline int MPI_Comm_size(MPI_Comm comm, int* size) {
    (void)comm;
    *size = 1;
    return MPI_SUCCESS;
}

/* ========================================================================== */
/* Synchronization                                                             */
/* ========================================================================== */

static inline int MPI_Barrier(MPI_Comm comm) {
    (void)comm;
    return MPI_SUCCESS;
}

/* ========================================================================== */
/* Timing                                                                      */
/* ========================================================================== */

static inline double MPI_Wtime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Error handling                                                              */
/* ========================================================================== */

static inline int MPI_Abort(MPI_Comm comm, int errorcode) {
    (void)comm;
    exit(errorcode);
    return MPI_SUCCESS; /* unreachable */
}

/* ========================================================================== */
/* Collective operations                                                       */
/* ========================================================================== */

/*
 * MPI_Allreduce: For a single PE, the result is the input itself.
 * If sendbuf == MPI_IN_PLACE, data is already in recvbuf -- no-op.
 * Otherwise, copy sendbuf to recvbuf.
 */
static inline int MPI_Allreduce(
    const void* sendbuf, void* recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    (void)op; (void)comm;
    if (sendbuf != MPI_IN_PLACE && sendbuf != recvbuf) {
        memcpy(recvbuf, sendbuf, (size_t)count * pseudo_mpi_sizeof(datatype));
    }
    return MPI_SUCCESS;
}

/*
 * MPI_Reduce: With 1 PE, root is always 0 (us). Same as Allreduce.
 */
static inline int MPI_Reduce(
    const void* sendbuf, void* recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
    (void)op; (void)root; (void)comm;
    if (sendbuf != MPI_IN_PLACE && sendbuf != recvbuf) {
        memcpy(recvbuf, sendbuf, (size_t)count * pseudo_mpi_sizeof(datatype));
    }
    return MPI_SUCCESS;
}

/*
 * MPI_Bcast: With 1 PE, data is already at the single process -- no-op.
 */
static inline int MPI_Bcast(
    void* buffer, int count, MPI_Datatype datatype,
    int root, MPI_Comm comm)
{
    (void)buffer; (void)count; (void)datatype; (void)root; (void)comm;
    return MPI_SUCCESS;
}

/*
 * MPI_Allgather: For 1 PE, if MPI_IN_PLACE the data is already in recvbuf.
 * Otherwise copy sendbuf to recvbuf.
 */
static inline int MPI_Allgather(
    const void* sendbuf, int sendcount, MPI_Datatype sendtype,
    void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    (void)comm;
    if (sendbuf != MPI_IN_PLACE) {
        memcpy(recvbuf, sendbuf, (size_t)sendcount * pseudo_mpi_sizeof(sendtype));
    }
    /* If MPI_IN_PLACE: data already at position [rank=0] in recvbuf */
    (void)recvcount; (void)recvtype;
    return MPI_SUCCESS;
}

/*
 * MPI_Allgatherv: For 1 PE, copy sendbuf to recvbuf at displacement[0].
 */
static inline int MPI_Allgatherv(
    const void* sendbuf, int sendcount, MPI_Datatype sendtype,
    void* recvbuf, const int* recvcounts, const int* displs,
    MPI_Datatype recvtype, MPI_Comm comm)
{
    (void)comm; (void)recvcounts;
    if (sendbuf != MPI_IN_PLACE) {
        size_t elem_size = pseudo_mpi_sizeof(recvtype);
        char* dst = (char*)recvbuf + (size_t)displs[0] * elem_size;
        memcpy(dst, sendbuf, (size_t)sendcount * pseudo_mpi_sizeof(sendtype));
    }
    return MPI_SUCCESS;
}

/*
 * MPI_Alltoall: For 1 PE, each PE sends sendcount elements to each PE.
 * With size=1, this is just a memcpy of sendcount elements.
 */
static inline int MPI_Alltoall(
    const void* sendbuf, int sendcount, MPI_Datatype sendtype,
    void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    (void)recvcount; (void)recvtype; (void)comm;
    memcpy(recvbuf, sendbuf, (size_t)sendcount * pseudo_mpi_sizeof(sendtype));
    return MPI_SUCCESS;
}

/*
 * MPI_Alltoallv: For 1 PE, copy from sendbuf[sdispls[0]] to recvbuf[rdispls[0]]
 * using sendcounts[0] elements.
 */
static inline int MPI_Alltoallv(
    const void* sendbuf, const int* sendcounts, const int* sdispls, MPI_Datatype sendtype,
    void* recvbuf, const int* recvcounts, const int* rdispls, MPI_Datatype recvtype,
    MPI_Comm comm)
{
    (void)recvcounts; (void)recvtype; (void)comm;
    size_t elem_size = pseudo_mpi_sizeof(sendtype);
    const char* src = (const char*)sendbuf + (size_t)sdispls[0] * elem_size;
    char* dst = (char*)recvbuf + (size_t)rdispls[0] * elem_size;
    memcpy(dst, src, (size_t)sendcounts[0] * elem_size);
    return MPI_SUCCESS;
}

/*
 * MPI_Exscan: With 1 PE, the result on rank 0 is undefined per the MPI
 * standard. KaGen code always pre-initializes the output to 0, so leaving
 * recvbuf untouched is correct.
 */
static inline int MPI_Exscan(
    const void* sendbuf, void* recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    (void)sendbuf; (void)recvbuf; (void)count;
    (void)datatype; (void)op; (void)comm;
    return MPI_SUCCESS;
}

/* ========================================================================== */
/* Derived datatype operations                                                 */
/* ========================================================================== */

/*
 * MPI_Type_contiguous: Create a derived type whose size is count * base_size.
 * We encode the total byte size as (size + 1024) so pseudo_mpi_sizeof can
 * decode it. Built-in type IDs are 0-10, so there's no collision.
 */
static inline int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype* newtype) {
    *newtype = (int)((size_t)count * pseudo_mpi_sizeof(oldtype)) + 1024;
    return MPI_SUCCESS;
}

static inline int MPI_Type_commit(MPI_Datatype* type) {
    (void)type;
    return MPI_SUCCESS;
}

static inline int MPI_Type_free(MPI_Datatype* type) {
    *type = MPI_DATATYPE_NULL;
    return MPI_SUCCESS;
}

/* ========================================================================== */
/* Point-to-point (stubs -- not used by KaGen in single-PE mode)               */
/* ========================================================================== */

static inline int MPI_Send(
    const void* buf, int count, MPI_Datatype datatype,
    int dest, int tag, MPI_Comm comm)
{
    (void)buf; (void)count; (void)datatype; (void)dest; (void)tag; (void)comm;
    return MPI_SUCCESS;
}

static inline int MPI_Recv(
    void* buf, int count, MPI_Datatype datatype,
    int source, int tag, MPI_Comm comm, MPI_Status* status)
{
    (void)buf; (void)count; (void)datatype; (void)source; (void)tag; (void)comm;
    if (status) {
        status->MPI_SOURCE = 0;
        status->MPI_TAG = tag;
        status->MPI_ERROR = MPI_SUCCESS;
    }
    return MPI_SUCCESS;
}

static inline int MPI_Gather(
    const void* sendbuf, int sendcount, MPI_Datatype sendtype,
    void* recvbuf, int recvcount, MPI_Datatype recvtype,
    int root, MPI_Comm comm)
{
    (void)recvcount; (void)recvtype; (void)root; (void)comm;
    if (sendbuf != MPI_IN_PLACE) {
        memcpy(recvbuf, sendbuf, (size_t)sendcount * pseudo_mpi_sizeof(sendtype));
    }
    return MPI_SUCCESS;
}

#ifdef __cplusplus
}
#endif

#endif /* PSEUDO_MPI_H */
