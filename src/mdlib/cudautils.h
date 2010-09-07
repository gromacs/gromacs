#ifndef CUDAUTILS_H
#define CUDAUTILS_H

#include "stdio.h"

#include "gmx_fatal.h"

#include "cuda.h"

#define CU_RET_ERR(status, s) \
    do { \
        if (status != cudaSuccess) \
        { \
            gmx_fatal(FARGS, "%s: %s\n", s, cudaGetErrorString(status)); \
        } \
    } while (0)


#define CU_LAUNCH_ERR(s) \
    do { \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            gmx_fatal(FARGS, "Error while launching kernel %s: %s\n", s, cudaGetErrorString(status)); \
        } \
    } while (0)


#define GRID_MAX_DIM        65535
#define NB_DEFAULT_THREADS  256


#ifdef __cplusplus
extern "C" {
#endif
int download_cudata(void * /*h_dest*/, void * /*d_src*/, size_t /*bytes*/);

int download_cudata_alloc(void ** /*h_dest*/, void * /*d_src*/, size_t /*bytes*/);

int upload_cudata(void * /*d_dest*/, void * /*h_src*/, size_t /*bytes*/);

int upload_cudata_alloc(void ** /*d_dest*/, void * /*h_src*/, size_t /*bytes*/);

void * pmalloc(size_t /*bytes*/, FILE * /*fplog*/);

void pfree(void * /*h_ptr*/, FILE * /*fplog*/);

#ifdef __cplusplus
}
#endif


#endif
