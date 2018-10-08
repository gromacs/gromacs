

#ifndef GMX_MDLIB_NBNXN_BUFFER_OPS_H
#define GMX_MDLIB_NBNXN_BUFFER_OPS_H


#include "gromacs/mdlib/nbnxn_gpu_types.h"
#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-grid.h"

#include "gromacs/mdlib/nbnxn_atomdata.h"

#define TPB 128 //CUDA threads per block

//to do fix this
#define BO_MAX_RANKS 32

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                                            \
    cudaError_t e=cudaGetLastError();                                                 \
    if(e!=cudaSuccess)                                                                \
    {                                                                                 \
        printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
        exit(0);                                                                      \
    }                                                                                 \
}

GPU_FUNC_QUALIFIER
void gpuBufferOpsCopyRvecToNbatReal(int ncxy, int g, int FillLocal,
                                    gmx_pme_t* pmedata, nbnxn_atomdata_t *nbat,
                                    gmx_nbnxn_gpu_t *gpu_nbv, const int* a,
                                    int a_nalloc, int* na_all, int* cxy_ind,
                                    int cell0, int na_sc, int iloc, int stride,
                                    const float x[][3], int xsize,t_commrec *cr) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void gpuBufferOpsAddNbatFToF(const int* cell,
                             const nbnxn_atomdata_t *nbat,
                             gmx_nbnxn_gpu_t *gpu_nbv,
                             nbnxn_atomdata_output_t *out,
                             const gmx_pme_t* pmedata,
                             int nfa,
                             int a0, int a1,
                             rvec *f) GPU_FUNC_TERM
  
GPU_FUNC_QUALIFIER
void gpuBufferOpsPMEGatherFromPP(float* dest, int* offset, int* size, int nsender) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void gpuBufferOpsPMEDistributeToPP(float* src, int* offset, int* copysize, int nsender) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void gpuBufferOpsPPRecvFromPME(float* src, int offset, t_commrec *cr) GPU_FUNC_TERM


GPU_FUNC_QUALIFIER
void gpuBufferOpsInitXSendBufIndexMap(int** mapptr_in, float** shiftptr_in, int size) GPU_FUNC_TERM

  GPU_FUNC_QUALIFIER
void gpuBufferOpsInitFRecvBufIndexMap(int** mapptr_in, float** shiftptr_in, int size) GPU_FUNC_TERM


GPU_FUNC_QUALIFIER
gmx_bool gpuBufferOpsTimestepInitFromPP(gmx_bool bNS,
                                        gmx_bool bUseGPU,
                                        gmx_bool bDutyPPAndPME,
                                        t_commrec *cr,
                                        int natoms_all) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
gmx_bool gpuBufferOpsTimestepInitFromPME(gmx_bool bNS,
                                         gmx_bool bUseGPU,
                                         t_commrec *cr) GPU_FUNC_TERM
    
GPU_FUNC_QUALIFIER
gmx_bool gpuBufferOpsActiveThisTimestep() GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
rvec* gpuBufferOpsGetXPtr() GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
rvec* gpuBufferOpsGetFPtr() GPU_FUNC_TERM
    
#if defined(__CUDACC__)
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_types.h"
#include "cuda.h"
#include "cuda_runtime.h"

void gpuBufferOpsFlagDataReady(cudaStream_t streamIn);


#endif

#endif
