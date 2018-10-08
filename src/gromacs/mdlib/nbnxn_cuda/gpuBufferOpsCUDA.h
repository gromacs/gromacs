/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_MDLIB_NBNXN_BUFFER_OPS_H
#define GMX_MDLIB_NBNXN_BUFFER_OPS_H


#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-grid.h"
#include "gromacs/ewald/pme-internal.h"
#include "gromacs/mdlib/nbnxn_atomdata.h"
#include "gromacs/mdlib/nbnxn_gpu_types.h"

#include "cuda.h"
#include "cuda_runtime.h"

#include "nbnxn_cuda.h"
#include "nbnxn_cuda_types.h"

#define TPB 128 //CUDA threads per block

//to do fix this
#define BO_MAX_RANKS 32

GPU_FUNC_QUALIFIER
void gpuBufferOpsCopyRvecToNbatReal(int ncxy, int g, int FillLocal,
                                    gmx_pme_t* pmedata, nbnxn_atomdata_t *nbat,
                                    gmx_nbnxn_gpu_t *gpu_nbv, const int* a,
                                    int a_nalloc, int* na_all, int* cxy_ind,
                                    int cell0, int na_sc, int iloc, int stride,
                                    const float x[][3], int xsize, t_commrec *cr) GPU_FUNC_TERM

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
gmx_bool gpuBufferOpsTimestepInitFromPP(gmx_bool   bNS,
                                        gmx_bool   bUseGPU,
                                        gmx_bool   bDutyPPAndPME,
                                        t_commrec *cr,
                                        int        natoms_all) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
gmx_bool gpuBufferOpsTimestepInitFromPME(gmx_bool   bNS,
                                         gmx_bool   bUseGPU,
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
