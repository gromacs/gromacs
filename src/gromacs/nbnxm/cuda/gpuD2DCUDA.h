/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019, by the GROMACS development team, led by
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


#ifndef GMX_D2D_H
#define GMX_D2D_H

#include "config.h"

#include "gromacs/domdec/domdec_constraints.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/domdec_specatomcomm.h"



typedef struct
{
    float xshift[3];  // X shift values
    float fshift[3];  // F shift values

    int   mapSize;    //size of map
    int   sendRank;   //rank to send data to
    int   recvRank;   //rank to recv data from
    int   bPBC;       // Period Boundary Conditions for this rank
    int   recvOffset; //offset for data recieved by this rank

} gpuD2Dparams_t;


GPU_FUNC_QUALIFIER
void dd_setup_gpu_d2d(gmx_domdec_t *dd, matrix box) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void dd_gpu_d2d_add_fshift_to_params(gmx_domdec_t *dd, rvec *fshift) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void gpuD2DSetCommrec(const t_commrec *cr) GPU_FUNC_TERM



#if defined(__CUDACC__)
#include "nbnxm_cuda.h"
#include "nbnxm_cuda_types.h"
#include "cuda.h"
#include "cuda_runtime.h"

GPU_FUNC_QUALIFIER
void gpuBufferOpsXCommCoord(rvec* x_d_ptr, cudaStream_t stream) GPU_FUNC_TERM

GPU_FUNC_QUALIFIER
void gpuBufferOpsFCommCoord(rvec* f_d_ptr, cudaStream_t stream) GPU_FUNC_TERM

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                          \
        cudaError_t e = cudaGetLastError();                                 \
        if (e != cudaSuccess) {                                              \
            printf("Cuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e));           \
            exit(0); \
        }                                                                 \
}


#endif


#endif
