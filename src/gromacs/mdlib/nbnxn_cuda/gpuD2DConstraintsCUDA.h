/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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


#ifndef GMX_D2D_CONSTRAINTS_H
#define GMX_D2D_CONSTRAINTS_H


#include "gromacs/domdec/domdec_constraints.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/domdec_specatomcomm.h"



typedef struct
{

    //Forward direction
    float shiftFW[3];   //shift values
    int   mapSizeFW;    //size of map
    int   sendRankFW;   //rank to send data to
    int   recvRankFW;   //rank to recv data from
    int   bPBCFW;       // Period Boundary Conditions for this rank
    int   recvOffsetFW; //offset for data recieved by this rank


    //Backward Direction direction
    float shiftBW[3];   //shift values
    int   mapSizeBW;    //size of map
    int   sendRankBW;   //rank to send data to
    int   recvRankBW;   //rank to recv data from
    int   bPBCBW;       // Period Boundary Conditions for this rank
    int   recvOffsetBW; //offset for data recieved by this rank

} gpuConstraintsD2Dparams_t;

enum dirFWBW {
    DIRFW,
    DIRBW
};


#if defined(__CUDACC__)
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_types.h"
#include "cuda.h"
#include "cuda_runtime.h"



GPU_FUNC_QUALIFIER
void gpuConstraintsD2D(rvec* x_d_ptr, const t_commrec *cr, cudaStream_t stream, matrix box, bool setupRequired) GPU_FUNC_TERM


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
