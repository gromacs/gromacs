/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "gromacs/domdec/domdec.h"

#include "config.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include <algorithm>
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/math/vec.h"
#include "gromacs/nbnxm/cuda/gpuUpdateConstraintsCUDA.h"
#include "gromacs/nbnxm/cuda/gpuD2DCUDA.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/domdec/domdec_constraints.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/domdec_specatomcomm.h"


struct gpuD2DConstraintsData{

/* pointers for inter-GPU communication index maps and data buffers */
    int                      * sendBufFWIndexMap;
    int                      * sendBufFWIndexMap_d;
    int                        sendBufFWIndexMapSize;
    float                    * sendBufFW_d;
    int                      * sendBufBWIndexMap;
    int                      * sendBufBWIndexMap_d;
    int                        sendBufBWIndexMapSize;
    float                    * sendBufBW_d;

    gpuConstraintsD2Dparams_t* d2dConstraintsParams;
    gpuConstraintsD2Dparams_t* d2dConstraintsParams_d;


};


//TEMPORARY Global variable to store a copy of above struct within this module.
//TODO move this to a central location (e.g. gpu_nbv) and pass through fn args.
static gpuD2DConstraintsData gpuD2DConsMod;

/*-------------------------------- CUDA kernels --------------------------------*/


/* pack non-local coordinate data buffer on the GPU using pre-populated "map" containing index information */
__global__ void pack_send_buf_kernel_constraints(float                    * data_packed,
                                                 float                    * data,
                                                 const int                * map,
                                                 gpuConstraintsD2Dparams_t* d2dParams,
                                                 int                        dir)
{

    int   idx = blockIdx.x*blockDim.x+threadIdx.x;

    int   data_size;
    int   bPBC;

    float shiftX;
    float shiftY;
    float shiftZ;

    if (dir == DIRFW)
    {
        data_size = d2dParams->mapSizeFW;
        bPBC      = d2dParams->bPBCFW;

        shiftX = d2dParams->shiftFW[0];
        shiftY = d2dParams->shiftFW[1];
        shiftZ = d2dParams->shiftFW[2];
    }
    else
    {
        data_size = d2dParams->mapSizeBW;
        bPBC      = d2dParams->bPBCBW;

        shiftX = d2dParams->shiftBW[0];
        shiftY = d2dParams->shiftBW[1];
        shiftZ = d2dParams->shiftBW[2];
    }

    if (idx < data_size)
    {

        if (bPBC)
        {
            data_packed[idx*3]   = data[3*map[idx]]+shiftX;
            data_packed[idx*3+1] = data[3*map[idx]+1]+shiftY;
            data_packed[idx*3+2] = data[3*map[idx]+2]+shiftZ;
        }
        else
        {
            data_packed[idx*3]   = data[3*map[idx]];
            data_packed[idx*3+1] = data[3*map[idx]+1];
            data_packed[idx*3+2] = data[3*map[idx]+2];
        }
    }
    return;
}

/*-------------------------------- End CUDA kernels-----------------------------*/


// Initialize map structure to allow gpu buffer packing for inter-GPU "move_x" style communication
static void gpuConstraintsInitD2D(int                      ** mapptrFW_in,
                                  int                         FWsize,
                                  int                      ** mapptrBW_in,
                                  int                         BWsize,
                                  gpuConstraintsD2Dparams_t** d2dParams_in)
{

    gpuD2DConstraintsData* gpuD2D = &gpuD2DConsMod;

    // Forward direction
    if (FWsize > gpuD2D->sendBufFWIndexMapSize)
    {

        gpuD2D->sendBufFWIndexMapSize = FWsize;

        if (gpuD2D->sendBufFWIndexMap)
        {
            free(gpuD2D->sendBufFWIndexMap);
            cudaFree(gpuD2D->sendBufFWIndexMap_d);
        }
        gpuD2D->sendBufFWIndexMap = (int*) malloc(FWsize*sizeof(int));
        cudaMalloc(&gpuD2D->sendBufFWIndexMap_d, FWsize*sizeof(int));
        cudaMalloc(&gpuD2D->sendBufFW_d, FWsize*sizeof(rvec));

        if (gpuD2D->d2dConstraintsParams)
        {
            free(gpuD2D->d2dConstraintsParams);
            cudaFree(gpuD2D->d2dConstraintsParams_d);
        }

        gpuD2D->d2dConstraintsParams = (gpuConstraintsD2Dparams_t*) malloc(sizeof(gpuConstraintsD2Dparams_t));
        cudaMalloc(&gpuD2D->d2dConstraintsParams_d, sizeof(gpuConstraintsD2Dparams_t));

    }

    *mapptrFW_in  = gpuD2D->sendBufFWIndexMap;
    *d2dParams_in = gpuD2D->d2dConstraintsParams;



    // Backward direction
    if (BWsize > gpuD2D->sendBufBWIndexMapSize)
    {
        gpuD2D->sendBufBWIndexMapSize = BWsize;

        if (gpuD2D->sendBufBWIndexMap)
        {
            free(gpuD2D->sendBufBWIndexMap);
            cudaFree(gpuD2D->sendBufBWIndexMap_d);
        }
        gpuD2D->sendBufBWIndexMap = (int*) malloc(BWsize*sizeof(int));
        cudaMalloc(&gpuD2D->sendBufBWIndexMap_d, BWsize*sizeof(int));
        cudaMalloc(&gpuD2D->sendBufBW_d, BWsize*sizeof(rvec));



    }

    *mapptrBW_in = gpuD2D->sendBufBWIndexMap;

    return;
}


// This funcion is based on move_x_specat, but instead of actually sending and recieving
// data a buffer is packed with index information on indexes to be sent/recieved.
// This buffer will later be passed to the GPU and used to pack a buffer, which will
// be sent directly between GPUs.
static void dd_setup_gpu_d2d_constraints(gmx_domdec_t *dd, gmx_domdec_specat_comm_t *spac,
                                         const matrix box)
{
    gmx_specatsend_t *spas;
    int               d, dim, dir;
    size_t            i;
    gmx_bool          bPBC, bScrew = FALSE;
    rvec              shift = {0, 0, 0};



    gpuD2DConstraintsData* gpuD2D = &gpuD2DConsMod;


    if (dd->ndim > 1)
    {
        gmx_fatal(FARGS, "Error: dd->ndim > 1 is not yet supported in GPU Constraints D2D");
    }


    //forward

    int             * sendBufFWIndexMap;
    gmx_specatsend_t *spasforward;
    spasforward = &spac->spas[0][1];
    int             * sendBufBWIndexMap;
    gmx_specatsend_t *spasbackward;
    spasbackward = &spac->spas[0][0];


    // Create and populate a data structure containing information on this communication pattern
    gpuConstraintsD2Dparams_t *d2dParams;
    gpuConstraintsInitD2D(&sendBufFWIndexMap, spasforward->a.size(), &sendBufBWIndexMap,
                          spasbackward->a.size(), &d2dParams);


    //size of map
    d2dParams->mapSizeFW = spasforward->a.size();
    //rank to send data to
    d2dParams->sendRankFW = dd->neighbor[0][0];
    //rank to recv data from
    d2dParams->recvRankFW = dd->neighbor[0][1];
    // Period Boundary Conditions for this rank
    d2dParams->bPBCFW = (dd->ci[0] == dd->nc[0]-1);
    //offset for data recieved by this rank
    d2dParams->recvOffsetFW = spac->at_start;


    d2dParams->shiftFW[0] = 0.;
    d2dParams->shiftFW[1] = 0.;
    d2dParams->shiftFW[2] = 0.;


    //backward

    //size of map
    d2dParams->mapSizeBW = spasbackward->a.size();
    //rank to send data to
    d2dParams->sendRankBW = dd->neighbor[0][1];
    //rank to recv data from
    d2dParams->recvRankBW = dd->neighbor[0][0];
    //offset for data recieved by this rank
    d2dParams->recvOffsetBW = spac->at_start+spasforward->nrecv;

    //three components of shift for PBC
    d2dParams->shiftBW[0] = 0.;
    d2dParams->shiftBW[1] = 0.;
    d2dParams->shiftBW[2] = 0.;


    // Populate index map including those indexes to be communicated
    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];
        if (dd->nc[dim] > 2)
        {
            /* Pulse the grid forward and backward */
            for (dir = 0; dir < 2; dir++)
            {
                int idx = 0;

                if (dir == 0 && dd->ci[dim] == 0)
                {
                    bPBC   = TRUE;
                    bScrew = (dd->bScrewPBC && dim == XX);
                    copy_rvec(box[dim], shift);

                }
                else if (dir == 1 && dd->ci[dim] == dd->nc[dim]-1)
                {
                    bPBC   = TRUE;
                    bScrew = (dd->bScrewPBC && dim == XX);
                    for (i = 0; i < DIM; i++)
                    {
                        shift[i] = -box[dim][i];
                    }
                }
                else
                {
                    bPBC   = FALSE;
                    bScrew = FALSE;
                }

                if (dir == 0)
                {

                    // Period Boundary Conditions for this rank
                    d2dParams->bPBCBW = bPBC;

                    //three components of shift for PBC
                    for (i = 0; i < DIM; i++)
                    {
                        d2dParams->shiftBW[i] = shift[i];
                    }
                }

                if (dir == 1)
                {

                    // Period Boundary Conditions for this rank
                    d2dParams->bPBCFW = bPBC;
                    //three components of shift for PBC
                    for (i = 0; i < DIM; i++)
                    {
                        d2dParams->shiftFW[i] = shift[i];
                    }
                }

                spas = &spac->spas[d][dir];
                /* Copy the required coordinates to the send buffer */
                if (!bPBC)
                {
                    /* Only copy */
                    for (i = 0; i < spas->a.size(); i++)
                    {

                        if (dir == 0)
                        {
                            gpuD2D->sendBufBWIndexMap[idx] = spas->a[i];
                            idx++;
                        }
                        if (dir == 1)
                        {
                            gpuD2D->sendBufFWIndexMap[idx] = spas->a[i];
                            idx++;
                        }


                    }
                }
                else if (!bScrew)
                {
                    /* Shift coordinates */
                    for (i = 0; i < spas->a.size(); i++)
                    {

                        if (dir == 0)
                        {
                            gpuD2D->sendBufBWIndexMap[idx] = spas->a[i];
                            idx++;
                        }
                        if (dir == 1)
                        {
                            gpuD2D->sendBufFWIndexMap[idx] = spas->a[i];
                            idx++;
                        }


                    }
                }
                else
                {
                    gmx_fatal(FARGS, "Error: bScrew=true not yet supported in gpu constraints D2D\n");
                }

            }
        }
        else
        {
            gmx_fatal(FARGS, "Error: This path not yet supported in gpu constraints d2d\n");
        }
    }
}




/* function to perform non-local X coordinate communication on GPU */
static void gpuConstraintsXCommCoordSendRecv_D2D(rvec* destbuf, float* buf_s, int dir)
{

    gpuD2DConstraintsData* gpuD2D = &gpuD2DConsMod;

    int                    copysize_s, rank_s, rank_r, offset;
    if (dir == DIRFW)
    {
        copysize_s = gpuD2D->d2dConstraintsParams->mapSizeFW;
        rank_s     = gpuD2D->d2dConstraintsParams->sendRankFW;
        rank_r     = gpuD2D->d2dConstraintsParams->recvRankFW;
        offset     = gpuD2D->d2dConstraintsParams->recvOffsetFW;
    }
    else
    {
        copysize_s = gpuD2D->d2dConstraintsParams->mapSizeBW;
        rank_s     = gpuD2D->d2dConstraintsParams->sendRankBW;
        rank_r     = gpuD2D->d2dConstraintsParams->recvRankBW;
        offset     = gpuD2D->d2dConstraintsParams->recvOffsetBW;
    }


    int copysize_r;
    //get remote copysize
    MPI_Sendrecv(&copysize_s, sizeof(copysize_s), MPI_BYTE, rank_s, 0,
                 &copysize_r, sizeof(copysize_r), MPI_BYTE, rank_r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    cudaDeviceSynchronize();
    MPI_Sendrecv(buf_s, copysize_s*3, MPI_FLOAT, rank_s, 0,
                 &destbuf[offset][0], copysize_r*3, MPI_FLOAT, rank_r, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    cudaCheckError();


    return;

}

void gpuConstraintsD2D(rvec* x_d_ptr, const t_commrec *cr, cudaStream_t stream, matrix box,
                       bool setupRequired)
{


    gpuD2DConstraintsData* gpuD2D = &gpuD2DConsMod;

    if (setupRequired)
    {
        dd_setup_gpu_d2d_constraints(cr->dd, cr->dd->constraint_comm, box);

        //perform non-local data exchange directly on GPU.
        //forward
        cudaMemcpyAsync(gpuD2D->sendBufFWIndexMap_d, gpuD2D->sendBufFWIndexMap,
                        gpuD2D->sendBufFWIndexMapSize*sizeof(int), cudaMemcpyHostToDevice, stream);
        //backward
        cudaMemcpyAsync(gpuD2D->sendBufBWIndexMap_d, gpuD2D->sendBufBWIndexMap,
                        gpuD2D->sendBufBWIndexMapSize*sizeof(int), cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(gpuD2D->d2dConstraintsParams_d, gpuD2D->d2dConstraintsParams,
                        sizeof(gpuConstraintsD2Dparams_t), cudaMemcpyHostToDevice, stream);

    }


    const int threadsPerBlock = 128;

    dim3 blocks(((gpuD2D->sendBufFWIndexMapSize)+threadsPerBlock-1)/threadsPerBlock, 1, 1);
    dim3 threads(threadsPerBlock, 1, 1);

    //pack forward buffer
    pack_send_buf_kernel_constraints<<< blocks, threads, 0, stream>>> (gpuD2D->sendBufFW_d, (float*) x_d_ptr,
                                                                       gpuD2D->sendBufFWIndexMap_d,
                                                                       gpuD2D->d2dConstraintsParams_d,
                                                                       DIRFW);

    //pack backward buffer
    blocks.x = (gpuD2D->sendBufBWIndexMapSize+threadsPerBlock-1)/threadsPerBlock;
    pack_send_buf_kernel_constraints<<< blocks, threads, 0, stream>>> (gpuD2D->sendBufBW_d, (float*) x_d_ptr,
                                                                       gpuD2D->sendBufBWIndexMap_d,
                                                                       gpuD2D->d2dConstraintsParams_d, DIRBW);

    //send forward buffer
    gpuConstraintsXCommCoordSendRecv_D2D(x_d_ptr, gpuD2D->sendBufFW_d, DIRFW);

    //send backward buffer
    gpuConstraintsXCommCoordSendRecv_D2D(x_d_ptr, gpuD2D->sendBufBW_d, DIRBW);

    return;
}
