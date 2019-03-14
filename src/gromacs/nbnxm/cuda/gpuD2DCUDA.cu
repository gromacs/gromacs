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
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/nbnxm/cuda/gpuUpdateConstraintsCUDA.h"
#include "gromacs/nbnxm/cuda/gpuD2DCUDA.h"

#include "gromacs/domdec/domdec_constraints.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/domdec_specatomcomm.h"



struct gpuD2DData{

    /* pointers for inter-GPU communication index maps and data buffers */
    int           * d2dBufIndexMap;
    int           * d2dBufIndexMap_d;
    int             d2dBufIndexMapSize;
    float         * sendBuf_d;
    gpuD2Dparams_t* d2dParams;
    gpuD2Dparams_t* d2dParams_d;
    float         * recvBuf_d;

    t_commrec     * commrecstatic;
    bool            copyparams;

    int             copysize_s;
    int             copysize_r;

};


//TEMPORARY Global variable to store a copy of above struct within this module.
//TODO move this to a central location (e.g. gpu_nbv) and pass through fn args.
static gpuD2DData gpuD2DMod;



/*-------------------------------- CUDA kernels --------------------------------*/

/* pack non-local coordinate data buffer on the GPU using pre-populated "map" containing index information */
__global__ void pack_send_buf_kernel(float               * data_packed,
                                     const float         * data,
                                     const int           * map,
                                     const gpuD2Dparams_t* d2dParams)
{

    int   idx       = blockIdx.x*blockDim.x+threadIdx.x;
    int   data_size = d2dParams->mapSize;
    int   bPBC      = d2dParams->bPBC;

    float shiftX = d2dParams->xshift[0];
    float shiftY = d2dParams->xshift[1];
    float shiftZ = d2dParams->xshift[2];



    if (idx < data_size)
    {

        if (bPBC)  //we have periodic boundary conditions so shift is needed
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


/* unpack non-local force data buffer on the GPU using pre-populated "map" containing index information */
__global__ void unpack_recv_buf_kernel(float               * data,
                                       const float         * data_packed,
                                       const int           * map,
                                       const gpuD2Dparams_t* d2dParams)
{

    int   idx       = blockIdx.x*blockDim.x+threadIdx.x;
    int   data_size = d2dParams->mapSize;
    int   bPBC      = d2dParams->bPBC;

    float shiftX = d2dParams->fshift[0];
    float shiftY = d2dParams->fshift[1];;
    float shiftZ = d2dParams->fshift[2];;

    if (idx < data_size)
    {

        if (bPBC)  //we have periodic boundary conditions so shift is needed
        {
            data[3*map[idx]]   += data_packed[idx*3]+shiftX;
            data[3*map[idx]+1] += data_packed[idx*3+1]+shiftY;
            data[3*map[idx]+2] += data_packed[idx*3+2]+shiftZ;
        }
        else
        {
            data[3*map[idx]]   += data_packed[idx*3];
            data[3*map[idx]+1] += data_packed[idx*3+1];
            data[3*map[idx]+2] += data_packed[idx*3+2];
        }

    }
    return;
}

/*-------------------------------- End CUDA kernels-----------------------------*/

/* this function performs the set-up for a buffer for information of on the indices for t\
   hose atoms involved in inter-GPU atom coordinate communications. It (re-)allocates the bu\
   ffers (if necessary) and passes pointers back to the calling routine, in which the buffer\
   s are populated. This index information is then used in gpuBufferOpsCopyRvecToNbatReal wh\
   ich packs a device-side buffer and performs the inter-GPU data communication using CUDA m\
   emory copies.
 */

static void gpuBufferOpsInitXSendBufIndexMap(int** mapptr_in,
                                             gpuD2Dparams_t** d2dParamptr_in, int size)
{

    gpuD2DData* gpuD2D = &gpuD2DMod;

    if (size > gpuD2D->d2dBufIndexMapSize)
    {

        gpuD2D->d2dBufIndexMapSize = size;

        if (gpuD2D->d2dBufIndexMap)
        {
            free(gpuD2D->d2dBufIndexMap);
            cudaFree(gpuD2D->d2dBufIndexMap_d);
        }
        gpuD2D->d2dBufIndexMap = (int*) malloc(size*sizeof(int));
        cudaMalloc(&gpuD2D->d2dBufIndexMap_d, size*sizeof(int));

        if (gpuD2D->sendBuf_d)
        {
            cudaFree(gpuD2D->sendBuf_d);
        }
        if (gpuD2D->recvBuf_d)
        {
            cudaFree(gpuD2D->recvBuf_d);
        }


        cudaMalloc(&gpuD2D->sendBuf_d, size*sizeof(rvec));
        cudaMalloc(&gpuD2D->recvBuf_d, size*sizeof(rvec));

        if (gpuD2D->d2dParams)
        {
            free(gpuD2D->d2dParams);
            cudaFree(gpuD2D->d2dParams_d);
        }

        gpuD2D->d2dParams = (gpuD2Dparams_t*) malloc(sizeof(gpuD2Dparams_t));
        cudaMalloc(&gpuD2D->d2dParams_d, sizeof(gpuD2Dparams_t));

    }

    *mapptr_in      = gpuD2D->d2dBufIndexMap;
    *d2dParamptr_in = gpuD2D->d2dParams;

    return;
}

static void gpuBufferOpsGetD2DParamsPtr(gpuD2Dparams_t** d2dParamptr_in)
{

    gpuD2DData* gpuD2D = &gpuD2DMod;

    *d2dParamptr_in = gpuD2D->d2dParams;

    return;
}



// This function is based on move_x/move_f, but instead of actually
// sending and recieving data a map structure is created containing
// information on indexes to be sent/recieved.  This map will later
// be passed to the GPU and used to pack a buffer, which will be sent
// directly between GPUs.

void dd_setup_gpu_d2d(gmx_domdec_t              *dd,
                      matrix                     box)
{


    int                    nzone, nat_tot;
    gmx_domdec_comm_t     *comm;
    gmx_domdec_comm_dim_t *cd;
    rvec                   shift = {0, 0, 0};
    gmx_bool               bPBC, bScrew;

    comm = dd->comm;

    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();


    if (dd->ndim > 1)
    {
        gmx_fatal(FARGS, "Error: dd->ndim > 1 is not yet supported in GPU D2D");
    }


    nzone   = 1;
    nat_tot = comm->atomRanges.numHomeAtoms();
    for (int d = 0; d < dd->ndim; d++)
    {
        bPBC   = (dd->ci[dd->dim[d]] == 0);
        bScrew = (bPBC && dd->bScrewPBC && dd->dim[d] == XX);
        if (bPBC)
        {
            copy_rvec(box[dd->dim[d]], shift);
        }
        cd = &comm->cd[d];

        if (!cd->receiveInPlace)
        {
            gmx_fatal(FARGS, "Error: out-of-place recieve is not yet supported in GPU D2D");
        }


        for (const gmx_domdec_ind_t &ind : cd->ind)
        {

            // Create and populate a data structure containing information on this communication pattern

            int           * sendBufIndexMap;
            gpuD2Dparams_t* d2dParams;

            gpuBufferOpsInitXSendBufIndexMap(&sendBufIndexMap, &d2dParams, ind.nsend[nzone+1]);


            //size of map
            d2dParams->mapSize = ind.nsend[nzone+1];
            //rank to send data to
            d2dParams->sendRank = dd->neighbor[d][1];
            //rank to recv data from
            d2dParams->recvRank = dd->neighbor[d][0];
            // Period Boundary Conditions for this rank
            d2dParams->bPBC = bPBC;
            //offset for data recieved by this rank
            d2dParams->recvOffset = nat_tot;

            //three components of shift for PBC
            d2dParams->xshift[0] = shift[0];
            d2dParams->xshift[1] = shift[1];
            d2dParams->xshift[2] = shift[2];


            int                        n          = 0;
            if (!bPBC)
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        sendBufIndexMap[n] = j;
                        n++;
                    }
                }
            }
            else if (!bScrew)
            {
                for (int g : ind.index)
                {
                    for (int j : atomGrouping.block(g))
                    {
                        /* We need to shift the coordinates */
                        for (int d = 0; d < DIM; d++)
                        {
                            sendBufIndexMap[n] = j;
                        }
                        n++;
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "Error: bScrew=true is not yet supported in DD GPU\n");

            }


        }

    }
}


// add force shift components to D2D communication data structure
void dd_gpu_d2d_add_fshift_to_params(gmx_domdec_t *dd, rvec *fshift)
{

    ivec                   vis;
    int                    is;

    /* Determine which shift vector we need */
    clear_ivec(vis);
    vis[dd->dim[0]] = 1;
    is              = IVEC2IS(vis);

    gpuD2Dparams_t* d2dParams;
    gpuBufferOpsGetD2DParamsPtr(&d2dParams);

    // //three components of shift for PBC
    d2dParams->fshift[0] = fshift[is][0];
    d2dParams->fshift[1] = fshift[is][1];
    d2dParams->fshift[2] = fshift[is][2];

}



/* function to perform non-local X coordinate communication on GPU */
static void gpuBufferOpsXCommCoordSendRecv_D2D(rvec* x_d_ptr,
                                               float* buf_s,
                                               int copysize_s,
                                               int rank_s,
                                               int offset, int rank_r)
{

    gpuD2DData* gpuD2D = &gpuD2DMod;

    if (gpuD2D->copyparams)
    {
        int copysize_r;
        //get remote copysize
        MPI_Sendrecv(&copysize_s, sizeof(copysize_s), MPI_BYTE, rank_s, 0,
                     &copysize_r, sizeof(copysize_r), MPI_BYTE, rank_r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        gpuD2D->copysize_s = copysize_s;
        gpuD2D->copysize_r = copysize_r;
    }


    MPI_Request request[2];
    MPI_Status  status[2];

    MPI_Irecv(&x_d_ptr[offset][0], gpuD2D->copysize_r*3, MPI_FLOAT,
              rank_r, 0,
              MPI_COMM_WORLD, &request[1]);

    MPI_Isend(buf_s, gpuD2D->copysize_s*3, MPI_FLOAT, rank_s, 0,
              MPI_COMM_WORLD, &request[0]);


    MPI_Waitall(2, request, status);

    cudaCheckError();


    return;

}


/* function to perform non-local F communication on GPU */
static void gpuBufferOpsFCommCoordSendRecv_D2D(rvec* f_d_ptr, float* buf_s,
                                               int copysize_s, int rank_s,
                                               int offset, int rank_r)
{


    gpuD2DData* gpuD2D = &gpuD2DMod;


    MPI_Request request[2];
    MPI_Status  status[2];

    MPI_Irecv(buf_s, gpuD2D->copysize_s*3, MPI_FLOAT, rank_s, 0,
              MPI_COMM_WORLD, &request[1]);

    MPI_Isend(&(f_d_ptr[offset][0]), gpuD2D->copysize_r*3, MPI_FLOAT, rank_r, 0,
              MPI_COMM_WORLD, &request[0]);


    MPI_Waitall(2, request, status);

    cudaCheckError();


    return;

}

void gpuBufferOpsXCommCoord(rvec* x_d_ptr, cudaStream_t stream)
{


    gpuD2DData* gpuD2D = &gpuD2DMod;

    if (gpuD2D->copyparams)
    {
        cudaMemcpyAsync(gpuD2D->d2dBufIndexMap_d, gpuD2D->d2dBufIndexMap,
                        gpuD2D->d2dBufIndexMapSize*sizeof(int), cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(gpuD2D->d2dParams_d, gpuD2D->d2dParams, sizeof(gpuD2Dparams_t),
                        cudaMemcpyHostToDevice, stream);
    }

    const int threadsPerBlock = 128;

    dim3 blocks((gpuD2D->d2dBufIndexMapSize+threadsPerBlock-1)/threadsPerBlock, 1, 1);
    dim3 threads(threadsPerBlock, 1, 1);

    pack_send_buf_kernel<<< blocks, threads, 0, stream>>> (gpuD2D->sendBuf_d, (float*) x_d_ptr,
                                                           gpuD2D->d2dBufIndexMap_d, gpuD2D->d2dParams_d);

    cudaStreamSynchronize(stream);

    int packeddatasize = gpuD2D->d2dParams->mapSize;
    int rank_s         = gpuD2D->d2dParams->sendRank;
    int rank_r         = gpuD2D->d2dParams->recvRank;
    int offset         = gpuD2D->d2dParams->recvOffset;

    //Push non-local data to rank_s
    gpuBufferOpsXCommCoordSendRecv_D2D(x_d_ptr, gpuD2D->sendBuf_d, packeddatasize, rank_s,
                                       offset, rank_r);

    //MPI_Barrier(gpuD2D->commrecstatic->dd->mpi_comm_all);

    if (gpuD2D->copyparams)
    {
        gpuD2D->copyparams = false;
    }
    return;
}


void gpuBufferOpsFCommCoord(rvec* f_d_ptr, cudaStream_t stream)
{

    //perform D2D comms


    gpuD2DData* gpuD2D = &gpuD2DMod;

    const int   threadsPerBlock = 128;



    int packeddatasize = gpuD2D->d2dParams->mapSize;
    int rank_s         = gpuD2D->d2dParams->sendRank;
    int rank_r         = gpuD2D->d2dParams->recvRank;
    int offset         = gpuD2D->d2dParams->recvOffset;

    cudaCheckError();

    cudaStreamSynchronize(stream);

    //MPI_Barrier(gpuD2D->commrecstatic->dd->mpi_comm_all);


    cudaCheckError();

    // this pulls data from rank_s
    gpuBufferOpsFCommCoordSendRecv_D2D(f_d_ptr, gpuD2D->recvBuf_d, packeddatasize, rank_s,
                                       offset, rank_r);




    cudaCheckError();
    cudaMemcpyAsync(gpuD2D->d2dParams_d, gpuD2D->d2dParams, sizeof(gpuD2Dparams_t),
                    cudaMemcpyHostToDevice, stream);
    dim3 blocks((gpuD2D->d2dBufIndexMapSize+threadsPerBlock-1)/threadsPerBlock, 1, 1);
    dim3 threads(threadsPerBlock, 1, 1);

    cudaCheckError();
    unpack_recv_buf_kernel<<< blocks, threads, 0, stream>>> ((float*) f_d_ptr,
                                                             gpuD2D->recvBuf_d,
                                                             gpuD2D->d2dBufIndexMap_d,
                                                             gpuD2D->d2dParams_d);

    cudaCheckError();

    return;

}


void gpuD2DSetCommrec(const t_commrec *cr)
{

    gpuD2DData* gpuD2D = &gpuD2DMod;

    gpuD2D->commrecstatic = (t_commrec*) cr;

    gpuD2D->copyparams = true;
}
