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
#include "gromacs/domdec/domdec_gpu.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/gpu_utils/cudautils.cuh"

/* float3 operator definitions */

static __device__ float3 operator+(const float3 &a, const float3 &b)
{
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

static __device__ float3 operator*(const bool &a, const float3 &b)
{
    return make_float3(a*b.x, a*b.y, a*b.z);
}

/* pack non-local coordinate data buffer on the GPU using pre-populated "map" containing index information */
__global__ void pack_send_buf_kernel(float               *data_packed,
                                     const float         *data,
                                     const int           *map,
                                     const int            mapSize,
                                     const bool           bPBC,
                                     const float         *xshift)
{
    int     idx       = blockIdx.x*blockDim.x+threadIdx.x;
    float3 *data_dest = (float3 *)&data_packed[idx*3];

    if (idx < mapSize)
    {
        *data_dest = *((float3*) &data[3*map[idx]]) + bPBC*(*((float3*)xshift));
    }

    return;
}

/* This function is based on move_x, but instead of actually sending and recieving data
 * a map structure is created containing information on indexes to be sent/recieved.  This
 * map will later be passed to the GPU and used to pack a buffer, which will be sent
 * directly between GPUs. */

void dd_init_move_x_gpu(gmx_domdec_t              *dd,
                        matrix                     box)
{

    ddGpu_t                      *ddGpu        = &(dd->ddGpu);
    int                           nzone        = 1;
    gmx_domdec_comm_t            *comm         = dd->comm;;
    int                           nat_tot      = comm->atomRanges.numHomeAtoms();
    gmx_domdec_comm_dim_t        *cd           = &comm->cd[0];
    bool                          bPBC         = (dd->ci[dd->dim[0]] == 0);
    const gmx_domdec_ind_t       &ind          = cd->ind[0];
    int                           size         = ind.nsend[nzone+1];
    const gmx::RangePartitioning &atomGrouping = dd->atomGrouping();

    if (dd->ndim > 1)
    {
        gmx_fatal(FARGS, "Error: dd->ndim > 1 is not yet supported in domdec_gpu");
    }

    if (!cd->receiveInPlace)
    {
        gmx_fatal(FARGS, "Error: out-of-place recieve is not yet supported in domdec_gpu");
    }

    if (bPBC && dd->bScrewPBC && dd->dim[0] == XX)
    {
        gmx_fatal(FARGS, "Error: screw is not yet supported in domdec_gpu\n");
    }

    if (size > ddGpu->indexMapSize) //(re-)allocate required structures
    {

        ddGpu->indexMapSize = size;

        if (ddGpu->indexMap)
        {
            free(ddGpu->indexMap);
            cudaFree(ddGpu->indexMap_d); //TODO use gmx API
        }
        ddGpu->indexMap = (int*) malloc(size*sizeof(int));
        cudaMalloc(&ddGpu->indexMap_d, size*sizeof(int)); //TODO use gmx API

        if (ddGpu->sendBuf_d)
        {
            cudaFree(ddGpu->sendBuf_d); //TODO use gmx API
        }

        cudaMalloc(&ddGpu->sendBuf_d, size*sizeof(rvec));         //TODO use gmx API

        cudaMalloc((void**) &(ddGpu->xshift_d), 3*sizeof(float)); //TODO use gmx API

    }

    ddGpu->indexMapSize = ind.nsend[nzone+1];  //size of map
    ddGpu->sendRank     = dd->neighbor[0][1];  //rank to send data to
    ddGpu->recvRank     = dd->neighbor[0][0];  //rank to recv data from
    ddGpu->bPBC         = bPBC;                // Period Boundary Conditions for this rank
    ddGpu->recvOffset   = nat_tot;             //offset for data recieved by this rank
    copy_rvec(box[dd->dim[0]], ddGpu->xshift); //three components of shift for PBC

    int n = 0;
    for (int g : ind.index)
    {
        for (int j : atomGrouping.block(g))
        {
            ddGpu->indexMap[n] = j; //add index to index map
            n++;
        }
    }

    // copy index map and shift to device
    //TODO use gmx API
    cudaMemcpy(ddGpu->indexMap_d, ddGpu->indexMap,
               ddGpu->indexMapSize*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(ddGpu->xshift_d, ddGpu->xshift, 3*sizeof(float), cudaMemcpyHostToDevice);

    // Setup send and recv sizes
    ddGpu->copysize_s = ddGpu->indexMapSize; //size sent by this rank
    // exchange sizes with remote ranks to get recv sizes
    MPI_Sendrecv(&ddGpu->indexMapSize, sizeof(ddGpu->indexMapSize), MPI_BYTE,
                 ddGpu->sendRank, 0, &ddGpu->copysize_r, sizeof(int), MPI_BYTE, ddGpu->recvRank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return;
}


void dd_move_x_gpu(gmx_domdec_t *dd, rvec* x_d_ptr, void *stream_ptr)
{
    ddGpu_t    * ddGpu  = &(dd->ddGpu);
    cudaStream_t stream = (cudaStream_t) stream_ptr;

    // launch kernel to pack send buffer

    const int          threadsPerBlock = 128;

    KernelLaunchConfig config;
    config.blockSize[0]     = threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (ddGpu->indexMapSize+threadsPerBlock-1)/threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;
    config.stream           = stream;

    float           *sendBuf  = ddGpu->sendBuf_d;
    const float     *x_d      = (float*) x_d_ptr;
    const int       *indexMap = ddGpu->indexMap_d;
    const int        size     = ddGpu->indexMapSize;
    const bool       bPBC     = ddGpu->bPBC;
    const float     *xshift   = ddGpu->xshift_d;

    auto             kernelFn                = pack_send_buf_kernel;

    const auto       kernelArgs   = prepareGpuKernelArguments(kernelFn, config, &sendBuf, &x_d, &indexMap,
                                                              &size, &bPBC, &xshift);

    launchGpuKernel(kernelFn, config, nullptr, "move_x_gpu", kernelArgs);
    cudaStreamSynchronize(stream);

    // perform halo exchange directly in device buffers

    MPI_Request request[2];
    MPI_Status  status[2];

    MPI_Irecv(&x_d_ptr[ddGpu->recvOffset][0], ddGpu->copysize_r*3, MPI_FLOAT, ddGpu->recvRank, 0,
              MPI_COMM_WORLD, &request[1]);

    MPI_Isend(ddGpu->sendBuf_d, ddGpu->copysize_s*3, MPI_FLOAT, ddGpu->sendRank, 0,
              MPI_COMM_WORLD, &request[0]);

    MPI_Waitall(2, request, status);

    return;
}
