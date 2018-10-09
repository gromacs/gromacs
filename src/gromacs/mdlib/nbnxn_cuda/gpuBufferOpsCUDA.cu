/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
/*! \file
 * \brief functionality for performing buffer operations on the GPU and  
 * inter-GPU data movement facilities. 
 * \author Alan Gray <alang@nvidia.com> and Jon Vincent <jvincent@nvidia.com> 
 */

#include "gmxpre.h"

#include <stdio.h>
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/mdlib/nb_verlet.h"

#include "gpuBufferOpsCUDA.h"

#if defined(_MSVC)
#include <limits>
#endif


/*-------------------------------- CUDA kernels-------------------------------- */
/*------------------------------------------------------------------------------*/

/* CUDA kernel for GPU version of copy_rvec_to_nbat_real (nbnxn_atomdata.cpp) */
__global__ void copy_rvec_to_nbat_real_cuda_kernel(int ncxy, float* xnb,int g,
                                                   gmx_bool FillLocal,const rvec* x,
                                                   const int* a,int* na_all,
                                                   int* cxy_ind, int cell0,int na_sc,
                                                   int stride){

  
    const float farAway = -1000000;

    /* map cell-level parallelism to y component of CUDA block index */
    int cxy =blockIdx.y;

    if ( cxy < ncxy){
    
        int na=na_all[cxy];
        int a0 = (cell0 + cxy_ind[cxy])*na_sc;
        int na_round;
        if (g == 0 && FillLocal)
        {
            na_round =
                (cxy_ind[cxy+1] - cxy_ind[cxy])*na_sc;
        }
        else
        {
            /* We fill only the real particle locations.
             * We assume the filling entries at the end have been
             * properly set before during pair-list generation.
             */
            na_round = na;
        }
    
        /* map parlallism within a cell to x component of CUDA block index linearized 
         * with threads within a block */
        int i, j0;
        i =blockIdx.x*blockDim.x+threadIdx.x;
    
        j0 =a0*stride;

        /* perform conversion of each element */
        if (i < na_round){
            if (i < na)
            {
    
                xnb[j0+4*i] = x[a[a0+i]][XX];
                xnb[j0+4*i+1] = x[a[a0+i]][YY];
                xnb[j0+4*i+2] = x[a[a0+i]][ZZ];
            }
            else{
                xnb[j0+4*i] = farAway;
                xnb[j0+4*i+1] = farAway;
                xnb[j0+4*i+2] = farAway;      
            }
        }
    }  
    
}


/*-------------------------------- End CUDA kernels-----------------------------*/



/*------------------------------------------------------------------------------*/



/* static variables internal to module. Separate instances are needed
 * for each MPI rank (because, when using threadMPI, we need to avoid 
 * interference and allow a global view of memory)
 * and in some cases for local/nonlocal parts which are in different
 * streams */

static int* a_d[BO_MAX_RANKS]; //input buffer index mapping for X-buffer ops on device
static int* na_all_d[BO_MAX_RANKS][2]; //number of atoms for X-buffer ops on device
static int* cxy_ind_d[BO_MAX_RANKS][2]; //cell index mapping for X-buffer ops on device
static rvec* x_d[BO_MAX_RANKS]; //X on device

/* allocation trackers */
static volatile int a1maxf[BO_MAX_RANKS]; //max size allocated for size-a1 data arrays
static int ncxy_max[BO_MAX_RANKS][2]; //max size allocated for size-ncxy data arrays
static int xsize_max[BO_MAX_RANKS];//max size allocated for size-xsize data arrays
static int a_nalloc_max[BO_MAX_RANKS];//max size allocated for size-a_nalloc data arrays


/* flags for initialization state */
static volatile bool x_already_alloced[BO_MAX_RANKS];
static volatile bool f_already_alloced[BO_MAX_RANKS];
static bool init_required_x[BO_MAX_RANKS][2];

/* local copy of parameters */
static int nAtomsAllThisStep[BO_MAX_RANKS];
static int nAtomsLocalThisStep[BO_MAX_RANKS];

/* boolean variable keeping track of whether or not GPU buffers active in a specific timestep. */
static bool bGpuBufferOpsActiveThisTimestep[BO_MAX_RANKS];



/* function to flag from external routine that buffer ops initialization is required */
static void gpuBufferOpsFlagInitRequired(){
    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    init_required_x[rank][0]=true;
    init_required_x[rank][1]=true;

    return;
}


/* function to allocate X and F */
static void gpuBufferOpsAllocXF(int a1)
{  

    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //allocate 2X the required space, to allow for expansion with re-distribution
    int mallocfactor=2;
  

        if (a1 > xsize_max[rank]){

            xsize_max[rank]=mallocfactor*a1;
            if(x_d[rank]) cudaFree(x_d[rank]);
            cudaMalloc(&x_d[rank],a1*mallocfactor*sizeof(rvec));	 

            x_already_alloced[rank]=true;
        }

  
    return;
}



/* GPU version of copy_rvec_to_nbat_real. See nbnxn_atomdata.cpp for CPU version */
void gpuBufferOpsCopyRvecToNbatReal(int ncxy,int g, gmx_bool FillLocal,
                                    gmx_nbnxn_gpu_t *gpu_nbv, const int* a,
                                    int a_nalloc,const int* na_all,
                                    const int* cxy_ind,
                                    int cell0,
                                    int na_sc, int iloc, int stride, rvec* x
                                    ){


    int rank=0;
    int size=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if (size > BO_MAX_RANKS){

        printf("Error, nbnxn_cuda_buffer_ops has been compiled with BO_MAX_RANKS=%d but MPI_Comm_size is %d\n",BO_MAX_RANKS,size);

        MPI_Abort(MPI_COMM_WORLD,1);

    }
    
    cu_atomdata_t       *adat    = gpu_nbv->atdat;
    cudaStream_t         stream  = gpu_nbv->stream[iloc];

  
    int na_round_max=0;
    for (int cxy = 0; cxy < ncxy; cxy++)
    {
      
      
        int na=na_all[cxy];
        int na_round;
        if (g == 0 && FillLocal)
        {
            na_round =
                (cxy_ind[cxy+1] - cxy_ind[cxy])*na_sc;
        }
        else
        {
            /* We fill only the real particle locations.
             * We assume the filling entries at the end have been
             * properly set before during pair-list generation.
             */
            na_round = na;
        }
      
      
        if (na_round > na_round_max) na_round_max=na_round;
      
    }


     
    // Perform allocations and memcopies if necessary
    if(init_required_x[rank][iloc])
    {

        init_required_x[rank][iloc]=false;

        if(iloc==0){	  
            if(a_nalloc > a_nalloc_max[rank]){
                a_nalloc_max[rank]=a_nalloc;
                if(a_d[rank])
                    cudaFree(a_d[rank]);
                cudaMalloc(&a_d[rank],a_nalloc*sizeof(int));
            }
        }
      
        if(ncxy > ncxy_max[rank][iloc]){
            ncxy_max[rank][iloc]=ncxy;
            if(na_all_d[rank][iloc])
                cudaFree(na_all_d[rank][iloc]);
            cudaMalloc(&na_all_d[rank][iloc],ncxy*sizeof(int));
	
            if(cxy_ind_d[rank][iloc])
                cudaFree(cxy_ind_d[rank][iloc]);
            cudaMalloc(&cxy_ind_d[rank][iloc],ncxy*sizeof(int));
	
        }
      
        cudaCheckError();
      
        if(iloc==0) cudaMemcpy(a_d[rank],a,a_nalloc*sizeof(int),cudaMemcpyHostToDevice);
        cudaMemcpy(na_all_d[rank][iloc],na_all,ncxy*sizeof(int),cudaMemcpyHostToDevice);
        cudaMemcpy(cxy_ind_d[rank][iloc],cxy_ind,ncxy*sizeof(int),cudaMemcpyHostToDevice);
	
	
        cudaCheckError();

      
    }

    
    
    int nCopyAtoms;
    int copyAtomStart;
    if (iloc==0){//local stream
        nCopyAtoms=nAtomsLocalThisStep[rank];
        copyAtomStart=0;
    }
    else{//non-local stream
        nCopyAtoms=nAtomsAllThisStep[rank]-nAtomsLocalThisStep[rank];
        copyAtomStart=nAtomsLocalThisStep[rank]; 
    }
    
    // TEMPORARY DATA COPY 
    // THIS WILL BE REMOVED FOR FULL NVIDIA OPTIMIZED VERSION
    // copy X-coordinate data to device
    cudaMemcpyAsync(&x_d[rank][copyAtomStart][0],&x[copyAtomStart][0],
                    nCopyAtoms*sizeof(rvec),
                    cudaMemcpyHostToDevice,stream);
    cudaCheckError();
    //END TEMPORARY DATA COPY
    
    // INTEGRATION FUNCTIONALITY COMMENTED OUT FOR THIS PARTIAL VERSION
    // if(iloc!=0){
        
    //     /* make sure x buffer has been copied by local stream */
    //     waitEvent(&localXDataReady[rank],&blocalXDataReady[rank],&stream);
        
    //     //perform non-local data exchange directly on GPU.
    //     gpuBufferOpsXCommCoord(x,stream);
        
    // }
    // END INTEGRATION FUNCTIONALITY 
    
  
    /* launch kernel on GPU */
    
    dim3 blocks(((na_round_max+1)+TPB-1)/TPB,ncxy,1);
    dim3 threads(TPB,1,1);    
    
    copy_rvec_to_nbat_real_cuda_kernel
        <<<blocks,threads,0,stream>>>
        (ncxy,(float*) adat->xq,g,FillLocal,
         x_d[rank],a_d[rank],na_all_d[rank][iloc],
         cxy_ind_d[rank][iloc],cell0,na_sc, stride);

    cudaCheckError();
        
    return;
}

/* This function is called from PP ranks and performs initialization
 of the new module by allocating internal buffers (if necessary) ,
 populating internal variables, and returning a boolean value
 indicating if the GPU buffer operations are active within the
 timestep. */    
gmx_bool gpuBufferOpsTimestepInitFromPP(gmx_bool bNS,
                                        gmx_bool bUseGPU,
                                        int natoms_all,
                                        int natoms_local,
                                        int xsize)
{
					  
    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    nAtomsAllThisStep[rank]=natoms_all;
    nAtomsLocalThisStep[rank]=natoms_local;
  
    if(bNS) gpuBufferOpsFlagInitRequired();
  
    if(bUseGPU){
        //allocate scratch buffers
        gpuBufferOpsAllocXF(xsize);
    }
  
    bGpuBufferOpsActiveThisTimestep[rank]=true;
  
    if(bNS || !bUseGPU) bGpuBufferOpsActiveThisTimestep[rank]=false;
    
    return bGpuBufferOpsActiveThisTimestep[rank];

}


/* function to return a bool indicating if the GPU buffer operations
   are active within the timestep */
gmx_bool gpuBufferOpsActiveThisTimestep(){

    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
    return bGpuBufferOpsActiveThisTimestep[rank];
};

