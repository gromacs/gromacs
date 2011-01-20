/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
*
*
* Gromacs                               Copyright (c) 1991-2005
* David van der Spoel, Erik Lindahl, University of Groningen.
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* To help us fund GROMACS development, we humbly ask that you cite
* the research papers on the package. Check out http://www.gromacs.org
* 
* And Hey:
* Gnomes, ROck Monsters And Chili Sauce
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef GMX_LIB_MPI 
#include <mpi.h>
#endif
#ifdef GMX_THREADS 
#include "tmpi.h"
#endif

#include "smalloc.h"
#include "gmx_parallel_3dfft.h"
#include "gmx_fft.h"
#include "gmxcomplex.h"
#include "gmx_fatal.h"

#include "fft5d.h"

struct gmx_parallel_3dfft  { 
    fft5d_plan p1,p2;    
};


static int *copy_int_array(int n,int *src)
{
    int *dest,i;

    dest = (int*)malloc(n*sizeof(int));
    for(i=0; i<n; i++)
        dest[i] = src[i];

    return dest;
}

static int *make_slab2grid(int nnodes,int ngrid)
{
    int *s2g,i;

    s2g = (int*)malloc((nnodes+1)*sizeof(int));
    for(i=0; i<nnodes+1; i++) {
        /* We always round up */
        s2g[i] = (i*ngrid + nnodes - 1)/nnodes;
    }

    return s2g;
}

int
gmx_parallel_3dfft_init   (gmx_parallel_3dfft_t *    pfft_setup,
                           ivec                      ndata,
						   real **                   real_data,
						   t_complex **              complex_data,
                           MPI_Comm                  comm[2],
                           int *                     slab2index_major,
                           int *                     slab2index_minor,
                           gmx_bool                      bReproducible)
{
    int rN=ndata[2],M=ndata[1],K=ndata[0];
    int flags = FFT5D_REALCOMPLEX | FFT5D_ORDER_YZ; /* FFT5D_DEBUG */
    MPI_Comm rcomm[]={comm[1],comm[0]};
    int Nb,Mb,Kb; /* dimension for backtransform (in starting order) */
    
    snew(*pfft_setup,1);
    if (bReproducible) flags |= FFT5D_NOMEASURE; 
    
    if (!(flags&FFT5D_ORDER_YZ)) { 
        Nb=M;Mb=K;Kb=rN;		
    } else {
        Nb=K;Mb=rN;Kb=M;  /* currently always true because ORDER_YZ always set */
    }
    
    (*pfft_setup)->p1 = fft5d_plan_3d(rN,M,K,rcomm, flags, (t_complex**)real_data, complex_data);
    
    (*pfft_setup)->p2 = fft5d_plan_3d(Nb,Mb,Kb,rcomm,
                                      (flags|FFT5D_BACKWARD|FFT5D_NOMALLOC)^FFT5D_ORDER_YZ, complex_data, (t_complex**)real_data);
    
    return (*pfft_setup)->p1 != 0 && (*pfft_setup)->p2 !=0;
}


static int
fft5d_limits(fft5d_plan p, 
             ivec                      local_ndata,
             ivec                      local_offset,
             ivec                      local_size) 
{
    int N1,M0,K0,K1,*coor;
    fft5d_local_size(p,&N1,&M0,&K0,&K1,&coor);  /* M0=MG/P[0], K1=KG/P[1], NG,MG,KG global sizes */
    
    local_offset[2]=0;
    local_offset[1]=p->oM[0];  /*=p->coor[0]*p->MG/p->P[0]; */
    local_offset[0]=p->oK[0];  /*=p->coor[1]*p->KG/p->P[1]; */
    
    local_ndata[2]=p->rC[0];
    local_ndata[1]=p->pM[0]; 
    local_ndata[0]=p->pK[0]; 
    
    if ((!(p->flags&FFT5D_BACKWARD)) && (p->flags&FFT5D_REALCOMPLEX)) {
        local_size[2]=p->C[0]*2;
    } else {
        local_size[2]=p->C[0];
    }
    local_size[1]=p->pM[0]; 
    local_size[0]=p->pK[0]; 
    return 0;
}

int
gmx_parallel_3dfft_real_limits(gmx_parallel_3dfft_t      pfft_setup,
                               ivec                      local_ndata,
                               ivec                      local_offset,
                               ivec                      local_size) {
    return fft5d_limits(pfft_setup->p1,local_ndata,local_offset,local_size);
}

static void reorder_ivec_yzx(ivec v)
{
    real tmp;

    tmp   = v[0];
    v[XX] = v[2];
    v[ZZ] = v[1];
    v[YY] = tmp;
}

int
gmx_parallel_3dfft_complex_limits(gmx_parallel_3dfft_t      pfft_setup,
                                  ivec                      complex_order,
                                  ivec                      local_ndata,
                                  ivec                      local_offset,
                                  ivec                      local_size) 
{
    int ret;

    /* For now everything is in-order, but prepare to save communication by avoiding transposes */
    complex_order[0] = 0;
    complex_order[1] = 1;
    complex_order[2] = 2;

    ret = fft5d_limits(pfft_setup->p2,local_ndata,local_offset,local_size);

    reorder_ivec_yzx(local_ndata);
    reorder_ivec_yzx(local_offset);
    reorder_ivec_yzx(local_size);

    return ret;
}


int
gmx_parallel_3dfft_execute(gmx_parallel_3dfft_t    pfft_setup,
                           enum gmx_fft_direction  dir,
                           void *                  in_data,
                           void *                  out_data) {
    if ((!(pfft_setup->p1->flags&FFT5D_REALCOMPLEX)) ^ (dir==GMX_FFT_FORWARD ||dir==GMX_FFT_BACKWARD)) { 
        gmx_fatal(FARGS,"Invalid transform. Plan and execution don't match regarding reel/complex");
    }
    if (dir==GMX_FFT_FORWARD || dir==GMX_FFT_REAL_TO_COMPLEX) {
        fft5d_execute(pfft_setup->p1,0);
    } else {
        fft5d_execute(pfft_setup->p2,0);
    }
    return 0;
}

int
gmx_parallel_3dfft_destroy(gmx_parallel_3dfft_t    pfft_setup) {
    fft5d_destroy(pfft_setup->p2);
    fft5d_destroy(pfft_setup->p1);
    sfree(pfft_setup);
    return 0;
}





