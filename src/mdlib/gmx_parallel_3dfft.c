/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
*
* $Id$
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

#ifdef GMX_MPI
 
#include <mpi.h>


#include "gmx_parallel_3dfft.h"
#include "gmx_fft.h"
#include "gmxcomplex.h"
#include "gmx_fatal.h"

typedef struct {
    int             *sdisps;
    int             *scounts;
    int             *rdisps;
    int             *rcounts;
} alltoallv_t;

struct gmx_parallel_3dfft 
{
    int             nx;
    int             ny;
    int             nz;
    int             nzc;
    int             local_slab;
    int             nnodes;
    gmx_fft_t       fft_yz;
    gmx_fft_t       fft_x;
	void *          work_rawptr;
	void *          work2_rawptr;
    t_complex *     work;
    t_complex *     work2;
    int             *node2slab;
    int             *slab2grid_x;
    int             *slab2grid_y;
    alltoallv_t     *aav;
    MPI_Comm        comm;
};

static int *copy_int_array(int n,int *src)
{
    int *dest,i;

    dest = malloc(n*sizeof(int));
    for(i=0; i<n; i++)
        dest[i] = src[i];

    return dest;
}

static int *make_slab2grid(int nnodes,int ngrid)
{
    int *s2g,i;

    s2g = malloc((nnodes+1)*sizeof(int));
    for(i=0; i<nnodes+1; i++) {
        /* We always round up */
        s2g[i] = (i*ngrid + nnodes - 1)/nnodes;
    }

    return s2g;
}

int
gmx_parallel_3dfft_init   (gmx_parallel_3dfft_t *    pfft_setup,
                           int                       ngridx,
                           int                       ngridy,
                           int                       ngridz,
                           int                       *node2slab,
                           int                       *slab2grid_x,
                           MPI_Comm                  comm,
                           bool                      bReproducible)
{
    gmx_parallel_3dfft_t p;
    int  nxy_n;
    void *p0;
    int  flags;
    
    flags = bReproducible ? GMX_FFT_FLAG_CONSERVATIVE : 0;
    
    p = malloc(sizeof(struct gmx_parallel_3dfft));
    
    if(p==NULL)
        return ENOMEM;
    
    p->nx  = ngridx;
    p->ny  = ngridy;
    p->nz  = ngridz;
    p->nzc = ngridz/2 + 1;

    MPI_Comm_rank( comm , &(p->local_slab) );
    
    MPI_Comm_dup( comm , &(p->comm) );

    MPI_Comm_size( p->comm , &p->nnodes);

    if (node2slab)
        p->node2slab = copy_int_array(p->nnodes,node2slab);
    else
        p->node2slab = NULL;

    if (p->node2slab)
        p->local_slab = p->node2slab[p->local_slab];
        
    MPI_Comm_dup( comm , &(p->comm) );

    MPI_Comm_size( p->comm , &p->nnodes);

    if (slab2grid_x)
        p->slab2grid_x = copy_int_array(p->nnodes+1,slab2grid_x);
    else
        p->slab2grid_x = make_slab2grid(p->nnodes,p->nx);
    p->slab2grid_y     = make_slab2grid(p->nnodes,p->ny);

    if (node2slab || p->nx % p->nnodes || p->ny % p->nnodes) {
        p->aav = malloc(sizeof(alltoallv_t));
        p->aav->sdisps  = malloc(p->nnodes*sizeof(int));
        p->aav->scounts = malloc(p->nnodes*sizeof(int));
        p->aav->rdisps  = malloc(p->nnodes*sizeof(int));
        p->aav->rcounts = malloc(p->nnodes*sizeof(int));
    } else {
        p->aav  = NULL;
    }

    /* initialize transforms */
    if ( ( gmx_fft_init_1d(&(p->fft_x),ngridx,flags) != 0 ) ||
         ( gmx_fft_init_2d_real(&(p->fft_yz),ngridy,ngridz,flags) != 0))
    {
        free(p);
        return -1;
    }

    /* Round up */
    nxy_n = p->nnodes*((p->nx + p->nnodes - 1)/p->nnodes)*
                      ((p->ny + p->nnodes - 1)/p->nnodes);

    
    p0               = malloc(sizeof(real)*2*(p->nzc)*nxy_n + 32);
    p->work_rawptr   = p0;
    p->work          = (void *) (((size_t) p0 + 32) & (~((size_t) 31)));

    p0               = malloc(sizeof(real)*2*(p->nzc)*nxy_n);
    p->work2_rawptr  = p0;
    p->work2         = (void *) (((size_t) p0 + 32) & (~((size_t) 31)));
    
    if(p->work == NULL || p->work2 == NULL)
    {
        if(p->work_rawptr != NULL)
            free(p->work_rawptr);
        if(p->work2_rawptr != NULL)
            free(p->work2_rawptr);
        free(p);
        return ENOMEM;
    }

    *pfft_setup = p;
    
    return 0;
}




int
gmx_parallel_3dfft_limits(gmx_parallel_3dfft_t      pfft_setup,
                          int *                     local_x_start,
                          int *                     local_nx,
                          int *                     local_y_start,
                          int *                     local_ny)
{
    int slab;

    slab = pfft_setup->local_slab;

    *local_x_start = pfft_setup->slab2grid_x[slab];
    *local_y_start = pfft_setup->slab2grid_y[slab];

    *local_nx = pfft_setup->slab2grid_x[slab+1] - (*local_x_start);
    *local_ny = pfft_setup->slab2grid_y[slab+1] - (*local_y_start);
    
    return 0;
}


                   
int
gmx_parallel_transpose_xy(t_complex *   data,
                          t_complex *   work,
                          int           nx,
                          int           ny,
                          int           local_slab,
                          int           *s2x,
                          int           *s2y,
                          int           nzc,
                          int           nnodes,
                          int           *node2slab,
                          alltoallv_t   *aav,
                          MPI_Comm      comm)
{
    int     i,j;
    int     local_nx,local_ny,blocksize,slab;
    
    local_nx = s2x[local_slab+1] - s2x[local_slab];
    local_ny = s2y[local_slab+1] - s2y[local_slab];
    
    /* A: Do a local transpose to get data continuous for communication.
    *     We can use NULL for the workarray since we do it out-of-place.
    */
    gmx_fft_transpose_2d_nelem(data,work,local_nx,ny,nzc,NULL);
    
    /* B: Parallel communication, exchange data blocks. */
    if (aav == NULL) {
        blocksize = local_nx*local_ny*nzc*2;
        MPI_Alltoall(work,
                     blocksize,
                     GMX_MPI_REAL,
                     data,
                     blocksize,
                     GMX_MPI_REAL,
                     comm);
    } else {
        for(i=0; i<nnodes; i++) {
            slab = (node2slab ? node2slab[i] : i);
            aav->sdisps [i] =                s2y[slab] *local_nx*nzc*2;
            aav->scounts[i] = (s2y[slab+1] - s2y[slab])*local_nx*nzc*2;
            aav->rdisps [i] = local_ny*               s2x[slab] *nzc*2;
            aav->rcounts[i] = local_ny*(s2x[slab+1] - s2x[slab])*nzc*2;
        }
        MPI_Alltoallv(work,
                      aav->scounts,aav->sdisps,
                      GMX_MPI_REAL,
                      data,
                      aav->rcounts,aav->rdisps,
                      GMX_MPI_REAL,
                      comm);
    }
        
    /* C: Copy entire blocks into place, so we have YXZ. */
    for(j=0;j<local_ny;j++)
    {
        for(i=0;i<nnodes;i++)
        {
            memcpy(work + j*nx*nzc + s2x[i]*nzc,
                   data + s2x[i]*local_ny*nzc + j*(s2x[i+1] - s2x[i])*nzc,
                   (s2x[i+1] - s2x[i])*nzc*sizeof(t_complex));
        }
    }
    return 0;
}

                       
int
gmx_parallel_3dfft(gmx_parallel_3dfft_t    pfft_setup,
                   enum gmx_fft_direction  dir,
                   void *                  in_data,
                   void *                  out_data)
{
    int          i,j,k;
    int          nx,ny,nz,nzc,nzr;
    int          local_x_start,local_nx;
    int          local_y_start,local_ny;    
    t_complex *  work;
    real *       rdata;
    t_complex *  cdata;
    t_complex *  ctmp;
    
    work    = pfft_setup->work;
    
    /* When we do in-place FFTs the data need to be embedded in the z-dimension,
     * so there is room for the complex data. This means the direct space
     * _grid_ (not data) dimensions will be nx*ny*(nzc*2), where nzc=nz/2+1.
     * If we do out-of-place transforms the direct space dimensions are simply
     * nx*ny*nz, and no embedding is used.
     * The complex dimensions are always ny*nx*nzc (note the transpose).
     *
     * The direct space _grid_ dimension is nzr.
     */
    
    nx  = pfft_setup->nx;
    ny  = pfft_setup->ny;
    nz  = pfft_setup->nz;
    nzc = pfft_setup->nzc;
    
    if(in_data == out_data)
    {
        nzr = 2*nzc;
    }
    else
    {
        nzr = nz;
    }

    gmx_parallel_3dfft_limits(pfft_setup,
                              &local_x_start,
                              &local_nx,
                              &local_y_start,
                              &local_ny);

    if(dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        rdata =      (real *)in_data  + local_x_start*ny*nzr;
        cdata = (t_complex *)out_data + local_x_start*ny*nzc;
        
        /* Perform nx local 2D real-to-complex FFTs in the yz slices.
         * When the input data is "embedded" for 3D-in-place transforms, this
         * must also be done in-place to get the data embedding right.
         * 
         * Note that rdata==cdata when we work in-place. 
         */
        for(i=0;i<local_nx;i++)
        {
            gmx_fft_2d_real(pfft_setup->fft_yz,
                            GMX_FFT_REAL_TO_COMPLEX,
                            rdata + i*ny*nzr,
                            cdata + i*ny*nzc);
        }
        
        /* Transpose to temporary work array */
        gmx_parallel_transpose_xy(cdata,
                                  work,
                                  nx,
                                  ny,
                                  pfft_setup->local_slab,
                                  pfft_setup->slab2grid_x,
                                  pfft_setup->slab2grid_y,
                                  nzc,
                                  pfft_setup->nnodes,
                                  pfft_setup->node2slab,
                                  pfft_setup->aav,
                                  pfft_setup->comm);

        /* Transpose from temporary work array in order YXZ to
         * the output array in order YZX. 
         */ 
        /* output cdata changes when nx or ny not divisible by nnodes */
        cdata = (t_complex *)out_data + local_y_start*nx*nzc;
        for(j=0;j<local_ny;j++)
        {
            gmx_fft_transpose_2d(work  + j*nzc*nx,
                                 cdata + j*nzc*nx,
                                 nx,
                                 nzc);
        }

        /* Perform local_ny*nzc complex FFTs along the x dimension */
        for(i=0;i<local_ny*nzc;i++)
        {
            gmx_fft_1d(pfft_setup->fft_x,
                       GMX_FFT_FORWARD,
                       cdata + i*nx,
                       work  + i*nx);
        }    
    
        /* Transpose back from YZX to YXZ. */
        for(j=0;j<local_ny;j++)
        {
            gmx_fft_transpose_2d(work  + j*nzc*nx,
                                 cdata + j*nzc*nx,
                                 nzc,
                                 nx);
        }
    }
    else if(dir == GMX_FFT_COMPLEX_TO_REAL)
    {
        cdata = (t_complex *)in_data  + local_y_start*nx*nzc;
        rdata =      (real *)out_data + local_x_start*ny*nzr;
        
        /* If we are working in-place it doesn't matter that we destroy
         * input data. Otherwise we use an extra temporary workspace array.
         */
        if(in_data == out_data)
        {
            ctmp = cdata;
        }
        else
        {
            ctmp = pfft_setup->work2;
        }
                
        /* Transpose from YXZ to YZX. */
        for(j=0;j<local_ny;j++)
        {
            gmx_fft_transpose_2d(cdata + j*nzc*nx,
                                 work  + j*nzc*nx,
                                 nx,
                                 nzc);
        }
        
        /* Perform local_ny*nzc complex FFTs along the x dimension */
        for(i=0;i<local_ny*nzc;i++)
        {
            gmx_fft_1d(pfft_setup->fft_x,
                       GMX_FFT_BACKWARD,
                       work + i*nx,
                       ctmp + i*nx);
        }    
        
        /* Transpose from YZX to YXZ. */
        for(j=0;j<local_ny;j++)
        {
            gmx_fft_transpose_2d(ctmp + j*nzc*nx,
                                 work + j*nzc*nx,
                                 nzc,
                                 nx);
        }
        
        if(in_data == out_data)
        {
            /* output cdata changes when nx or ny not divisible by nnodes */
           ctmp = (t_complex *)in_data + local_x_start*ny*nzc;
        }
        gmx_parallel_transpose_xy(work,
                                  ctmp,
                                  ny,
                                  nx,
                                  pfft_setup->local_slab,
                                  pfft_setup->slab2grid_y,
                                  pfft_setup->slab2grid_x,
                                  nzc,
                                  pfft_setup->nnodes,
                                  pfft_setup->node2slab,
                                  pfft_setup->aav,
                                  pfft_setup->comm);
        
        
        /* Perform nx local 2D complex-to-real FFTs in the yz slices.
         * The 3D FFT is done in-place, so we need to do this in-place too in order
         * to get the data organization right.
         */
        for(i=0;i<local_nx;i++)
        {
            gmx_fft_2d_real(pfft_setup->fft_yz,
                            GMX_FFT_COMPLEX_TO_REAL,
                            ctmp  + i*ny*nzc,
                            rdata + i*ny*nzr);
        }
    }
    else
    {
        gmx_fatal(FARGS,"Incorrect FFT direction.");
    }
    
    /* Skip the YX backtranspose to save communication! Grid is now YXZ */
    return 0;
}




int
gmx_parallel_3dfft_complex2real(gmx_parallel_3dfft_t    pfft_setup,
                                void *                  data)
{
    int          i,j,k;
    int          nx,ny,nzc;
    int          local_x_start,local_nx;
    int          local_y_start,local_ny;    
    t_complex *  work;
    t_complex *  cdata;
    
    work    = pfft_setup->work;
    cdata   = data;

    nx  = pfft_setup->nx;
    ny  = pfft_setup->ny;
    nzc = pfft_setup->nzc;
    
    gmx_parallel_3dfft_limits(pfft_setup,
                              &local_x_start,
                              &local_nx,
                              &local_y_start,
                              &local_ny);

    
    
    return 0;    
}



int
gmx_parallel_3dfft_destroy(gmx_parallel_3dfft_t    pfft_setup)
{
    gmx_fft_destroy(pfft_setup->fft_x);
    gmx_fft_destroy(pfft_setup->fft_yz);

    free(pfft_setup->slab2grid_x);
    free(pfft_setup->slab2grid_y);
    if (pfft_setup->aav) {
        free(pfft_setup->aav->sdisps);
        free(pfft_setup->aav->scounts);
        free(pfft_setup->aav->rdisps);
        free(pfft_setup->aav->rcounts);
        free(pfft_setup->aav);
    }
    free(pfft_setup->work_rawptr);
    free(pfft_setup->work2_rawptr);

    return 0;
}

#else
/* Dummy function to avoid warnings without MPI enabled */
void
gmx_parallel_3dfft_dummy()
{
}

#endif /* GMX_MPI */

