/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
/* IMPORTANT FOR DEVELOPERS:
 *
 * Triclinic pme stuff isn't entirely trivial, and we've experienced
 * some bugs during development (many of them due to me). To avoid
 * this in the future, please check the following things if you make
 * changes in this file:
 *
 * 1. You should obtain identical (at least to the PME precision)
 *    energies, forces, and virial for
 *    a rectangular box and a triclinic one where the z (or y) axis is
 *    tilted a whole box side. For instance you could use these boxes:
 *
 *    rectangular       triclinic
 *     2  0  0           2  0  0
 *     0  2  0           0  2  0
 *     0  0  6           2  2  6
 *
 * 2. You should check the energy conservation in a triclinic box.
 *
 * It might seem an overkill, but better safe than sorry.
 * /Erik 001109
 */ 

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif


#include <stdio.h>
#include <string.h>
#include <math.h>
#include "typedefs.h"
#include "txtdump.h"
#include "vec.h"
#include "gmxcomplex.h"
#include "smalloc.h"
#include "futil.h"
#include "coulomb.h"
#include "gmx_fatal.h"
#include "pme.h"
#include "network.h"
#include "physics.h"
#include "nrnb.h"
#include "copyrite.h"
#include "gmx_wallcycle.h"
#include "gmx_parallel_3dfft.h"
#include "pdbio.h"

#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
#include "gmx_sse2_single.h"
#endif

#include "mpelogging.h"

#define DFT_TOL 1e-7
/* #define PRT_FORCE */
/* conditions for on the fly time-measurement */
/* #define TAKETIME (step > 1 && timesteps < 10) */
#define TAKETIME FALSE

#ifdef GMX_DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif

/* Internal datastructures */
typedef struct {
    int send_index0;
    int send_nindex;
    int recv_index0;
    int recv_nindex;
} pme_grid_comm_t;

typedef struct {
#ifdef GMX_MPI
    MPI_Comm mpi_comm;
#endif
    int  nnodes,nodeid;
    int  *s2g0;
    int  *s2g1;
    int  noverlap_nodes;
    int  *send_id,*recv_id;
    pme_grid_comm_t *comm_data;
} pme_overlap_t;

typedef struct {
    int  dimind;            /* The index of the dimension, 0=x, 1=y */
    int  nslab;
    int  nodeid;
#ifdef GMX_MPI
    MPI_Comm mpi_comm;
#endif

    int  *node_dest;        /* The nodes to send x and q to with DD */
    int  *node_src;         /* The nodes to receive x and q from with DD */
    int  *buf_index;        /* Index for commnode into the buffers */

    int  maxshift;

    int  npd;
    int  pd_nalloc;
    int  *pd;
    int  *count;            /* The number of atoms to send to each node */
    int  *rcount;           /* The number of atoms to receive */

    int  n;
    int  nalloc;
    rvec *x;
    real *q;
    rvec *f;
    gmx_bool bSpread;           /* These coordinates are used for spreading */
    int  pme_order;
    splinevec theta,dtheta;
    ivec *idx;
    rvec *fractx;            /* Fractional coordinate relative to the
                              * lower cell boundary 
                              */
} pme_atomcomm_t;

typedef struct gmx_pme {
    int  ndecompdim;         /* The number of decomposition dimensions */
    int  nodeid;             /* Our nodeid in mpi->mpi_comm */
    int  nodeid_major;
    int  nodeid_minor;
    int  nnodes;             /* The number of nodes doing PME */
    int  nnodes_major;
    int  nnodes_minor;

    MPI_Comm mpi_comm;
    MPI_Comm mpi_comm_d[2];  /* Indexed on dimension, 0=x, 1=y */
#ifdef GMX_MPI
    MPI_Datatype  rvec_mpi;  /* the pme vector's MPI type */
#endif

    gmx_bool bPPnode;            /* Node also does particle-particle forces */
    gmx_bool bFEP;               /* Compute Free energy contribution */
    int nkx,nky,nkz;         /* Grid dimensions */
    int pme_order;
    real epsilon_r;           
    
    real *  pmegridA;  /* Grids on which we do spreading/interpolation, includes overlap */
    real *  pmegridB;
    int     pmegrid_nx,pmegrid_ny,pmegrid_nz;
    int     pmegrid_start_ix,pmegrid_start_iy,pmegrid_start_iz;    
    
    real *  pmegrid_sendbuf;
    real *  pmegrid_recvbuf;
    
    real *fftgridA;             /* Grids for FFT. With 1D FFT decomposition this can be a pointer */
    real *fftgridB;             /* inside the interpolation grid, but separate for 2D PME decomp. */
    int   fftgrid_nx,fftgrid_ny,fftgrid_nz;
    
    t_complex *cfftgridA;             /* Grids for complex FFT data */
    t_complex *cfftgridB;            
    int   cfftgrid_nx,cfftgrid_ny,cfftgrid_nz;
    
    gmx_parallel_3dfft_t  pfft_setupA;
    gmx_parallel_3dfft_t  pfft_setupB;
    
    int  *nnx,*nny,*nnz;
    real *fshx,*fshy,*fshz;
    
    pme_atomcomm_t atc[2];  /* Indexed on decomposition index */
    matrix    recipbox;
    splinevec bsp_mod;
    
    pme_overlap_t overlap[2]; /* Indexed on dimension, 0=x, 1=y */


    pme_atomcomm_t atc_energy; /* Only for gmx_pme_calc_energy */
    
    rvec *bufv;             /* Communication buffer */
    real *bufr;             /* Communication buffer */
    int  buf_nalloc;        /* The communication buffer size */

    /* work data for solve_pme */
    int      work_nalloc;
    real *   work_mhx;
    real *   work_mhy;
    real *   work_mhz;
    real *   work_m2;
    real *   work_denom;
    real *   work_tmp1_alloc;
    real *   work_tmp1;
    real *   work_m2inv;

    /* Work data for PME_redist */
    gmx_bool     redist_init;
    int *    scounts; 
    int *    rcounts;
    int *    sdispls;
    int *    rdispls;
    int *    sidx;
    int *    idxa;    
    real *   redist_buf;
    int      redist_buf_nalloc;
    
    /* Work data for sum_qgrid */
    real *   sum_qgrid_tmp;
    real *   sum_qgrid_dd_tmp;
} t_gmx_pme;


static void calc_interpolation_idx(gmx_pme_t pme,pme_atomcomm_t *atc)
{
    int  i;
    int  *idxptr,tix,tiy,tiz;
    real *xptr,*fptr,tx,ty,tz;
    real rxx,ryx,ryy,rzx,rzy,rzz;
    int  nx,ny,nz;
    int  start_ix,start_iy,start_iz;
    
    nx  = pme->nkx;
    ny  = pme->nky;
    nz  = pme->nkz;
    
    start_ix = pme->pmegrid_start_ix;
    start_iy = pme->pmegrid_start_iy;
    start_iz = pme->pmegrid_start_iz;
    
    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];
    
    for(i=0; (i<atc->n); i++) {
        xptr   = atc->x[i];
        idxptr = atc->idx[i];
        fptr   = atc->fractx[i];
        
        /* Fractional coordinates along box vectors, add 2.0 to make 100% sure we are positive for triclinic boxes */
        tx = nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx + 2.0 );
        ty = ny * (                  xptr[YY] * ryy + xptr[ZZ] * rzy + 2.0 );
        tz = nz * (                                   xptr[ZZ] * rzz + 2.0 );
        
        tix = (int)(tx);
        tiy = (int)(ty);
        tiz = (int)(tz);
        
        /* Because decomposition only occurs in x and y,
         * we never have a fraction correction in z.
         */
        fptr[XX] = tx - tix + pme->fshx[tix];
        fptr[YY] = ty - tiy + pme->fshy[tiy];
        fptr[ZZ] = tz - tiz;   

        idxptr[XX] = pme->nnx[tix];
        idxptr[YY] = pme->nny[tiy];
        idxptr[ZZ] = pme->nnz[tiz];

#ifdef DEBUG
        range_check(idxptr[XX],0,pme->pmegrid_nx);
        range_check(idxptr[YY],0,pme->pmegrid_ny);
        range_check(idxptr[ZZ],0,pme->pmegrid_nz);
#endif
  }  
}

static void pme_calc_pidx(int natoms, matrix recipbox, rvec x[],
                          pme_atomcomm_t *atc)
{
    int  nslab,i;
    int  si;
    real *xptr,s;
    real rxx,ryx,rzx,ryy,rzy;
    int *pd,*count;

    /* Calculate PME task index (pidx) for each grid index.
     * Here we always assign equally sized slabs to each node
     * for load balancing reasons (the PME grid spacing is not used).
     */
    
    nslab = atc->nslab;
    pd    = atc->pd;
    count = atc->count;

    /* Reset the count */
    for(i=0; i<nslab; i++)
    {
        count[i] = 0;
    }
    
    if (atc->dimind == 0)
    {
        rxx = recipbox[XX][XX];
        ryx = recipbox[YY][XX];
        rzx = recipbox[ZZ][XX];
        /* Calculate the node index in x-dimension */
        for(i=0; (i<natoms); i++)
        {
            xptr   = x[i];
            /* Fractional coordinates along box vectors */
            s = nslab*(xptr[XX]*rxx + xptr[YY]*ryx + xptr[ZZ]*rzx);
            si = (int)(s + 2*nslab) % nslab;
            pd[i] = si;
            count[si]++;
        }
    }
    else
    {
        ryy = recipbox[YY][YY];
        rzy = recipbox[ZZ][YY];
        /* Calculate the node index in y-dimension */
        for(i=0; (i<natoms); i++)
        {
            xptr   = x[i];
            /* Fractional coordinates along box vectors */
            s = nslab*(xptr[YY]*ryy + xptr[ZZ]*rzy);
            si = (int)(s + 2*nslab) % nslab;
            pd[i] = si;
            count[si]++;
        }
    }
}

static void pme_realloc_atomcomm_things(pme_atomcomm_t *atc)
{
    int nalloc_old,i;
    
    /* We have to avoid a NULL pointer for atc->x to avoid
     * possible fatal errors in MPI routines.
     */
    if (atc->n > atc->nalloc || atc->nalloc == 0)
    {
        nalloc_old = atc->nalloc;
        atc->nalloc = over_alloc_dd(max(atc->n,1));
        
        if (atc->nslab > 1) {
            srenew(atc->x,atc->nalloc);
            srenew(atc->q,atc->nalloc);
            srenew(atc->f,atc->nalloc);
            for(i=nalloc_old; i<atc->nalloc; i++)
            {
                clear_rvec(atc->f[i]);
            }
        }
        if (atc->bSpread) {
            for(i=0;i<DIM;i++) {
                srenew(atc->theta[i] ,atc->pme_order*atc->nalloc); 
                srenew(atc->dtheta[i],atc->pme_order*atc->nalloc);
            }
            srenew(atc->fractx,atc->nalloc); 
            srenew(atc->idx   ,atc->nalloc);
        }
    }
}

static void pmeredist_pd(gmx_pme_t pme, gmx_bool forw,
                         int n, gmx_bool bXF, rvec *x_f, real *charge,
                         pme_atomcomm_t *atc)
/* Redistribute particle data for PME calculation */
/* domain decomposition by x coordinate           */
{
    int *idxa;
    int i, ii;
    
    if(FALSE == pme->redist_init) {
        snew(pme->scounts,atc->nslab);
        snew(pme->rcounts,atc->nslab);
        snew(pme->sdispls,atc->nslab);
        snew(pme->rdispls,atc->nslab);
        snew(pme->sidx,atc->nslab);
        pme->redist_init = TRUE;
    }
    if (n > pme->redist_buf_nalloc) {
        pme->redist_buf_nalloc = over_alloc_dd(n);
        srenew(pme->redist_buf,pme->redist_buf_nalloc*DIM);
    }
    
    pme->idxa = atc->pd;

#ifdef GMX_MPI
    if (forw && bXF) {
        /* forward, redistribution from pp to pme */ 
        
        /* Calculate send counts and exchange them with other nodes */
        for(i=0; (i<atc->nslab); i++) pme->scounts[i]=0;
        for(i=0; (i<n); i++) pme->scounts[pme->idxa[i]]++;
        MPI_Alltoall( pme->scounts, 1, MPI_INT, pme->rcounts, 1, MPI_INT, atc->mpi_comm);
        
        /* Calculate send and receive displacements and index into send 
           buffer */
        pme->sdispls[0]=0;
        pme->rdispls[0]=0;
        pme->sidx[0]=0;
        for(i=1; i<atc->nslab; i++) {
            pme->sdispls[i]=pme->sdispls[i-1]+pme->scounts[i-1];
            pme->rdispls[i]=pme->rdispls[i-1]+pme->rcounts[i-1];
            pme->sidx[i]=pme->sdispls[i];
        }
        /* Total # of particles to be received */
        atc->n = pme->rdispls[atc->nslab-1] + pme->rcounts[atc->nslab-1];
        
        pme_realloc_atomcomm_things(atc);
        
        /* Copy particle coordinates into send buffer and exchange*/
        for(i=0; (i<n); i++) {
            ii=DIM*pme->sidx[pme->idxa[i]];
            pme->sidx[pme->idxa[i]]++;
            pme->redist_buf[ii+XX]=x_f[i][XX];
            pme->redist_buf[ii+YY]=x_f[i][YY];
            pme->redist_buf[ii+ZZ]=x_f[i][ZZ];
        }
        MPI_Alltoallv(pme->redist_buf, pme->scounts, pme->sdispls, 
                      pme->rvec_mpi, atc->x, pme->rcounts, pme->rdispls, 
                      pme->rvec_mpi, atc->mpi_comm);
    }
    if (forw) {
        /* Copy charge into send buffer and exchange*/
        for(i=0; i<atc->nslab; i++) pme->sidx[i]=pme->sdispls[i];
        for(i=0; (i<n); i++) {
            ii=pme->sidx[pme->idxa[i]];
            pme->sidx[pme->idxa[i]]++;
            pme->redist_buf[ii]=charge[i];
        }
        MPI_Alltoallv(pme->redist_buf, pme->scounts, pme->sdispls, mpi_type,
                      atc->q, pme->rcounts, pme->rdispls, mpi_type,
                      atc->mpi_comm);
    }
    else { /* backward, redistribution from pme to pp */ 
        MPI_Alltoallv(atc->f, pme->rcounts, pme->rdispls, pme->rvec_mpi,
                      pme->redist_buf, pme->scounts, pme->sdispls, 
                      pme->rvec_mpi, atc->mpi_comm);
        
        /* Copy data from receive buffer */
        for(i=0; i<atc->nslab; i++)
            pme->sidx[i] = pme->sdispls[i];
        for(i=0; (i<n); i++) {
            ii = DIM*pme->sidx[pme->idxa[i]];
            x_f[i][XX] += pme->redist_buf[ii+XX];
            x_f[i][YY] += pme->redist_buf[ii+YY];
            x_f[i][ZZ] += pme->redist_buf[ii+ZZ];
            pme->sidx[pme->idxa[i]]++;
        }
    }
#endif 
}

static void pme_dd_sendrecv(pme_atomcomm_t *atc,
                            gmx_bool bBackward,int shift,
                            void *buf_s,int nbyte_s,
                            void *buf_r,int nbyte_r)
{
#ifdef GMX_MPI
    int dest,src;
    MPI_Status stat;
    
    if (bBackward == FALSE) {
        dest = atc->node_dest[shift];
        src  = atc->node_src[shift];
    } else {
        dest = atc->node_src[shift];
        src  = atc->node_dest[shift];
    }
    
    if (nbyte_s > 0 && nbyte_r > 0) {
        MPI_Sendrecv(buf_s,nbyte_s,MPI_BYTE,
                     dest,shift,
                     buf_r,nbyte_r,MPI_BYTE,
                     src,shift,
                     atc->mpi_comm,&stat);
    } else if (nbyte_s > 0) {
        MPI_Send(buf_s,nbyte_s,MPI_BYTE,
                 dest,shift,
                 atc->mpi_comm);
    } else if (nbyte_r > 0) {
        MPI_Recv(buf_r,nbyte_r,MPI_BYTE,
                 src,shift,
                 atc->mpi_comm,&stat);
    }
#endif
}

static void dd_pmeredist_x_q(gmx_pme_t pme, 
                             int n, gmx_bool bX, rvec *x, real *charge,
                             pme_atomcomm_t *atc)
{
    int *commnode,*buf_index;
    int nnodes_comm,i,nsend,local_pos,buf_pos,node,scount,rcount;
    
    commnode  = atc->node_dest;
    buf_index = atc->buf_index;
    
    nnodes_comm = min(2*atc->maxshift,atc->nslab-1);
    
    nsend = 0;
    for(i=0; i<nnodes_comm; i++) {
        buf_index[commnode[i]] = nsend;
        nsend += atc->count[commnode[i]];
    }
    if (bX) {
        if (atc->count[atc->nodeid] + nsend != n)
            gmx_fatal(FARGS,"%d particles communicated to PME node %d are more than 2/3 times the cut-off out of the domain decomposition cell of their charge group in dimension %c.\n"
                      "This usually means that your system is not well equilibrated.",
                      n - (atc->count[atc->nodeid] + nsend),
                      pme->nodeid,'x'+atc->dimind);
        
        if (nsend > pme->buf_nalloc) {
            pme->buf_nalloc = over_alloc_dd(nsend);
            srenew(pme->bufv,pme->buf_nalloc);
            srenew(pme->bufr,pme->buf_nalloc);
        }
        
        atc->n = atc->count[atc->nodeid];
        for(i=0; i<nnodes_comm; i++) {
            scount = atc->count[commnode[i]];
            /* Communicate the count */
            if (debug)
                fprintf(debug,"dimind %d PME node %d send to node %d: %d\n",
                        atc->dimind,atc->nodeid,commnode[i],scount);
            pme_dd_sendrecv(atc,FALSE,i,
                            &scount,sizeof(int),
                            &atc->rcount[i],sizeof(int));
            atc->n += atc->rcount[i];
        }
        
        pme_realloc_atomcomm_things(atc);
    }
    
    local_pos = 0;
    for(i=0; i<n; i++) {
        node = atc->pd[i];
        if (node == atc->nodeid) {
            /* Copy direct to the receive buffer */
            if (bX) {
                copy_rvec(x[i],atc->x[local_pos]);
            }
            atc->q[local_pos] = charge[i];
            local_pos++;
        } else {
            /* Copy to the send buffer */
            if (bX) {
                copy_rvec(x[i],pme->bufv[buf_index[node]]);
            }
            pme->bufr[buf_index[node]] = charge[i];
            buf_index[node]++;
        }
    }
    
    buf_pos = 0;
    for(i=0; i<nnodes_comm; i++) {
        scount = atc->count[commnode[i]];
        rcount = atc->rcount[i];
        if (scount > 0 || rcount > 0) {
            if (bX) {
                /* Communicate the coordinates */
                pme_dd_sendrecv(atc,FALSE,i,
                                pme->bufv[buf_pos],scount*sizeof(rvec),
                                atc->x[local_pos],rcount*sizeof(rvec));
            }
            /* Communicate the charges */
            pme_dd_sendrecv(atc,FALSE,i,
                            pme->bufr+buf_pos,scount*sizeof(real),
                            atc->q+local_pos,rcount*sizeof(real));
            buf_pos   += scount;
            local_pos += atc->rcount[i];
        }
    }
}

static void dd_pmeredist_f(gmx_pme_t pme, pme_atomcomm_t *atc,
                           int n, rvec *f,
                           gmx_bool bAddF)
{
  int *commnode,*buf_index;
  int nnodes_comm,local_pos,buf_pos,i,scount,rcount,node;

  commnode  = atc->node_dest;
  buf_index = atc->buf_index;

  nnodes_comm = min(2*atc->maxshift,atc->nslab-1);

  local_pos = atc->count[atc->nodeid];
  buf_pos = 0;
  for(i=0; i<nnodes_comm; i++) {
    scount = atc->rcount[i];
    rcount = atc->count[commnode[i]];
    if (scount > 0 || rcount > 0) {
      /* Communicate the forces */
      pme_dd_sendrecv(atc,TRUE,i,
                      atc->f[local_pos],scount*sizeof(rvec),
                      pme->bufv[buf_pos],rcount*sizeof(rvec));
      local_pos += scount;
    }
    buf_index[commnode[i]] = buf_pos;
    buf_pos   += rcount;
  }

    local_pos = 0;
    if (bAddF)
    {
        for(i=0; i<n; i++)
        {
            node = atc->pd[i];
            if (node == atc->nodeid)
            {
                /* Add from the local force array */
                rvec_inc(f[i],atc->f[local_pos]);
                local_pos++;
            }
            else
            {
                /* Add from the receive buffer */
                rvec_inc(f[i],pme->bufv[buf_index[node]]);
                buf_index[node]++;
            }
        }
    }
    else
    {
        for(i=0; i<n; i++)
        {
            node = atc->pd[i];
            if (node == atc->nodeid)
            {
                /* Copy from the local force array */
                copy_rvec(atc->f[local_pos],f[i]);
                local_pos++;
            }
            else
            {
                /* Copy from the receive buffer */
                copy_rvec(pme->bufv[buf_index[node]],f[i]);
                buf_index[node]++;
            }
        }
    }
}

#ifdef GMX_MPI
static void 
gmx_sum_qgrid_dd(gmx_pme_t pme, real *grid, int direction)
{
    pme_overlap_t *overlap;
    int send_index0,send_nindex;
    int recv_index0,recv_nindex;
    MPI_Status stat;
    int i,j,k,ix,iy,iz,icnt;
    int ipulse,send_id,recv_id,datasize;
    real *p;
    real *sendptr,*recvptr;
    
    /* Start with minor-rank communication. This is a bit of a pain since it is not contiguous */
    overlap = &pme->overlap[1];
    
    for(ipulse=0;ipulse<overlap->noverlap_nodes;ipulse++)
    {
        /* Since we have already (un)wrapped the overlap in the z-dimension,
         * we only have to communicate 0 to nkz (not pmegrid_nz).
         */
        if (direction==GMX_SUM_QGRID_FORWARD)
        {
            send_id = overlap->send_id[ipulse];
            recv_id = overlap->recv_id[ipulse];
            send_index0   = overlap->comm_data[ipulse].send_index0;
            send_nindex   = overlap->comm_data[ipulse].send_nindex;
            recv_index0   = overlap->comm_data[ipulse].recv_index0;
            recv_nindex   = overlap->comm_data[ipulse].recv_nindex;
        }
        else
        {
            send_id = overlap->recv_id[ipulse];
            recv_id = overlap->send_id[ipulse];
            send_index0   = overlap->comm_data[ipulse].recv_index0;
            send_nindex   = overlap->comm_data[ipulse].recv_nindex;            
            recv_index0   = overlap->comm_data[ipulse].send_index0;
            recv_nindex   = overlap->comm_data[ipulse].send_nindex;
        }

        /* Copy data to contiguous send buffer */
        if (debug)
        {
            fprintf(debug,"PME send node %d %d -> %d grid start %d Communicating %d to %d\n",
                    pme->nodeid,overlap->nodeid,send_id,
                    pme->pmegrid_start_iy,
                    send_index0-pme->pmegrid_start_iy,
                    send_index0-pme->pmegrid_start_iy+send_nindex);
        }
        icnt = 0;
        for(i=0;i<pme->pmegrid_nx;i++)
        {
            ix = i;
            for(j=0;j<send_nindex;j++)
            {
                iy = j + send_index0 - pme->pmegrid_start_iy;
                for(k=0;k<pme->nkz;k++)
                {
                    iz = k;
                    pme->pmegrid_sendbuf[icnt++] = grid[ix*(pme->pmegrid_ny*pme->pmegrid_nz)+iy*(pme->pmegrid_nz)+iz];
                }
            }
        }
            
        datasize      = pme->pmegrid_nx * pme->nkz;
        
        MPI_Sendrecv(pme->pmegrid_sendbuf,send_nindex*datasize,GMX_MPI_REAL,
                     send_id,ipulse,
                     pme->pmegrid_recvbuf,recv_nindex*datasize,GMX_MPI_REAL,
                     recv_id,ipulse,
                     overlap->mpi_comm,&stat);
        
        /* Get data from contiguous recv buffer */
        if (debug)
        {
            fprintf(debug,"PME recv node %d %d <- %d grid start %d Communicating %d to %d\n",
                    pme->nodeid,overlap->nodeid,recv_id,
                    pme->pmegrid_start_iy,
                    recv_index0-pme->pmegrid_start_iy,
                    recv_index0-pme->pmegrid_start_iy+recv_nindex);
        }
        icnt = 0;
        for(i=0;i<pme->pmegrid_nx;i++)
        {
            ix = i;
            for(j=0;j<recv_nindex;j++)
            {
                iy = j + recv_index0 - pme->pmegrid_start_iy;
                for(k=0;k<pme->nkz;k++)
                {
                    iz = k;
                    if(direction==GMX_SUM_QGRID_FORWARD)
                    {
                        grid[ix*(pme->pmegrid_ny*pme->pmegrid_nz)+iy*(pme->pmegrid_nz)+iz] += pme->pmegrid_recvbuf[icnt++];
                    }
                    else
                    {
                        grid[ix*(pme->pmegrid_ny*pme->pmegrid_nz)+iy*(pme->pmegrid_nz)+iz]  = pme->pmegrid_recvbuf[icnt++];
                    }
                }
            }
        }
    }
    
    /* Major dimension is easier, no copying required,
     * but we might have to sum to separate array.
     * Since we don't copy, we have to communicate up to pmegrid_nz,
     * not nkz as for the minor direction.
     */
    overlap = &pme->overlap[0];
    
    for(ipulse=0;ipulse<overlap->noverlap_nodes;ipulse++)
    {
        if(direction==GMX_SUM_QGRID_FORWARD)
        {
            send_id = overlap->send_id[ipulse];
            recv_id = overlap->recv_id[ipulse];
            send_index0   = overlap->comm_data[ipulse].send_index0;
            send_nindex   = overlap->comm_data[ipulse].send_nindex;
            recv_index0   = overlap->comm_data[ipulse].recv_index0;
            recv_nindex   = overlap->comm_data[ipulse].recv_nindex;
            recvptr   = pme->pmegrid_recvbuf;
        }
        else
        {
            send_id = overlap->recv_id[ipulse];
            recv_id = overlap->send_id[ipulse];
            send_index0   = overlap->comm_data[ipulse].recv_index0;
            send_nindex   = overlap->comm_data[ipulse].recv_nindex;            
            recv_index0   = overlap->comm_data[ipulse].send_index0;
            recv_nindex   = overlap->comm_data[ipulse].send_nindex;
            recvptr   = grid + (recv_index0-pme->pmegrid_start_ix)*(pme->pmegrid_ny*pme->pmegrid_nz);
        }
                
        sendptr       = grid + (send_index0-pme->pmegrid_start_ix)*(pme->pmegrid_ny*pme->pmegrid_nz);
        datasize      = pme->pmegrid_ny * pme->pmegrid_nz;

        if (debug)
        {
            fprintf(debug,"PME send node %d %d -> %d grid start %d Communicating %d to %d\n",
                    pme->nodeid,overlap->nodeid,send_id,
                    pme->pmegrid_start_ix,
                    send_index0-pme->pmegrid_start_ix,
                    send_index0-pme->pmegrid_start_ix+send_nindex);
            fprintf(debug,"PME recv node %d %d <- %d grid start %d Communicating %d to %d\n",
                    pme->nodeid,overlap->nodeid,recv_id,
                    pme->pmegrid_start_ix,
                    recv_index0-pme->pmegrid_start_ix,
                    recv_index0-pme->pmegrid_start_ix+recv_nindex);
        }

        MPI_Sendrecv(sendptr,send_nindex*datasize,GMX_MPI_REAL,
                     send_id,ipulse,
                     recvptr,recv_nindex*datasize,GMX_MPI_REAL,
                     recv_id,ipulse,
                     overlap->mpi_comm,&stat);
        
        /* ADD data from contiguous recv buffer */
        if(direction==GMX_SUM_QGRID_FORWARD)
        {        
            p = grid + (recv_index0-pme->pmegrid_start_ix)*(pme->pmegrid_ny*pme->pmegrid_nz);
            for(i=0;i<recv_nindex*datasize;i++)
            {
                p[i] += pme->pmegrid_recvbuf[i];
            }
        }
    }
}
#endif


static int
copy_pmegrid_to_fftgrid(gmx_pme_t pme, real *pmegrid, real *fftgrid)
{
    ivec    local_fft_ndata,local_fft_offset,local_fft_size;
    ivec    local_pme_size;
    int     i,ix,iy,iz;
    int     pmeidx,fftidx;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_real_limits(pme->pfft_setupA,
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);
    
    local_pme_size[0] = pme->pmegrid_nx;
    local_pme_size[1] = pme->pmegrid_ny;
    local_pme_size[2] = pme->pmegrid_nz;
    
    /* The fftgrid is always 'justified' to the lower-left corner of the PME grid, 
     the offset is identical, and the PME grid always has more data (due to overlap)
     */
    {
#ifdef DEBUG_PME
        FILE *fp,*fp2;
        char fn[STRLEN],format[STRLEN];
        real val;
        sprintf(fn,"pmegrid%d.pdb",pme->nodeid);
        fp = ffopen(fn,"w");
        sprintf(fn,"pmegrid%d.txt",pme->nodeid);
        fp2 = ffopen(fn,"w");
     sprintf(format,"%s%s\n",pdbformat,"%6.2f%6.2f");
#endif
    for(ix=0;ix<local_fft_ndata[XX];ix++)
    {
        for(iy=0;iy<local_fft_ndata[YY];iy++)
        {
            for(iz=0;iz<local_fft_ndata[ZZ];iz++)
            {
                pmeidx = ix*(local_pme_size[YY]*local_pme_size[ZZ])+iy*(local_pme_size[ZZ])+iz;
                fftidx = ix*(local_fft_size[YY]*local_fft_size[ZZ])+iy*(local_fft_size[ZZ])+iz;
                fftgrid[fftidx] = pmegrid[pmeidx];
#ifdef DEBUG_PME
                val = 100*pmegrid[pmeidx];
                if (pmegrid[pmeidx] != 0)
                fprintf(fp,format,"ATOM",pmeidx,"CA","GLY",' ',pmeidx,' ',
                        5.0*ix,5.0*iy,5.0*iz,1.0,val);
                if (pmegrid[pmeidx] != 0)
                    fprintf(fp2,"%-12s  %5d  %5d  %5d  %12.5e\n",
                            "qgrid",
                            pme->pmegrid_start_ix + ix,
                            pme->pmegrid_start_iy + iy,
                            pme->pmegrid_start_iz + iz,
                            pmegrid[pmeidx]);
#endif
            }
        }
    }
#ifdef DEBUG_PME
    fclose(fp);
    fclose(fp2);
#endif
    }
    return 0;
}


static int
copy_fftgrid_to_pmegrid(gmx_pme_t pme, real *fftgrid, real *pmegrid)
{
    ivec    local_fft_ndata,local_fft_offset,local_fft_size;
    ivec    local_pme_size;
    int     i,ix,iy,iz;
    int     pmeidx,fftidx;
    
    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_real_limits(pme->pfft_setupA,
                                   local_fft_ndata,
                                   local_fft_offset,
                                   local_fft_size);

    local_pme_size[0] = pme->pmegrid_nx;
    local_pme_size[1] = pme->pmegrid_ny;
    local_pme_size[2] = pme->pmegrid_nz;
    
    /* The fftgrid is always 'justified' to the lower-left corner of the PME grid, 
     the offset is identical, and the PME grid always has more data (due to overlap)
     */
    for(ix=0;ix<local_fft_ndata[XX];ix++)
    {
        for(iy=0;iy<local_fft_ndata[YY];iy++)
        {
            for(iz=0;iz<local_fft_ndata[ZZ];iz++)
            {
                pmeidx = ix*(local_pme_size[YY]*local_pme_size[ZZ])+iy*(local_pme_size[ZZ])+iz;
                fftidx = ix*(local_fft_size[YY]*local_fft_size[ZZ])+iy*(local_fft_size[ZZ])+iz;
                pmegrid[pmeidx] = fftgrid[fftidx];
            }
        }
    }   
    return 0;
}


static void
wrap_periodic_pmegrid(gmx_pme_t pme, real *pmegrid)
{
    int     nx,ny,nz,pnx,pny,pnz,ny_x,overlap,ix,iy,iz;

    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    pnx = pme->pmegrid_nx;
    pny = pme->pmegrid_ny;
    pnz = pme->pmegrid_nz;

    overlap = pme->pme_order - 1;

    /* Add periodic overlap in z */
    for(ix=0; ix<pnx; ix++)
    {
        for(iy=0; iy<pny; iy++)
        {
            for(iz=0; iz<overlap; iz++)
            {
                pmegrid[(ix*pny+iy)*pnz+iz] +=
                    pmegrid[(ix*pny+iy)*pnz+nz+iz];
            }
        }
    }

    if (pme->nnodes_minor == 1)
    {
       for(ix=0; ix<pnx; ix++)
       {
           for(iy=0; iy<overlap; iy++)
           {
               for(iz=0; iz<nz; iz++)
               {
                   pmegrid[(ix*pny+iy)*pnz+iz] +=
                       pmegrid[(ix*pny+ny+iy)*pnz+iz];
               }
           }
       }
    }
     
    if (pme->nnodes_major == 1)
    {
        ny_x = (pme->nnodes_minor == 1 ? ny : pny);

        for(ix=0; ix<overlap; ix++)
        {
            for(iy=0; iy<ny_x; iy++)
            {
                for(iz=0; iz<nz; iz++)
                {
                    pmegrid[(ix*pny+iy)*pnz+iz] +=
                        pmegrid[((nx+ix)*pny+iy)*pnz+iz];
                }
            }
        }
    }
}


static void
unwrap_periodic_pmegrid(gmx_pme_t pme, real *pmegrid)
{
    int     nx,ny,nz,pnx,pny,pnz,ny_x,overlap,ix,iy,iz;

    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;

    pnx = pme->pmegrid_nx;
    pny = pme->pmegrid_ny;
    pnz = pme->pmegrid_nz;

    overlap = pme->pme_order - 1;

    if (pme->nnodes_major == 1)
    {
        ny_x = (pme->nnodes_minor == 1 ? ny : pny);

        for(ix=0; ix<overlap; ix++)
        {
            for(iy=0; iy<ny_x; iy++)
            {
                for(iz=0; iz<nz; iz++)
                {
                    pmegrid[((nx+ix)*pny+iy)*pnz+iz] =
                        pmegrid[(ix*pny+iy)*pnz+iz];
                }
            }
        }
    }

    if (pme->nnodes_minor == 1)
    {
       for(ix=0; ix<pnx; ix++)
       {
           for(iy=0; iy<overlap; iy++)
           {
               for(iz=0; iz<nz; iz++)
               {
                   pmegrid[(ix*pny+ny+iy)*pnz+iz] =
                       pmegrid[(ix*pny+iy)*pnz+iz];
               }
           }
       }
    }

    /* Copy periodic overlap in z */
    for(ix=0; ix<pnx; ix++)
    {
        for(iy=0; iy<pny; iy++)
        {
            for(iz=0; iz<overlap; iz++)
            {
                pmegrid[(ix*pny+iy)*pnz+nz+iz] =
                    pmegrid[(ix*pny+iy)*pnz+iz];
            }
        }
    }
}


/* This has to be a macro to enable full compiler optimization with xlC (and probably others too) */
#define DO_BSPLINE(order)                            \
for(ithx=0; (ithx<order); ithx++)                    \
{                                                    \
    index_x = (i0+ithx)*pny*pnz;                     \
    valx    = qn*thx[ithx];                          \
                                                     \
    for(ithy=0; (ithy<order); ithy++)                \
    {                                                \
        valxy    = valx*thy[ithy];                   \
        index_xy = index_x+(j0+ithy)*pnz;            \
                                                     \
        for(ithz=0; (ithz<order); ithz++)            \
        {                                            \
            index_xyz        = index_xy+(k0+ithz);   \
            grid[index_xyz] += valxy*thz[ithz];      \
        }                                            \
    }                                                \
}


static void spread_q_bsplines(gmx_pme_t pme, pme_atomcomm_t *atc, 
                              real *grid)
{

    /* spread charges from home atoms to local grid */
    pme_overlap_t *ol;
    int      b,i,nn,n,ithx,ithy,ithz,i0,j0,k0;
    int *    idxptr;
    int      order,norder,index_x,index_xy,index_xyz;
    real     valx,valxy,qn;
    real     *thx,*thy,*thz;
    int      localsize, bndsize;
  
    int      pnx,pny,pnz,ndatatot;
  
    pnx = pme->pmegrid_nx;
    pny = pme->pmegrid_ny;
    pnz = pme->pmegrid_nz;
    ndatatot = pnx*pny*pnz;
    
    for(i=0;i<ndatatot;i++)
    {
        grid[i] = 0;
    }

    order = pme->pme_order;

    for(nn=0; (nn<atc->n);nn++) 
    {
        n      = nn;
        qn     = atc->q[n];

        if (qn != 0) 
        {
            idxptr = atc->idx[n];
            norder = n*order;
            
            i0   = idxptr[XX]; 
            j0   = idxptr[YY];
            k0   = idxptr[ZZ];
            thx = atc->theta[XX] + norder;
            thy = atc->theta[YY] + norder;
            thz = atc->theta[ZZ] + norder;
            
            switch (order) {
            case 4:  DO_BSPLINE(4);     break;
            case 5:  DO_BSPLINE(5);     break;
            default: DO_BSPLINE(order); break;
            }
        }
    }	
}


#if ( !defined(GMX_DOUBLE) && ( defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE) || defined(GMX_X86_64_SSE2) ) )
    /* Calculate exponentials through SSE in float precision */
#define CALC_EXPONENTIALS(start,end,r_aligned)      \
    {                                               \
        __m128 tmp_sse;                             \
        for(kx=0; kx<end; kx+=4)                    \
        {                                           \
            tmp_sse = _mm_load_ps(r_aligned+kx);    \
            tmp_sse = gmx_mm_exp_ps(tmp_sse);       \
            _mm_store_ps(r_aligned+kx,tmp_sse);     \
        }                                           \
    }
#else
#define CALC_EXPONENTIALS(start,end,r)          \
    for(kx=start; kx<end; kx++)                 \
    {                                           \
        r[kx] = exp(r[kx]);                     \
    }
#endif


static int solve_pme_yzx(gmx_pme_t pme,t_complex *grid,
                         real ewaldcoeff,real vol,
                         gmx_bool bEnerVir,real *mesh_energy,matrix vir)
{
    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    t_complex *p0;
    int     kx,ky,kz,maxkx,maxky,maxkz;
    int     nx,ny,nz,iy,iz,kxstart,kxend;
    real    mx,my,mz;
    real    factor=M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    real    ets2,struct2,vfactor,ets2vf;
    real    eterm,d1,d2,energy=0;
    real    by,bz;
    real    virxx=0,virxy=0,virxz=0,viryy=0,viryz=0,virzz=0;
    real    rxx,ryx,ryy,rzx,rzy,rzz;
	real    *mhx,*mhy,*mhz,*m2,*denom,*tmp1,*m2inv;
    real    mhxk,mhyk,mhzk,m2k;
    real    corner_fac;
    ivec    complex_order;
    ivec    local_ndata,local_offset,local_size;
    
    nx = pme->nkx;
    ny = pme->nky;
    nz = pme->nkz;
    
    /* Dimensions should be identical for A/B grid, so we just use A here */
    gmx_parallel_3dfft_complex_limits(pme->pfft_setupA,
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);
    
    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];
    
    maxkx = (nx+1)/2;
    maxky = (ny+1)/2;
    maxkz = nz/2+1;
	
	mhx   = pme->work_mhx;
	mhy   = pme->work_mhy;
	mhz   = pme->work_mhz;
	m2    = pme->work_m2;
	denom = pme->work_denom;
	tmp1  = pme->work_tmp1;
	m2inv = pme->work_m2inv;	

    for(iy=0;iy<local_ndata[YY];iy++)
    {
        ky = iy + local_offset[YY];
        
        if (ky < maxky) 
        {
            my = ky;
        }
        else 
        {
            my = (ky - ny);
        }
        
        by = M_PI*vol*pme->bsp_mod[YY][ky];

        for(iz=0;iz<local_ndata[ZZ];iz++)
        {
            kz = iz + local_offset[ZZ];
            
            mz = kz;

            bz = pme->bsp_mod[ZZ][kz];
            
            /* 0.5 correction for corner points */
			corner_fac = 1;
            if (kz == 0)
                corner_fac = 0.5;
            if (kz == (nz+1)/2)
                corner_fac = 0.5;
                      
            p0 = grid + iy*local_size[ZZ]*local_size[XX] + iz*local_size[XX];
            
            /* We should skip the k-space point (0,0,0) */
            if (local_offset[XX] > 0 ||
                local_offset[YY] > 0 || ky > 0 ||
                kz > 0)
            {
                kxstart = local_offset[XX];
            }
            else
            {
                kxstart = local_offset[XX] + 1;
                p0++;
            }
            kxend = local_offset[XX] + local_ndata[XX];
			
            if (bEnerVir)
            {
                /* More expensive inner loop, especially because of the storage
                 * of the mh elements in array's.
                 * Because x is the minor grid index, all mh elements
                 * depend on kx for triclinic unit cells.
                 */

                /* Two explicit loops to avoid a conditional inside the loop */
                for(kx=kxstart; kx<maxkx; kx++)
                {
                    mx = kx;
                    
                    mhxk      = mx * rxx;
                    mhyk      = mx * ryx + my * ryy;
                    mhzk      = mx * rzx + my * rzy + mz * rzz;
                    m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                    mhx[kx]   = mhxk;
                    mhy[kx]   = mhyk;
                    mhz[kx]   = mhzk;
                    m2[kx]    = m2k;
                    denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                    tmp1[kx]  = -factor*m2k;
                }
                
                for(kx=maxkx; kx<kxend; kx++)
                {
                    mx = (kx - nx);

                    mhxk      = mx * rxx;
                    mhyk      = mx * ryx + my * ryy;
                    mhzk      = mx * rzx + my * rzy + mz * rzz;
                    m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                    mhx[kx]   = mhxk;
                    mhy[kx]   = mhyk;
                    mhz[kx]   = mhzk;
                    m2[kx]    = m2k;
                    denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                    tmp1[kx]  = -factor*m2k;
                }
                
                for(kx=kxstart; kx<kxend; kx++)
                {
                    m2inv[kx] = 1.0/m2[kx];
                }
                for(kx=kxstart; kx<kxend; kx++)
                {
                    denom[kx] = 1.0/denom[kx];
                }

                CALC_EXPONENTIALS(kxstart,kxend,tmp1);

                for(kx=kxstart; kx<kxend; kx++,p0++)
                {
                    d1      = p0->re;
                    d2      = p0->im;
                    
                    eterm    = ONE_4PI_EPS0/pme->epsilon_r*tmp1[kx]*denom[kx];
                    
                    p0->re  = d1*eterm;
                    p0->im  = d2*eterm;
                    
                    struct2 = 2.0*(d1*d1+d2*d2);
                    
                    tmp1[kx] = eterm*struct2;
                }
                
                for(kx=kxstart; kx<kxend; kx++)
                {
                    ets2     = corner_fac*tmp1[kx];
                    vfactor  = (factor*m2[kx] + 1.0)*2.0*m2inv[kx];
                    energy  += ets2;
                    
                    ets2vf   = ets2*vfactor;
                    virxx   += ets2vf*mhx[kx]*mhx[kx] - ets2;
                    virxy   += ets2vf*mhx[kx]*mhy[kx];
                    virxz   += ets2vf*mhx[kx]*mhz[kx];
                    viryy   += ets2vf*mhy[kx]*mhy[kx] - ets2;
                    viryz   += ets2vf*mhy[kx]*mhz[kx];
                    virzz   += ets2vf*mhz[kx]*mhz[kx] - ets2;
                }
            }
            else
            {
                /* We don't need to calculate the energy and the virial.
                 * In this case the triclinic overhead is small.
                 */

                /* Two explicit loops to avoid a conditional inside the loop */

                for(kx=kxstart; kx<maxkx; kx++)
                {
                    mx = kx;
                    
                    mhxk      = mx * rxx;
                    mhyk      = mx * ryx + my * ryy;
                    mhzk      = mx * rzx + my * rzy + mz * rzz;
                    m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                    denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                    tmp1[kx]  = -factor*m2k;
                }
                
                for(kx=maxkx; kx<kxend; kx++)
                {
                    mx = (kx - nx);
                    
                    mhxk      = mx * rxx;
                    mhyk      = mx * ryx + my * ryy;
                    mhzk      = mx * rzx + my * rzy + mz * rzz;
                    m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                    denom[kx] = m2k*bz*by*pme->bsp_mod[XX][kx];
                    tmp1[kx]  = -factor*m2k;
                }
                
                for(kx=kxstart; kx<kxend; kx++)
                {
                    denom[kx] = 1.0/denom[kx];
                }

                CALC_EXPONENTIALS(kxstart,kxend,tmp1);
               
                for(kx=kxstart; kx<kxend; kx++,p0++)
                {
                    d1      = p0->re;
                    d2      = p0->im;
                    
                    eterm    = ONE_4PI_EPS0/pme->epsilon_r*tmp1[kx]*denom[kx];
                    
                    p0->re  = d1*eterm;
                    p0->im  = d2*eterm;
                }
            }
        }
    }
    
    if (bEnerVir)
    {
        /* Update virial with local values.
         * The virial is symmetric by definition.
         * this virial seems ok for isotropic scaling, but I'm
         * experiencing problems on semiisotropic membranes.
         * IS THAT COMMENT STILL VALID??? (DvdS, 2001/02/07).
         */
        vir[XX][XX] = 0.25*virxx;
        vir[YY][YY] = 0.25*viryy;
        vir[ZZ][ZZ] = 0.25*virzz;
        vir[XX][YY] = vir[YY][XX] = 0.25*virxy;
        vir[XX][ZZ] = vir[ZZ][XX] = 0.25*virxz;
        vir[YY][ZZ] = vir[ZZ][YY] = 0.25*viryz;
        
        /* This energy should be corrected for a charged system */
        *mesh_energy = 0.5*energy;
    }

    /* Return the loop count */
    return local_ndata[YY]*local_ndata[ZZ]*local_ndata[XX];
}


#define DO_FSPLINE(order)                      \
for(ithx=0; (ithx<order); ithx++)              \
{           								   \
    index_x = (i0+ithx)*pny*pnz;               \
    tx      = thx[ithx];                       \
    dx      = dthx[ithx];                      \
                                               \
    for(ithy=0; (ithy<order); ithy++)          \
    {										   \
        index_xy = index_x+(j0+ithy)*pnz;      \
        ty       = thy[ithy];                  \
        dy       = dthy[ithy];                 \
        fxy1     = fz1 = 0;                    \
                                               \
        for(ithz=0; (ithz<order); ithz++)      \
        {     								   \
            gval  = grid[index_xy+(k0+ithz)];  \
            fxy1 += thz[ithz]*gval;            \
            fz1  += dthz[ithz]*gval;           \
        }                                      \
        fx += dx*ty*fxy1;                      \
        fy += tx*dy*fxy1;                      \
        fz += tx*ty*fz1;                       \
    }                                          \
}


void gather_f_bsplines(gmx_pme_t pme,real *grid,
                       gmx_bool bClearF,pme_atomcomm_t *atc,real scale)
{
    /* sum forces for local particles */  
    int     nn,n,ithx,ithy,ithz,i0,j0,k0;
    int     index_x,index_xy;
    int     nx,ny,nz,pnx,pny,pnz;
    int *   idxptr;
    real    tx,ty,dx,dy,qn;
    real    fx,fy,fz,gval;
    real    fxy1,fz1;
    real    *thx,*thy,*thz,*dthx,*dthy,*dthz;
    int     norder;
    real    rxx,ryx,ryy,rzx,rzy,rzz;
    int     order;
    
    order = pme->pme_order;
    thx   = atc->theta[XX];
    thy   = atc->theta[YY];
    thz   = atc->theta[ZZ];
    dthx  = atc->dtheta[XX];
    dthy  = atc->dtheta[YY];
    dthz  = atc->dtheta[ZZ];
    nx    = pme->nkx;
    ny    = pme->nky;
    nz    = pme->nkz;
    pnx   = pme->pmegrid_nx;
    pny   = pme->pmegrid_ny;
    pnz   = pme->pmegrid_nz;
    
    rxx   = pme->recipbox[XX][XX];
    ryx   = pme->recipbox[YY][XX];
    ryy   = pme->recipbox[YY][YY];
    rzx   = pme->recipbox[ZZ][XX];
    rzy   = pme->recipbox[ZZ][YY];
    rzz   = pme->recipbox[ZZ][ZZ];

    for(nn=0; (nn<atc->n); nn++) 
    {
        n = nn;
        qn      = scale*atc->q[n];
        
        if (bClearF) 
        {
            atc->f[n][XX] = 0;
            atc->f[n][YY] = 0;
            atc->f[n][ZZ] = 0;
        }
        if (qn != 0) 
        {
            fx     = 0;
            fy     = 0;
            fz     = 0;
            idxptr = atc->idx[n];
            norder = n*order;
            
            i0   = idxptr[XX]; 
            j0   = idxptr[YY];
            k0   = idxptr[ZZ];
            
            /* Pointer arithmetic alert, next six statements */
            thx  = atc->theta[XX] + norder;
            thy  = atc->theta[YY] + norder;
            thz  = atc->theta[ZZ] + norder;
            dthx = atc->dtheta[XX] + norder;
            dthy = atc->dtheta[YY] + norder;
            dthz = atc->dtheta[ZZ] + norder;
            
            switch (order) {
            case 4:  DO_FSPLINE(4);     break;
            case 5:  DO_FSPLINE(5);     break;
            default: DO_FSPLINE(order); break;
            }

            atc->f[n][XX] += -qn*( fx*nx*rxx );
            atc->f[n][YY] += -qn*( fx*nx*ryx + fy*ny*ryy );
            atc->f[n][ZZ] += -qn*( fx*nx*rzx + fy*ny*rzy + fz*nz*rzz );
        }
    }
    /* Since the energy and not forces are interpolated
     * the net force might not be exactly zero.
     * This can be solved by also interpolating F, but
     * that comes at a cost.
     * A better hack is to remove the net force every
     * step, but that must be done at a higher level
     * since this routine doesn't see all atoms if running
     * in parallel. Don't know how important it is?  EL 990726
     */
}

static real gather_energy_bsplines(gmx_pme_t pme,real *grid,
                                   pme_atomcomm_t *atc)
{
    int     n,ithx,ithy,ithz,i0,j0,k0;
    int     index_x,index_xy;
    int *   idxptr;
    real    energy,pot,tx,ty,qn,gval;
    real    *thx,*thy,*thz;
    int     norder;
    int     order;
    
    
    order = pme->pme_order;
    
    energy = 0;
    for(n=0; (n<atc->n); n++) {
        qn      = atc->q[n];
        
        if (qn != 0) {
            idxptr = atc->idx[n];
            norder = n*order;
            
            i0   = idxptr[XX]; 
            j0   = idxptr[YY];
            k0   = idxptr[ZZ];
            
            /* Pointer arithmetic alert, next three statements */
            thx  = atc->theta[XX] + norder;
            thy  = atc->theta[YY] + norder;
            thz  = atc->theta[ZZ] + norder;

            pot = 0;
            for(ithx=0; (ithx<order); ithx++)
            {
                index_x = (i0+ithx)*pme->pmegrid_ny*pme->pmegrid_nz;
                tx      = thx[ithx];

                for(ithy=0; (ithy<order); ithy++)
                {
                    index_xy = index_x+(j0+ithy)*pme->pmegrid_nz;
                    ty       = thy[ithy];

                    for(ithz=0; (ithz<order); ithz++)
                    {
                        gval  = grid[index_xy+(k0+ithz)];
                        pot  += tx*ty*thz[ithz]*gval;
                    }

                }
            }

            energy += pot*qn;
        }
    }

    return energy;
}

void make_bsplines(splinevec theta,splinevec dtheta,int order,
                   rvec fractx[],int nr,real charge[],
                   gmx_bool bFreeEnergy)
{
    /* construct splines for local atoms */
    int  i,j,k,l;
    real dr,div;
    real *data,*ddata,*xptr;
    
    for(i=0; (i<nr); i++) {
        /* With free energy we do not use the charge check.
         * In most cases this will be more efficient than calling make_bsplines
         * twice, since usually more than half the particles have charges.
         */
        if (bFreeEnergy || charge[i] != 0.0) {
            xptr = fractx[i];
            for(j=0; (j<DIM); j++) {
                dr  = xptr[j];
                
                /* dr is relative offset from lower cell limit */
                data=&(theta[j][i*order]);
                data[order-1]=0;
                data[1]=dr;
                data[0]=1-dr;
                
                for(k=3; (k<order); k++) {
                    div=1.0/(k-1.0);    
                    data[k-1]=div*dr*data[k-2];
                    for(l=1; (l<(k-1)); l++) {
                        data[k-l-1]=div*((dr+l)*data[k-l-2]+(k-l-dr)*
                                         data[k-l-1]);
                    }
                    data[0]=div*(1-dr)*data[0];
                }
                /* differentiate */
                ddata    = &(dtheta[j][i*order]);
                ddata[0] = -data[0];
                for(k=1; (k<order); k++) {
                    ddata[k]=data[k-1]-data[k];
                }
                
                div=1.0/(order-1);
                data[order-1]=div*dr*data[order-2];
                for(l=1; (l<(order-1)); l++) {
                    data[order-l-1]=div*((dr+l)*data[order-l-2]+
                                         (order-l-dr)*data[order-l-1]);
                }
                data[0]=div*(1-dr)*data[0]; 
            }
        }
    }
}


void make_dft_mod(real *mod,real *data,int ndata)
{
  int i,j;
  real sc,ss,arg;
    
  for(i=0;i<ndata;i++) {
    sc=ss=0;
    for(j=0;j<ndata;j++) {
      arg=(2.0*M_PI*i*j)/ndata;
      sc+=data[j]*cos(arg);
      ss+=data[j]*sin(arg);
    }
    mod[i]=sc*sc+ss*ss;
  }
  for(i=0;i<ndata;i++)
    if(mod[i]<1e-7)
      mod[i]=(mod[i-1]+mod[i+1])*0.5;
}



void make_bspline_moduli(splinevec bsp_mod,int nx,int ny,int nz,int order)
{
  int nmax=max(nx,max(ny,nz));
  real *data,*ddata,*bsp_data;
  int i,k,l;
  real div;
    
  snew(data,order);
  snew(ddata,order);
  snew(bsp_data,nmax);

  data[order-1]=0;
  data[1]=0;
  data[0]=1;
	    
  for(k=3;k<order;k++) {
    div=1.0/(k-1.0);
    data[k-1]=0;
    for(l=1;l<(k-1);l++)
      data[k-l-1]=div*(l*data[k-l-2]+(k-l)*data[k-l-1]);
    data[0]=div*data[0];
  }
  /* differentiate */
  ddata[0]=-data[0];
  for(k=1;k<order;k++)
    ddata[k]=data[k-1]-data[k];
  div=1.0/(order-1);
  data[order-1]=0;
  for(l=1;l<(order-1);l++)
    data[order-l-1]=div*(l*data[order-l-2]+(order-l)*data[order-l-1]);
  data[0]=div*data[0]; 

  for(i=0;i<nmax;i++)
    bsp_data[i]=0;
  for(i=1;i<=order;i++)
    bsp_data[i]=data[i-1];
    
  make_dft_mod(bsp_mod[XX],bsp_data,nx);
  make_dft_mod(bsp_mod[YY],bsp_data,ny);
  make_dft_mod(bsp_mod[ZZ],bsp_data,nz);

  sfree(data);
  sfree(ddata);
  sfree(bsp_data);
}

static void setup_coordinate_communication(pme_atomcomm_t *atc)
{
  int nslab,n,i;
  int fw,bw;

  nslab = atc->nslab;

  n = 0;
  for(i=1; i<=nslab/2; i++) {
    fw = (atc->nodeid + i) % nslab;
    bw = (atc->nodeid - i + nslab) % nslab;
    if (n < nslab - 1) {
      atc->node_dest[n] = fw;
      atc->node_src[n]  = bw;
      n++;
    } 
    if (n < nslab - 1) {
      atc->node_dest[n] = bw;
      atc->node_src[n]  = fw;
      n++;
    }
  }
}

int gmx_pme_destroy(FILE *log,gmx_pme_t *pmedata)
{
    if(NULL != log)
    {
        fprintf(log,"Destroying PME data structures.\n");
    }

    sfree((*pmedata)->nnx);
    sfree((*pmedata)->nny);
    sfree((*pmedata)->nnz);
	
    sfree((*pmedata)->pmegridA);
    sfree((*pmedata)->fftgridA);
    sfree((*pmedata)->cfftgridA);
    gmx_parallel_3dfft_destroy((*pmedata)->pfft_setupA);
    
    if((*pmedata)->pmegridB)
    {
        sfree((*pmedata)->pmegridB);
        sfree((*pmedata)->fftgridB);
        sfree((*pmedata)->cfftgridB);
        gmx_parallel_3dfft_destroy((*pmedata)->pfft_setupB);
    }
    sfree((*pmedata)->work_mhz);
    sfree((*pmedata)->work_m2);
    sfree((*pmedata)->work_denom);
    sfree((*pmedata)->work_tmp1_alloc);
    sfree((*pmedata)->work_m2inv);
	
    sfree(*pmedata);
    *pmedata = NULL;
  
  return 0;
}

static int mult_up(int n,int f)
{
    return ((n + f - 1)/f)*f;
}


static double pme_load_imbalance(gmx_pme_t pme)
{
    int    nma,nmi;
    double n1,n2,n3;

    nma = pme->nnodes_major;
    nmi = pme->nnodes_minor;

    n1 = mult_up(pme->nkx,nma)*mult_up(pme->nky,nmi)*pme->nkz;
    n2 = mult_up(pme->nkx,nma)*mult_up(pme->nkz,nmi)*pme->nky;
    n3 = mult_up(pme->nky,nma)*mult_up(pme->nkz,nmi)*pme->nkx;

    /* pme_solve is roughly double the cost of an fft */

    return (n1 + n2 + 3*n3)/(double)(6*pme->nkx*pme->nky*pme->nkz);
}

static void init_atomcomm(gmx_pme_t pme,pme_atomcomm_t *atc, t_commrec *cr,
                          int dimind,gmx_bool bSpread)
{
    int nk,k,s;

    atc->dimind = dimind;
    atc->nslab  = 1;
    atc->nodeid = 0;
    atc->pd_nalloc = 0;
#ifdef GMX_MPI
    if (pme->nnodes > 1)
    {
        atc->mpi_comm = pme->mpi_comm_d[dimind];
        MPI_Comm_size(atc->mpi_comm,&atc->nslab);
        MPI_Comm_rank(atc->mpi_comm,&atc->nodeid);
    }
    if (debug)
    {
        fprintf(debug,"For PME atom communication in dimind %d: nslab %d rank %d\n",atc->dimind,atc->nslab,atc->nodeid);
    }
#endif

    atc->bSpread   = bSpread;
    atc->pme_order = pme->pme_order;

    if (atc->nslab > 1)
    {
        /* These three allocations are not required for particle decomp. */
        snew(atc->node_dest,atc->nslab);
        snew(atc->node_src,atc->nslab);
        setup_coordinate_communication(atc);
        
        snew(atc->count,atc->nslab);
        snew(atc->rcount,atc->nslab);
        snew(atc->buf_index,atc->nslab);
    }
}

static void 
init_overlap_comm(pme_overlap_t *  ol,
                  int              norder,
#ifdef GMX_MPI
                  MPI_Comm         comm,  
#endif
                  int              nnodes, 
                  int              nodeid,
                  int              ndata)
{
    int lbnd,rbnd,maxlr,b,i;
    int exten;
    int nn,nk;
    pme_grid_comm_t *pgc;
    gmx_bool bCont;
    int fft_start,fft_end,send_index1,recv_index1;
    
#ifdef GMX_MPI
    ol->mpi_comm = comm;
#endif
    
    ol->nnodes = nnodes;
    ol->nodeid = nodeid;

    /* Linear translation of the PME grid wo'nt affect reciprocal space
     * calculations, so to optimize we only interpolate "upwards",
     * which also means we only have to consider overlap in one direction.
     * I.e., particles on this node might also be spread to grid indices
     * that belong to higher nodes (modulo nnodes)
     */

    snew(ol->s2g0,ol->nnodes+1);
    snew(ol->s2g1,ol->nnodes);
    if (debug) { fprintf(debug,"PME slab boundaries:"); }
    for(i=0; i<nnodes; i++) 
    {
        /* s2g0 the local interpolation grid start.
         * s2g1 the local interpolation grid end.
         * Because grid overlap communication only goes forward,
         * the grid the slabs for fft's should be rounded down.
         */
        ol->s2g0[i] = ( i   *ndata + 0       )/nnodes;
        ol->s2g1[i] = ((i+1)*ndata + nnodes-1)/nnodes + norder - 1;

        if (debug)
        {
            fprintf(debug,"  %3d %3d",ol->s2g0[i],ol->s2g1[i]);
        }
    }
    ol->s2g0[nnodes] = ndata;
    if (debug) { fprintf(debug,"\n"); }

    /* Determine with how many nodes we need to communicate the grid overlap */
    b = 0;
    do
    {
        b++;
        bCont = FALSE;
        for(i=0; i<nnodes; i++)
        {
            if ((i+b <  nnodes && ol->s2g1[i] > ol->s2g0[i+b]) ||
                (i+b >= nnodes && ol->s2g1[i] > ol->s2g0[i+b-nnodes] + ndata))
            {
                bCont = TRUE;
            }
        }
    }
    while (bCont && b < nnodes);
    ol->noverlap_nodes = b - 1;

    snew(ol->send_id,ol->noverlap_nodes);
    snew(ol->recv_id,ol->noverlap_nodes);
    for(b=0; b<ol->noverlap_nodes; b++)
    {
        ol->send_id[b] = (ol->nodeid + (b + 1)) % ol->nnodes;
        ol->recv_id[b] = (ol->nodeid - (b + 1) + ol->nnodes) % ol->nnodes;
    }
    snew(ol->comm_data, ol->noverlap_nodes);

    for(b=0; b<ol->noverlap_nodes; b++)
    {
        pgc = &ol->comm_data[b];
        /* Send */
        fft_start        = ol->s2g0[ol->send_id[b]];
        fft_end          = ol->s2g0[ol->send_id[b]+1];
        if (ol->send_id[b] < nodeid)
        {
            fft_start += ndata;
            fft_end   += ndata;
        }
        send_index1      = ol->s2g1[nodeid];
        send_index1      = min(send_index1,fft_end);
        pgc->send_index0 = fft_start;
        pgc->send_nindex = max(0,send_index1 - pgc->send_index0);

        /* We always start receiving to the first index of our slab */
        fft_start        = ol->s2g0[ol->nodeid];
        fft_end          = ol->s2g0[ol->nodeid+1];
        recv_index1      = ol->s2g1[ol->recv_id[b]];
        if (ol->recv_id[b] > nodeid)
        {
            recv_index1 -= ndata;
        }
        recv_index1      = min(recv_index1,fft_end);
        pgc->recv_index0 = fft_start;
        pgc->recv_nindex = max(0,recv_index1 - pgc->recv_index0);
    }
}

static void
make_gridindex5_to_localindex(int n,int local_start,int local_range,
                              int **global_to_local,
                              real **fraction_shift)
{
    int i;
    int * gtl;
    real * fsh;

    snew(gtl,5*n);
    snew(fsh,5*n);
    for(i=0; (i<5*n); i++)
    {
        /* Determine the global to local grid index */
        gtl[i] = (i - local_start + n) % n;
        /* For coordinates that fall within the local grid the fraction
         * is correct, we don't need to shift it.
         */
        fsh[i] = 0;
        if (local_range < n)
        {
            /* Due to rounding issues i could be 1 beyond the lower or
             * upper boundary of the local grid. Correct the index for this.
             * If we shift the index, we need to shift the fraction by
             * the same amount in the other direction to not affect
             * the weights.
             * Note that due to this shifting the weights at the end of
             * the spline might change, but that will only involve values
             * between zero and values close to the precision of a real,
             * which is anyhow the accuracy of the whole mesh calculation.
             */
            /* With local_range=0 we should not change i=local_start */
            if (i % n != local_start)
            {
                if (gtl[i] == n-1)
                {
                    gtl[i] = 0;
                    fsh[i] = -1; 
                }
                else if (gtl[i] == local_range)
                {
                    gtl[i] = local_range - 1;
                    fsh[i] = 1;
                }
            }
        }
    }

    *global_to_local = gtl;
    *fraction_shift  = fsh;
}

static void
gmx_pme_check_grid_restrictions(FILE *fplog,char dim,int nnodes,int *nk)
{
    int nk_new;

    if (*nk % nnodes != 0)
    {
        nk_new = nnodes*(*nk/nnodes + 1);

        if (2*nk_new >= 3*(*nk))
        {
            gmx_fatal(FARGS,"The PME grid size in dim %c (%d) is not divisble by the number of nodes doing PME in dim %c (%d). The grid size would have to be increased by more than 50%% to make the grid divisible. Change the total number of nodes or the number of domain decomposition cells in x or the PME grid %c dimension (and the cut-off).",
                      dim,*nk,dim,nnodes,dim);
        }
        
        if (fplog != NULL)
        {
            fprintf(fplog,"\nNOTE: The PME grid size in dim %c (%d) is not divisble by the number of nodes doing PME in dim %c (%d). Increasing the PME grid size in dim %c to %d. This will increase the accuracy and will not decrease the performance significantly on this number of nodes. For optimal performance change the total number of nodes or the number of domain decomposition cells in x or the PME grid %c dimension (and the cut-off).\n\n",
                    dim,*nk,dim,nnodes,dim,nk_new,dim);
        }
            
        *nk = nk_new;
    }
}

int gmx_pme_init(gmx_pme_t *         pmedata,
                 t_commrec *         cr,
                 int                 nnodes_major,
                 int                 nnodes_minor,
                 t_inputrec *        ir,
                 int                 homenr,
                 gmx_bool                bFreeEnergy,
                 gmx_bool                bReproducible)
{
    gmx_pme_t pme=NULL;
    
    pme_atomcomm_t *atc;
    int bufsizex,bufsizey,bufsize;
    ivec ndata;
    
    if (debug)
        fprintf(debug,"Creating PME data structures.\n");
    snew(pme,1);
        
    pme->redist_init         = FALSE;
    pme->sum_qgrid_tmp       = NULL;
    pme->sum_qgrid_dd_tmp    = NULL;
    pme->buf_nalloc          = 0;
    pme->redist_buf_nalloc   = 0;
    
    pme->nnodes              = 1;
    pme->bPPnode             = TRUE;
    
    pme->nnodes_major        = nnodes_major;
    pme->nnodes_minor        = nnodes_minor;

#ifdef GMX_MPI
    if (nnodes_major*nnodes_minor > 1 && PAR(cr)) 
    {
        pme->mpi_comm = cr->mpi_comm_mygroup;
        
        MPI_Comm_rank(pme->mpi_comm,&pme->nodeid);
        MPI_Comm_size(pme->mpi_comm,&pme->nnodes);
        if (pme->nnodes != nnodes_major*nnodes_minor)
        {
            gmx_incons("PME node count mismatch");
        }
    }
    else
    {
        pme->mpi_comm = MPI_COMM_NULL;
    }
#endif

    if (pme->nnodes == 1)
    {
#ifdef GMX_MPI
        pme->mpi_comm_d[0] = MPI_COMM_NULL;
        pme->mpi_comm_d[1] = MPI_COMM_NULL;
#endif
        pme->ndecompdim = 0;
        pme->nodeid_major = 0;
        pme->nodeid_minor = 0;
    }
    else
    {
        if (nnodes_minor == 1)
        {
#ifdef GMX_MPI
            pme->mpi_comm_d[0] = pme->mpi_comm;
            pme->mpi_comm_d[1] = MPI_COMM_NULL;
#endif
            pme->ndecompdim = 1;
            pme->nodeid_major = pme->nodeid;
            pme->nodeid_minor = 0;
            
        }
        else if (nnodes_major == 1)
        {
#ifdef GMX_MPI
            pme->mpi_comm_d[0] = MPI_COMM_NULL;
            pme->mpi_comm_d[1] = pme->mpi_comm;
#endif
            pme->ndecompdim = 1;
            pme->nodeid_major = 0;
            pme->nodeid_minor = pme->nodeid;
        }
        else 
        {
            if (pme->nnodes % nnodes_major != 0)
            {
                gmx_incons("For 2D PME decomposition, #PME nodes must be divisible by the number of nodes in the major dimension");
            }
            pme->ndecompdim = 2;
            
#ifdef GMX_MPI
            MPI_Comm_split(pme->mpi_comm,pme->nodeid % nnodes_minor,
                           pme->nodeid,&pme->mpi_comm_d[0]);  /* My communicator along major dimension */
            MPI_Comm_split(pme->mpi_comm,pme->nodeid/nnodes_minor,
                           pme->nodeid,&pme->mpi_comm_d[1]);  /* My communicator along minor dimension */
            
            MPI_Comm_rank(pme->mpi_comm_d[0],&pme->nodeid_major);
            MPI_Comm_size(pme->mpi_comm_d[0],&pme->nnodes_major);
            MPI_Comm_rank(pme->mpi_comm_d[1],&pme->nodeid_minor);
            MPI_Comm_size(pme->mpi_comm_d[1],&pme->nnodes_minor);
#endif
        }
        pme->bPPnode = (cr->duty & DUTY_PP);
    }
    
    if (ir->ePBC == epbcSCREW)
    {
        gmx_fatal(FARGS,"pme does not (yet) work with pbc = screw");
    }
    
    pme->bFEP        = ((ir->efep != efepNO) && bFreeEnergy);
    pme->nkx         = ir->nkx;
    pme->nky         = ir->nky;
    pme->nkz         = ir->nkz;
    pme->pme_order   = ir->pme_order;
    pme->epsilon_r   = ir->epsilon_r;
    
    /* Currently pme.c supports only the fft5d FFT code.
     * Therefore the grid always needs to be divisible by nnodes.
     * When the old 1D code is also supported again, change this check.
     *
     * This check should be done before calling gmx_pme_init
     * and fplog should be passed iso stderr.
     *
    if (pme->ndecompdim >= 2)
    */
    if (pme->ndecompdim >= 1)
    {
        /*
        gmx_pme_check_grid_restrictions(pme->nodeid==0 ? stderr : NULL,
                                        'x',nnodes_major,&pme->nkx);
        gmx_pme_check_grid_restrictions(pme->nodeid==0 ? stderr : NULL,
                                        'y',nnodes_minor,&pme->nky);
        */
    }

    if (pme->nkx <= pme->pme_order*(pme->nnodes_major > 1 ? 2 : 1) ||
        pme->nky <= pme->pme_order*(pme->nnodes_minor > 1 ? 2 : 1) ||
        pme->nkz <= pme->pme_order)
    {
        gmx_fatal(FARGS,"The pme grid dimensions need to be larger than pme_order (%d) and in parallel larger than 2*pme_order for x and/or y",pme->pme_order);
    }

    if (pme->nnodes > 1) {
        double imbal;

#ifdef GMX_MPI
        MPI_Type_contiguous(DIM, mpi_type, &(pme->rvec_mpi));
        MPI_Type_commit(&(pme->rvec_mpi));
#endif
        
        /* Note that the charge spreading and force gathering, which usually
         * takes about the same amount of time as FFT+solve_pme,
         * is always fully load balanced
         * (unless the charge distribution is inhomogeneous).
         */
        
        imbal = pme_load_imbalance(pme);
        if (imbal >= 1.2 && pme->nodeid_major == 0 && pme->nodeid_minor == 0)
        {
            fprintf(stderr,
                    "\n"
                    "NOTE: The load imbalance in PME FFT and solve is %d%%.\n"
                    "      For optimal PME load balancing\n"
                    "      PME grid_x (%d) and grid_y (%d) should be divisible by #PME_nodes_x (%d)\n"
                    "      and PME grid_y (%d) and grid_z (%d) should be divisible by #PME_nodes_y (%d)\n"
                    "\n",
                    (int)((imbal-1)*100 + 0.5),
                    pme->nkx,pme->nky,pme->nnodes_major,
                    pme->nky,pme->nkz,pme->nnodes_minor);
        }
    }

    init_overlap_comm(&pme->overlap[0],pme->pme_order,
#ifdef GMX_MPI
                      pme->mpi_comm_d[0],
#endif
                      pme->nnodes_major,pme->nodeid_major,pme->nkx);
    
    init_overlap_comm(&pme->overlap[1],pme->pme_order,
#ifdef GMX_MPI
                      pme->mpi_comm_d[1],
#endif
                      pme->nnodes_minor,pme->nodeid_minor,pme->nky);
    
    snew(pme->bsp_mod[XX],pme->nkx);
    snew(pme->bsp_mod[YY],pme->nky);
    snew(pme->bsp_mod[ZZ],pme->nkz);
    
    /* Allocate data for the interpolation grid, including overlap */
    pme->pmegrid_nx = pme->overlap[0].s2g1[pme->nodeid_major] -
                      pme->overlap[0].s2g0[pme->nodeid_major];
    pme->pmegrid_ny = pme->overlap[1].s2g1[pme->nodeid_minor] - 
                      pme->overlap[1].s2g0[pme->nodeid_minor];
    pme->pmegrid_nz = pme->nkz + pme->pme_order - 1;
    
    pme->pmegrid_start_ix = pme->overlap[0].s2g0[pme->nodeid_major];
    pme->pmegrid_start_iy = pme->overlap[1].s2g0[pme->nodeid_minor];
    pme->pmegrid_start_iz = 0;
    
    make_gridindex5_to_localindex(pme->nkx,
                                  pme->pmegrid_start_ix,
                                  pme->pmegrid_nx - (pme->pme_order-1),
                                  &pme->nnx,&pme->fshx);
    make_gridindex5_to_localindex(pme->nky,
                                  pme->pmegrid_start_iy,
                                  pme->pmegrid_ny - (pme->pme_order-1),
                                  &pme->nny,&pme->fshy);
    make_gridindex5_to_localindex(pme->nkz,
                                  pme->pmegrid_start_iz,
                                  pme->pmegrid_nz - (pme->pme_order-1),
                                  &pme->nnz,&pme->fshz);
    
    snew(pme->pmegridA,pme->pmegrid_nx*pme->pmegrid_ny*pme->pmegrid_nz);
    
    /* For non-divisible grid we need pme_order iso pme_order-1 */
    /* x overlap is copied in place: take padding into account.
     * y is always copied through a buffer: we don't need padding in z,
     * but we do need the overlap in x because of the communication order.
     */
    bufsizex = pme->pme_order*pme->pmegrid_ny*pme->pmegrid_nz;
    bufsizey = pme->pme_order*pme->pmegrid_nx*pme->nkz;
    bufsize  = (bufsizex>bufsizey) ? bufsizex : bufsizey;
    
    snew(pme->pmegrid_sendbuf,bufsize);
    snew(pme->pmegrid_recvbuf,bufsize);
    
    ndata[0] = pme->nkx;
    ndata[1] = pme->nky;
    ndata[2] = pme->nkz;
    
    /* This routine will allocate the grid data to fit the FFTs */
    gmx_parallel_3dfft_init(&pme->pfft_setupA,ndata,
                            &pme->fftgridA,&pme->cfftgridA,
                            pme->mpi_comm_d,
                            pme->overlap[0].s2g0,pme->overlap[1].s2g0,
                            bReproducible);
    
    if (bFreeEnergy)
    {
        snew(pme->pmegridB,pme->pmegrid_nx*pme->pmegrid_ny*pme->pmegrid_nz);    
        gmx_parallel_3dfft_init(&pme->pfft_setupB,ndata,
                                &pme->fftgridB,&pme->cfftgridB,
                                pme->mpi_comm_d,
                                pme->overlap[0].s2g0,pme->overlap[1].s2g0,
                                bReproducible);
    } else 
    {
        pme->pmegridB    = NULL;
        pme->fftgridB    = NULL;
        pme->cfftgridB   = NULL;
    }
    
    make_bspline_moduli(pme->bsp_mod,pme->nkx,pme->nky,pme->nkz,pme->pme_order);
    
    /* Use atc[0] for spreading */
    init_atomcomm(pme,&pme->atc[0],cr,nnodes_major > 1 ? 0 : 1,TRUE);
    if (pme->ndecompdim >= 2)
    {
        init_atomcomm(pme,&pme->atc[1],cr,1,FALSE);
    }
    
    if (pme->nnodes == 1) {
        pme->atc[0].n = homenr;
        pme_realloc_atomcomm_things(&pme->atc[0]);
    }
    
    /* Use fft5d, order after FFT is y major, z, x minor */
    pme->work_nalloc = pme->nkx;
    snew(pme->work_mhx,pme->work_nalloc);
    snew(pme->work_mhy,pme->work_nalloc);
    snew(pme->work_mhz,pme->work_nalloc);
    snew(pme->work_m2,pme->work_nalloc);
    snew(pme->work_denom,pme->work_nalloc);
    /* Allocate an aligned pointer for SSE operations, including 3 extra
     * elements at the end since SSE operates on 4 elements at a time.
     */
    snew(pme->work_tmp1_alloc,pme->work_nalloc+8);
    pme->work_tmp1 = (real *) (((size_t) pme->work_tmp1_alloc + 16) & (~((size_t) 15)));
    snew(pme->work_m2inv,pme->work_nalloc);

    *pmedata = pme;
    
    return 0;
}

static void spread_on_grid(gmx_pme_t pme,
                           pme_atomcomm_t *atc,real *grid,
                           gmx_bool bCalcSplines,gmx_bool bSpread)
{    
    if (bCalcSplines)
    {
    
        /* Compute fftgrid index for all atoms,
         * with help of some extra variables.
         */
        calc_interpolation_idx(pme,atc);
        
        /* make local bsplines  */
        make_bsplines(atc->theta,atc->dtheta,pme->pme_order,
                      atc->fractx,atc->n,atc->q,pme->bFEP);
    }    
    
    if (bSpread)
    {
        /* put local atoms on grid. */
        spread_q_bsplines(pme,atc,grid);
    }
}

void gmx_pme_calc_energy(gmx_pme_t pme,int n,rvec *x,real *q,real *V)
{
    pme_atomcomm_t *atc;
    real *grid;

    if (pme->nnodes > 1)
    {
        gmx_incons("gmx_pme_calc_energy called in parallel");
    }
    if (pme->bFEP > 1)
    {
        gmx_incons("gmx_pme_calc_energy with free energy");
    }

    atc = &pme->atc_energy;
    atc->nslab     = 1;
    atc->bSpread   = TRUE;
    atc->pme_order = pme->pme_order;
    atc->n         = n;
    pme_realloc_atomcomm_things(atc);
    atc->x         = x;
    atc->q         = q;
    
    /* We only use the A-charges grid */
    grid = pme->pmegridA;

    spread_on_grid(pme,atc,NULL,TRUE,FALSE);

    *V = gather_energy_bsplines(pme,grid,atc);
}


static void reset_pmeonly_counters(t_commrec *cr,gmx_wallcycle_t wcycle,
        t_nrnb *nrnb,t_inputrec *ir, gmx_large_int_t step_rel)
{
    /* Reset all the counters related to performance over the run */
    wallcycle_stop(wcycle,ewcRUN);
    wallcycle_reset_all(wcycle);
    init_nrnb(nrnb);
    ir->init_step += step_rel;
    ir->nsteps    -= step_rel;
    wallcycle_start(wcycle,ewcRUN);
}


int gmx_pmeonly(gmx_pme_t pme,
                t_commrec *cr,    t_nrnb *nrnb,
                gmx_wallcycle_t wcycle,
                real ewaldcoeff,  gmx_bool bGatherOnly,
                t_inputrec *ir)
{
    gmx_pme_pp_t pme_pp;
    int  natoms;
    matrix box;
    rvec *x_pp=NULL,*f_pp=NULL;
    real *chargeA=NULL,*chargeB=NULL;
    real lambda=0;
    int  maxshift_x=0,maxshift_y=0;
    real energy,dvdlambda;
    matrix vir;
    float cycles;
    int  count;
    gmx_bool bEnerVir;
    gmx_large_int_t step,step_rel;
    
    
    pme_pp = gmx_pme_pp_init(cr);
    
    init_nrnb(nrnb);
    
    count = 0;
    do /****** this is a quasi-loop over time steps! */
    {
        /* Domain decomposition */
        natoms = gmx_pme_recv_q_x(pme_pp,
                                  &chargeA,&chargeB,box,&x_pp,&f_pp,
                                  &maxshift_x,&maxshift_y,
                                  &pme->bFEP,&lambda,
                                  &bEnerVir,
                                  &step);
        
        if (natoms == -1) {
            /* We should stop: break out of the loop */
            break;
        }
        
        step_rel = step - ir->init_step;
        
        if (count == 0)
            wallcycle_start(wcycle,ewcRUN);
        
        wallcycle_start(wcycle,ewcPMEMESH);
        
        dvdlambda = 0;
        clear_mat(vir);
        gmx_pme_do(pme,0,natoms,x_pp,f_pp,chargeA,chargeB,box,
                   cr,maxshift_x,maxshift_y,nrnb,wcycle,vir,ewaldcoeff,
                   &energy,lambda,&dvdlambda,
                   GMX_PME_DO_ALL_F | (bEnerVir ? GMX_PME_CALC_ENER_VIR : 0));
        
        cycles = wallcycle_stop(wcycle,ewcPMEMESH);
        
        gmx_pme_send_force_vir_ener(pme_pp,
                                    f_pp,vir,energy,dvdlambda,
                                    cycles);
        
        count++;

        if (step_rel == wcycle_get_reset_counters(wcycle))
        {
            /* Reset all the counters related to performance over the run */
            reset_pmeonly_counters(cr,wcycle,nrnb,ir,step_rel);
            wcycle_set_reset_counters(wcycle, 0);
        }
        
    } /***** end of quasi-loop, we stop with the break above */
    while (TRUE);
    
    return 0;
}

int gmx_pme_do(gmx_pme_t pme,
               int start,       int homenr,
               rvec x[],        rvec f[],
               real *chargeA,   real *chargeB,
               matrix box,	t_commrec *cr,
               int  maxshift_x, int maxshift_y,
               t_nrnb *nrnb,    gmx_wallcycle_t wcycle,
               matrix vir,      real ewaldcoeff,
               real *energy,    real lambda, 
               real *dvdlambda, int flags)
{
    int     q,d,i,j,ntot,npme;
    int     nx,ny,nz;
    int     n_d,local_ny;
    int     loop_count;
    pme_atomcomm_t *atc=NULL;
    real *  grid=NULL;
    real    *ptr;
    rvec    *x_d,*f_d;
    real    *charge=NULL,*q_d,vol;
    real    energy_AB[2];
    matrix  vir_AB[2];
    gmx_bool    bClearF;
    gmx_parallel_3dfft_t pfft_setup;
    real *  fftgrid;
    t_complex * cfftgrid;

    if (pme->nnodes > 1) {
        atc = &pme->atc[0];
        atc->npd = homenr;
        if (atc->npd > atc->pd_nalloc) {
            atc->pd_nalloc = over_alloc_dd(atc->npd);
            srenew(atc->pd,atc->pd_nalloc);
        }
        atc->maxshift = (atc->dimind==0 ? maxshift_x : maxshift_y);
    }
    else
    {
        /* This could be necessary for TPI */
        pme->atc[0].n = homenr;
    }
    
    for(q=0; q<(pme->bFEP ? 2 : 1); q++) {
        if (q == 0) {
            grid = pme->pmegridA;
            fftgrid = pme->fftgridA;
            cfftgrid = pme->cfftgridA;
            pfft_setup = pme->pfft_setupA;
            charge = chargeA+start;
        } else {
            grid = pme->pmegridB;
            fftgrid = pme->fftgridB;
            cfftgrid = pme->cfftgridB;
            pfft_setup = pme->pfft_setupB;
            charge = chargeB+start;
        }
        /* Unpack structure */
        if (debug) {
            fprintf(debug,"PME: nnodes = %d, nodeid = %d\n",
                    cr->nnodes,cr->nodeid);
            fprintf(debug,"Grid = %p\n",(void*)grid);
            if (grid == NULL)
                gmx_fatal(FARGS,"No grid!");
        }
        where();
        
        m_inv_ur0(box,pme->recipbox); 

        if (pme->nnodes == 1) {
            atc = &pme->atc[0];
            if (DOMAINDECOMP(cr)) {
                atc->n = homenr;
                pme_realloc_atomcomm_things(atc);
            }
            atc->x = x;
            atc->q = charge;
            atc->f = f;
        } else {
            wallcycle_start(wcycle,ewcPME_REDISTXF);
            for(d=pme->ndecompdim-1; d>=0; d--)
            {
                if (d == pme->ndecompdim-1)
                {
                    n_d = homenr;
                    x_d = x + start;
                    q_d = charge;
                }
                else
                {
                    n_d = pme->atc[d+1].n;
                    x_d = atc->x;
                    q_d = atc->q;
                }
                atc = &pme->atc[d];
                atc->npd = n_d;
                if (atc->npd > atc->pd_nalloc) {
                    atc->pd_nalloc = over_alloc_dd(atc->npd);
                    srenew(atc->pd,atc->pd_nalloc);
                }
                atc->maxshift = (atc->dimind==0 ? maxshift_x : maxshift_y);
                pme_calc_pidx(n_d,pme->recipbox,x_d,atc);
                where();
                
                GMX_BARRIER(cr->mpi_comm_mygroup);
                /* Redistribute x (only once) and qA or qB */
                if (DOMAINDECOMP(cr)) {
                    dd_pmeredist_x_q(pme, n_d, q==0, x_d, q_d, atc);
                } else {
                    pmeredist_pd(pme, TRUE, n_d, q==0, x_d, q_d, atc);
                }
            }
            where();

            wallcycle_stop(wcycle,ewcPME_REDISTXF);
        }
        
        if (debug)
            fprintf(debug,"Node= %6d, pme local particles=%6d\n",
                    cr->nodeid,atc->n);

        if (flags & GMX_PME_SPREAD_Q)
        {
            wallcycle_start(wcycle,ewcPME_SPREADGATHER);

            /* Spread the charges on a grid */
            GMX_MPE_LOG(ev_spread_on_grid_start);
            
            /* Spread the charges on a grid */
            spread_on_grid(pme,&pme->atc[0],grid,q==0,TRUE);
            GMX_MPE_LOG(ev_spread_on_grid_finish);

            if (q == 0)
            {
                inc_nrnb(nrnb,eNR_WEIGHTS,DIM*atc->n);
            }
            inc_nrnb(nrnb,eNR_SPREADQBSP,
                     pme->pme_order*pme->pme_order*pme->pme_order*atc->n);

            wrap_periodic_pmegrid(pme,grid);

            /* sum contributions to local grid from other nodes */
#ifdef GMX_MPI
            if (pme->nnodes > 1) {
                GMX_BARRIER(cr->mpi_comm_mygroup);
                gmx_sum_qgrid_dd(pme,grid,GMX_SUM_QGRID_FORWARD);
                where();
            }
#endif
            where();

            copy_pmegrid_to_fftgrid(pme,grid,fftgrid);

            wallcycle_stop(wcycle,ewcPME_SPREADGATHER);
        }
         
        if (flags & GMX_PME_SOLVE)
        {
            /* do 3d-fft */ 
            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_gmxfft3d_start);
            wallcycle_start(wcycle,ewcPME_FFT);
            gmx_parallel_3dfft_execute(pfft_setup,GMX_FFT_REAL_TO_COMPLEX,fftgrid,cfftgrid);
            wallcycle_stop(wcycle,ewcPME_FFT);
            GMX_MPE_LOG(ev_gmxfft3d_finish);
            where();
            
            /* solve in k-space for our local cells */
            vol = det(box);
            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_solve_pme_start);
            wallcycle_start(wcycle,ewcPME_SOLVE);
            loop_count =
                solve_pme_yzx(pme,cfftgrid,ewaldcoeff,vol,
                              flags & GMX_PME_CALC_ENER_VIR,
                              &energy_AB[q],vir_AB[q]);
            wallcycle_stop(wcycle,ewcPME_SOLVE);
            where();
            GMX_MPE_LOG(ev_solve_pme_finish);
            inc_nrnb(nrnb,eNR_SOLVEPME,loop_count);
        }

        if ((flags & GMX_PME_CALC_F) ||
            (flags & GMX_PME_CALC_POT))
        {
            
            /* do 3d-invfft */
            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_gmxfft3d_start);
            where();
            wallcycle_start(wcycle,ewcPME_FFT);
            gmx_parallel_3dfft_execute(pfft_setup,GMX_FFT_COMPLEX_TO_REAL,cfftgrid,fftgrid);
            wallcycle_stop(wcycle,ewcPME_FFT);

            where();
            GMX_MPE_LOG(ev_gmxfft3d_finish);

            if (pme->nodeid == 0)
            {
                ntot = pme->nkx*pme->nky*pme->nkz;
                npme  = ntot*log((real)ntot)/log(2.0);
                inc_nrnb(nrnb,eNR_FFT,2*npme);
            }

            wallcycle_start(wcycle,ewcPME_SPREADGATHER);

            copy_fftgrid_to_pmegrid(pme,fftgrid,grid);

            /* distribute local grid to all nodes */
#ifdef GMX_MPI
            if (pme->nnodes > 1) {
                GMX_BARRIER(cr->mpi_comm_mygroup);
                gmx_sum_qgrid_dd(pme,grid,GMX_SUM_QGRID_BACKWARD);
            }
#endif
            where();

            unwrap_periodic_pmegrid(pme,grid);
        }

        if (flags & GMX_PME_CALC_F)
        {
            /* interpolate forces for our local atoms */
            GMX_BARRIER(cr->mpi_comm_mygroup);
            GMX_MPE_LOG(ev_gather_f_bsplines_start);

            where();
            
            /* If we are running without parallelization,
             * atc->f is the actual force array, not a buffer,
             * therefore we should not clear it.
             */
            bClearF = (q == 0 && PAR(cr));
            gather_f_bsplines(pme,grid,bClearF,&pme->atc[0],
                              pme->bFEP ? (q==0 ? 1.0-lambda : lambda) : 1.0);
            where();
            
            GMX_MPE_LOG(ev_gather_f_bsplines_finish);
            
            inc_nrnb(nrnb,eNR_GATHERFBSP,
                     pme->pme_order*pme->pme_order*pme->pme_order*pme->atc[0].n);
            wallcycle_stop(wcycle,ewcPME_SPREADGATHER);
       }
    } /* of q-loop */
    
    if ((flags & GMX_PME_CALC_F) && pme->nnodes > 1) {
        wallcycle_start(wcycle,ewcPME_REDISTXF);
        for(d=0; d<pme->ndecompdim; d++)
        {
            atc = &pme->atc[d];
            if (d == pme->ndecompdim - 1)
            {
                n_d = homenr;
                f_d = f + start;
            }
            else
            {
                n_d = pme->atc[d+1].n;
                f_d = pme->atc[d+1].f;
            }
            GMX_BARRIER(cr->mpi_comm_mygroup);
            if (DOMAINDECOMP(cr)) {
                dd_pmeredist_f(pme,atc,n_d,f_d,
                               d==pme->ndecompdim-1 && pme->bPPnode);
            } else {
                pmeredist_pd(pme, FALSE, n_d, TRUE, f_d, NULL, atc);
            }
        }

        wallcycle_stop(wcycle,ewcPME_REDISTXF);
    }
    where();
    
    if (!pme->bFEP) {
        *energy = energy_AB[0];
        m_add(vir,vir_AB[0],vir);
    } else {
        *energy = (1.0-lambda)*energy_AB[0] + lambda*energy_AB[1];
        *dvdlambda += energy_AB[1] - energy_AB[0];
        for(i=0; i<DIM; i++)
            for(j=0; j<DIM; j++)
                vir[i][j] += (1.0-lambda)*vir_AB[0][i][j] + lambda*vir_AB[1][i][j];
    }

    if (debug)
        fprintf(debug,"PME mesh energy: %g\n",*energy);
    
    return 0;
}
