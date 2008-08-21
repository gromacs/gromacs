/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * $Id$
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
#include "fftgrid.h"
#include "gmx_fatal.h"
#include "pme.h"
#include "network.h"
#include "physics.h"
#include "nrnb.h"
#include "copyrite.h"
#include "gmx_wallcycle.h"
#ifdef GMX_MPI
#include <mpi.h>
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

#ifdef GMX_MPI
MPI_Datatype  rvec_mpi;
#endif

/* Internal datastructures */
typedef struct {
    int snd0;
    int snds;
    int rcv0;
    int rcvs;
} pme_grid_comm_t;

typedef struct {
#ifdef GMX_MPI
    MPI_Comm mpi_comm;
#endif
    int  nslab;
    int  *s2g;
    int  nleftbnd,nrightbnd;  /* The number of nodes to communicate with */
    int  nodeid,*leftid,*rightid;
    pme_grid_comm_t *leftc,*rightc;
} pme_overlap_t;

typedef struct {
#ifdef GMX_MPI
    MPI_Comm mpi_comm;
#endif
    int  nslab;
    int  nodeid;

    int  *node_dest;        /* The nodes to send x and q to with DD */
    int  *node_src;         /* The nodes to receive x and q from with DD */
    int  *buf_index;        /* Index for commnode into the buffers */

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
    bool bSpread;           /* These coordinates are used for spreading */
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
    int  nnodes;             /* The number of nodes doing PME */
#ifdef GMX_MPI
    MPI_Comm mpi_comm;
#endif

    bool bPPnode;            /* Node also does particle-particle forces */
    bool bFEP;               /* Compute Free energy contribution */
    int nkx,nky,nkz;         /* Grid dimensions */
    int pme_order;
    real epsilon_r;           
    t_fftgrid *gridA,*gridB;
    int  *nnx,*nny,*nnz;
    
    pme_atomcomm_t atc[2];
    matrix    recipbox;
    splinevec bsp_mod;
    
    pme_overlap_t overlap[2];
    
    rvec *bufv;             /* Communication buffer */
    real *bufr;             /* Communication buffer */
    int  buf_nalloc;        /* The communication buffer size */
} t_gmx_pme;

/* The following stuff is needed for signal handling on the PME nodes. 
 * signal_handler needs to be defined in md.c, the bGot..Signal variables
 * here */ 
extern RETSIGTYPE signal_handler(int n);

volatile bool bGotTermSignal = FALSE, bGotUsr1Signal = FALSE; 

/* #define SORTPME */

static void pr_grid_dist(FILE *fp,char *title,t_fftgrid *grid)
{
  int     i,j,k,l,ntoti,ntot=0;
  int     nx,ny,nz,nx2,ny2,nz2,la12,la2;
  real *  ptr;

  /* Unpack structure */
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  for(i=0; (i<nx); i++) {
    ntoti=0;
    for(j=0; (j<ny); j++)
      for(k=0; (k<nz); k++) {
        l= INDEX(i,j,k);
        if (ptr[l] != 0) {
          ntoti++;
          ntot++;
        }
      }
          fprintf(fp,"%-12s  %5d  %5d\n",title,i,ntoti);
  }
  fprintf(fp,"%d non zero elements in %s\n",ntot,title);
}
/* test */

static void calc_recipbox(matrix box,matrix recipbox)
{
  /* Save some time by assuming upper right part is zero */

  real tmp=1.0/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);

  recipbox[XX][XX]=box[YY][YY]*box[ZZ][ZZ]*tmp;
  recipbox[XX][YY]=0;
  recipbox[XX][ZZ]=0;
  recipbox[YY][XX]=-box[YY][XX]*box[ZZ][ZZ]*tmp;
  recipbox[YY][YY]=box[XX][XX]*box[ZZ][ZZ]*tmp;
  recipbox[YY][ZZ]=0;
  recipbox[ZZ][XX]=(box[YY][XX]*box[ZZ][YY]-box[YY][YY]*box[ZZ][XX])*tmp;
  recipbox[ZZ][YY]=-box[ZZ][YY]*box[XX][XX]*tmp;
  recipbox[ZZ][ZZ]=box[XX][XX]*box[YY][YY]*tmp;
}

static void calc_idx(gmx_pme_t pme,rvec x[])
{
    pme_atomcomm_t *atc;
    int  i;
    int  *idxptr,tix,tiy,tiz;
    real *xptr,*fptr,tx,ty,tz;
    real rxx,ryx,ryy,rzx,rzy,rzz;
    int  nx,ny,nz,nx2,ny2,nz2,la12,la2;
    real *ptr;

#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined GMX_DOUBLE && defined GMX_X86TRUNC)
    int x86_cw,x86_cwsave;

    asm("fnstcw %0" : "=m" (*&x86_cwsave));
    x86_cw = x86_cwsave | 3072;
    asm("fldcw %0" : : "m" (*&x86_cw));
#define x86trunc(a,b) asm("fld %1\nfistpl %0\n" : "=m" (*&b) : "f" (a));
#endif
 
    atc = &pme->atc[0];

    /* Unpack structure */
    unpack_fftgrid(pme->gridA,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
    
    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];
    
    for(i=0; (i<atc->n); i++) {
        xptr   = x[i];
        idxptr = atc->idx[i];
        fptr   = atc->fractx[i];
        
        /* Fractional coordinates along box vectors */
        tx = nx2 + nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx );
        ty = ny2 + ny * (                  xptr[YY] * ryy + xptr[ZZ] * rzy );
        tz = nz2 + nz * (                                   xptr[ZZ] * rzz );
        
#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined GMX_DOUBLE && defined GMX_X86TRUNC)
        x86trunc(tx,tix);
        x86trunc(ty,tiy);
        x86trunc(tz,tiz);
#else
        tix = (int)(tx);
        tiy = (int)(ty);
        tiz = (int)(tz);
#endif
        
        fptr[XX] = tx - tix;
        fptr[YY] = ty - tiy;
        fptr[ZZ] = tz - tiz;   
        
        idxptr[XX] = pme->nnx[tix];
        idxptr[YY] = pme->nny[tiy];
        idxptr[ZZ] = pme->nnz[tiz];

#ifdef DEBUG
        range_check(idxptr[XX],0,nx);
        range_check(idxptr[YY],0,ny);
        range_check(idxptr[ZZ],0,nz);
#endif
  }  
#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined GMX_DOUBLE && defined GMX_X86TRUNC)  
    asm("fldcw %0" : : "m" (*&x86_cwsave));
#endif
}

static void pme_calc_pidx(int npmenodes,
                          int natoms,matrix box, rvec x[],
                          pme_atomcomm_t *atc)
{
  int  i;
  int  tix;
  real *xptr,tx;
  real rxx,ryx,rzx;
  matrix recipbox;
    
  /* Calculate PME task index (pidx) for each grid index.
   * Here we always assign equally sized slabs to each node
   * for load balancing reasons (the PME grid spacing is not used).
   */

  /* Reset the count */
  for(i=0; i<npmenodes; i++)
    atc->count[i] = 0;

  calc_recipbox(box,recipbox);
  rxx = recipbox[XX][XX];
  ryx = recipbox[YY][XX];
  rzx = recipbox[ZZ][XX];
  /* Calculate the node index in x-dimension */
  for(i=0; (i<natoms); i++) {
    xptr   = x[i];
    /* Fractional coordinates along box vectors */
    tx = npmenodes * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx );
    tix = (int)(tx + npmenodes) - npmenodes;
    if (tix < 0)
      tix += npmenodes;
    else if (tix >= npmenodes)
      tix -= npmenodes;
    atc->pd[i] = tix;
    atc->count[tix]++;
  }
}

static void pme_realloc_atomcomm_things(pme_atomcomm_t *atc)
{
    int i;
    
    if (atc->n > atc->nalloc) {
        atc->nalloc = over_alloc_dd(atc->n);
        
        if (atc->nslab > 1) {
            srenew(atc->x,atc->nalloc);
            srenew(atc->q,atc->nalloc);
            srenew(atc->f,atc->nalloc);
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

static void pmeredist(gmx_pme_t pme, bool forw,
                      int n, bool bXF, rvec *x_f, real *charge,
                      pme_atomcomm_t *atc)
/* Redistribute particle data for PME calculation */
/* domain decomposition by x coordinate           */
{
    static int *scounts, *rcounts,*sdispls, *rdispls, *sidx;
    static real *buf=NULL;
    static int buf_nalloc=0;
    int *idxa;
    int i, ii;
    static bool bFirst=TRUE;
    
    if(bFirst) {
        snew(scounts,atc->nslab);
        snew(rcounts,atc->nslab);
        snew(sdispls,atc->nslab);
        snew(rdispls,atc->nslab);
        snew(sidx,atc->nslab);
        bFirst=FALSE;
    }
    if (n > buf_nalloc) {
        buf_nalloc = over_alloc_dd(n);
        srenew(buf,buf_nalloc*DIM);
    }
    
    idxa = atc->pd;

#ifdef GMX_MPI
    if (forw && bXF) {
        /* forward, redistribution from pp to pme */ 
        
        /* Calculate send counts and exchange them with other nodes */
        for(i=0; (i<atc->nslab); i++) scounts[i]=0;
        for(i=0; (i<n); i++) scounts[idxa[i]]++;
        MPI_Alltoall( scounts, 1, MPI_INT, rcounts, 1, MPI_INT, atc->mpi_comm);
        
        /* Calculate send and receive displacements and index into send buffer */
        sdispls[0]=0;
        rdispls[0]=0;
        sidx[0]=0;
        for(i=1; i<atc->nslab; i++) {
            sdispls[i]=sdispls[i-1]+scounts[i-1];
            rdispls[i]=rdispls[i-1]+rcounts[i-1];
            sidx[i]=sdispls[i];
        }
        /* Total # of particles to be received */
        atc->n = rdispls[atc->nslab-1] + rcounts[atc->nslab-1];
        
        pme_realloc_atomcomm_things(atc);
        
        /* Copy particle coordinates into send buffer and exchange*/
        for(i=0; (i<n); i++) {
            ii=DIM*sidx[idxa[i]];
            sidx[idxa[i]]++;
            buf[ii+XX]=x_f[i][XX];
            buf[ii+YY]=x_f[i][YY];
            buf[ii+ZZ]=x_f[i][ZZ];
        }
        MPI_Alltoallv(buf, scounts, sdispls, rvec_mpi,
                      atc->x, rcounts, rdispls, rvec_mpi,
                      atc->mpi_comm);
    }
    if (forw) {
        /* Copy charge into send buffer and exchange*/
        for(i=0; i<atc->nslab; i++) sidx[i]=sdispls[i];
        for(i=0; (i<n); i++) {
            ii=sidx[idxa[i]];
            sidx[idxa[i]]++;
            buf[ii]=charge[i];
        }
        MPI_Alltoallv(buf, scounts, sdispls, mpi_type,
                      atc->q, rcounts, rdispls, mpi_type,
                      atc->mpi_comm);
    }
    else { /* backward, redistribution from pme to pp */ 
        MPI_Alltoallv(atc->f, rcounts, rdispls, rvec_mpi,
                      buf, scounts, sdispls, rvec_mpi,
                      atc->mpi_comm);
        
        /* Copy data from receive buffer */
        for(i=0; i<atc->nslab; i++)
            sidx[i] = sdispls[i];
        for(i=0; (i<n); i++) {
            ii = DIM*sidx[idxa[i]];
            x_f[i][XX] += buf[ii+XX];
            x_f[i][YY] += buf[ii+YY];
            x_f[i][ZZ] += buf[ii+ZZ];
            sidx[idxa[i]]++;
        }
    }
#endif 
}

static void pme_dd_sendrecv(pme_atomcomm_t *atc,
                            bool bBackward,int shift,
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

static void dd_pmeredist_x_q(gmx_pme_t pme, int maxshift,
                             int n, bool bX, rvec *x, real *charge,
                             pme_atomcomm_t *atc)
{
    int *commnode,*buf_index;
    int nnodes_comm,i,nsend,local_pos,buf_pos,node,scount,rcount;
    
    commnode  = atc->node_dest;
    buf_index = atc->buf_index;
    
    nnodes_comm = min(2*maxshift,atc->nslab-1);
    
    nsend = 0;
    for(i=0; i<nnodes_comm; i++) {
        buf_index[commnode[i]] = nsend;
        nsend += atc->count[commnode[i]];
    }
    if (bX) {
        if (atc->count[pme->nodeid] + nsend != n)
            gmx_fatal(FARGS,"%d particles communicated to PME node %d are more than a cell length out of the domain decomposition cell of their charge group",
                      n - (atc->count[pme->nodeid] + nsend),pme->nodeid);
        
        if (nsend > pme->buf_nalloc) {
            pme->buf_nalloc = over_alloc_dd(nsend);
            srenew(pme->bufv,pme->buf_nalloc);
            srenew(pme->bufr,pme->buf_nalloc);
        }
        
        atc->n = atc->count[pme->nodeid];
        for(i=0; i<nnodes_comm; i++) {
            scount = atc->count[commnode[i]];
            /* Communicate the count */
            if (debug)
                fprintf(debug,"PME node %d send to node %d: %d\n",
                        pme->nodeid,commnode[i],scount);
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
        if (node == pme->nodeid) {
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
                           int maxshift,
                           int n, rvec *f,
                           bool bAddF)
{
  int *commnode,*buf_index;
  int nnodes_comm,local_pos,buf_pos,i,scount,rcount,node;

  commnode  = atc->node_dest;
  buf_index = atc->buf_index;

  nnodes_comm = min(2*maxshift,atc->nslab-1);

  local_pos = atc->count[pme->nodeid];
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
            if (node == pme->nodeid)
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
            if (node == pme->nodeid)
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

static void gmx_sum_qgrid_dd(pme_overlap_t *ol,t_fftgrid *grid,int direction)
{
    static real *tmp=NULL;
    int b,i;
    int la12r,localsize;
    pme_grid_comm_t *pgc;
    real *from, *to;
#ifdef GMX_MPI
    MPI_Status stat;
#endif
    
    GMX_MPE_LOG(ev_sum_qgrid_start);
    
#ifdef GMX_MPI
    
    la12r      = grid->la12r;
    localsize  = la12r*grid->pfft.local_nx;
    
    if (grid->workspace) {
        tmp = grid->workspace;
    } else {
        if (tmp == NULL) {
            snew(tmp,localsize);
        }
    }
    
    if (direction == GMX_SUM_QGRID_FORWARD) { 
        /* sum contributions to local grid */
        /* Send left boundaries */
        for(b=0; b<ol->nleftbnd; b++) {
            pgc = &ol->leftc[b];
            from = grid->ptr + la12r*pgc->snd0;
            to   = grid->ptr + la12r*pgc->rcv0;
            MPI_Sendrecv(from,la12r*pgc->snds,mpi_type,
                         ol->leftid[b], ol->nodeid,
                         tmp, la12r*pgc->rcvs,mpi_type,
                         ol->rightid[b],ol->rightid[b],
                         ol->mpi_comm,&stat);
            GMX_MPE_LOG(ev_test_start); 
            for(i=0; (i<la12r*pgc->rcvs); i++) {
                to[i] += tmp[i];
            }
        }
        GMX_MPE_LOG(ev_test_finish);
        /* Send right boundaries */
        for(b=0; b<ol->nrightbnd; b++) {
            pgc = &ol->rightc[b];
            from = grid->ptr + la12r*pgc->snd0;
            to   = grid->ptr + la12r*pgc->rcv0;
            MPI_Sendrecv(from,la12r*pgc->snds,mpi_type,
                         ol->rightid[b],ol->nodeid,
                         tmp, la12r*pgc->rcvs,mpi_type,
                         ol->leftid[b], ol->leftid[b],
                         ol->mpi_comm,&stat);
            GMX_MPE_LOG(ev_test_start); 
            for(i=0; (i<la12r*pgc->rcvs); i++) {
                to[i] += tmp[i];
            }
        }
        GMX_MPE_LOG(ev_test_finish);
    }
    else if (direction  == GMX_SUM_QGRID_BACKWARD) { 
        /* distribute local grid to all processors */
        /* Send right boundaries */
        for(b=0; b<ol->nrightbnd; b++) {
            pgc = &ol->rightc[b];
            from = grid->ptr + la12r*pgc->rcv0;
            to   = grid->ptr + la12r*pgc->snd0;
            MPI_Sendrecv(from,la12r*pgc->rcvs,mpi_type,
                         ol->leftid[b], ol->nodeid,
                         to,  la12r*pgc->snds,mpi_type,
                         ol->rightid[b],ol->rightid[b],
                         ol->mpi_comm,&stat);
        }
        /* Send left boundaries */
        for(b=0; b<ol->nleftbnd; b++) {
            pgc = &ol->leftc[b];
            from = grid->ptr + la12r*pgc->rcv0;
            to   = grid->ptr + la12r*pgc->snd0;
            MPI_Sendrecv(from,la12r*pgc->rcvs,mpi_type,
                         ol->rightid[b],ol->nodeid,
                         to,  la12r*pgc->snds,mpi_type,
                         ol->leftid[b], ol->leftid[b],
                         ol->mpi_comm,&stat);
        }
    }
    else {
        gmx_fatal(FARGS,"Invalid direction %d for summing qgrid",direction);
    }
    
#else
    gmx_fatal(FARGS,"Parallel grid summation requires MPI and FFTW.\n");    
#endif
    GMX_MPE_LOG(ev_sum_qgrid_finish);
}

void gmx_sum_qgrid(gmx_pme_t gmx,t_commrec *cr,t_fftgrid *grid,int direction)
{
    static bool bFirst=TRUE;
    static real *tmp;
    int i;
    static int localsize;
    static int maxproc;
    
#ifdef GMX_MPI
    if(bFirst) {
        localsize=grid->la12r*grid->pfft.local_nx;
        if(!grid->workspace) {
            snew(tmp,localsize);
        }
        maxproc=grid->nx/grid->pfft.local_nx;
    }
    /* NOTE: FFTW doesnt necessarily use all processors for the fft;
     * above I assume that the ones that do have equal amounts of data.
     * this is bad since its not guaranteed by fftw, but works for now...
     * This will be fixed in the next release.
     */
    bFirst=FALSE;
    if (grid->workspace) {
        tmp=grid->workspace;
    }
    if (direction == GMX_SUM_QGRID_FORWARD) { 
        /* sum contributions to local grid */
        
        GMX_BARRIER(cr->mpi_comm_mygroup);
        GMX_MPE_LOG(ev_reduce_start);
        for(i=0;i<maxproc;i++) {
            MPI_Reduce(grid->ptr+i*localsize, /*ptr arithm.     */
                       tmp,localsize,      
                       GMX_MPI_REAL,MPI_SUM,i,cr->mpi_comm_mygroup);
        }
        GMX_MPE_LOG(ev_reduce_finish);
        
        if(cr->nodeid<maxproc) {
            memcpy(grid->ptr+cr->nodeid*localsize,tmp,localsize*sizeof(real));
        }
    }
    else if (direction == GMX_SUM_QGRID_BACKWARD) { 
    /* distribute local grid to all processors */
        for(i=0;i<maxproc;i++)
            MPI_Bcast(grid->ptr+i*localsize, /* ptr arithm     */
                      localsize,       
                      GMX_MPI_REAL,i,cr->mpi_comm_mygroup);
    }
    else {
        gmx_fatal(FARGS,"Invalid direction %d for summing qgrid",direction);
    }
    
#else
    gmx_fatal(FARGS,"Parallel grid summation requires MPI and FFTW.\n");    
#endif
}

static void spread_q_bsplines(gmx_pme_t pme, pme_atomcomm_t *atc, 
                              t_fftgrid *grid)
{
  /* spread charges from home atoms to local grid */
    real     *ptr;
    pme_overlap_t *ol;
    int      b,i,nn,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
    int      nx,ny,nz,nx2,ny2,nz2,la2,la12;
    int      order,norder,*idxptr,index_x,index_xy,index_xyz;
    real     valx,valxy,qn;
    real     *thx,*thy,*thz;
    int localsize, bndsize;
  
    if (pme->ndecompdim == 0) {
        clear_fftgrid(grid); 
#ifdef GMX_MPI
    } else {
        localsize = grid->la12r*grid->pfft.local_nx;
        ptr = grid->ptr + grid->la12r*grid->pfft.local_x_start;
        for (i=0; (i<localsize); i++) {
            ptr[i] = 0;
        }
        ol = &pme->overlap[0];
        /* clear left boundary area */
        for(b=0; b<ol->nleftbnd; b++) {
            ptr     = grid->ptr + ol->leftc[b].snd0*grid->la12r;
            bndsize =             ol->leftc[b].snds*grid->la12r;
            for (i=0; (i<bndsize); i++) {
                ptr[i] = 0;
            }
        }
        /* clear right boundary area */
        for(b=0; b<ol->nrightbnd; b++) {
            ptr     = grid->ptr + ol->rightc[b].snd0*grid->la12r;
            bndsize =             ol->rightc[b].snds*grid->la12r;
            for (i=0; (i<bndsize); i++) {
                ptr[i] = 0;
            }
        }
#endif
    }

    unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);

    order = pme->pme_order;

    ii0   = pme->nnx + nx2 + 1 - order/2;
    jj0   = pme->nny + ny2 + 1 - order/2;
    kk0   = pme->nnz + nz2 + 1 - order/2;
    
    for(nn=0; (nn<atc->n); nn++) {
        n = nn;
        qn     = atc->q[n];
        idxptr = atc->idx[n];
        
        if (qn != 0) {
            norder  = n*pme->pme_order;
            
            /* Pointer arithmetic alert, next six statements */
            i0  = ii0 + idxptr[XX]; 
            j0  = jj0 + idxptr[YY];
            k0  = kk0 + idxptr[ZZ];
            thx = atc->theta[XX] + norder;
            thy = atc->theta[YY] + norder;
            thz = atc->theta[ZZ] + norder;
            
            for(ithx=0; (ithx<order); ithx++) {
                index_x = la12*i0[ithx];
                valx    = qn*thx[ithx];
                
                for(ithy=0; (ithy<order); ithy++) {
                    valxy    = valx*thy[ithy];
                    index_xy = index_x+la2*j0[ithy];
                    
                    for(ithz=0; (ithz<order); ithz++) {
                        index_xyz       = index_xy+k0[ithz];
                        ptr[index_xyz] += valxy*thz[ithz];
                    }
                }
            }
        }
    }
}

real solve_pme(gmx_pme_t pme,t_fftgrid *grid,
               real ewaldcoeff,real vol,matrix vir,t_commrec *cr)
{
    /* do recip sum over local cells in grid */
    t_complex *ptr,*p0;
    int     nx,ny,nz,nx2,ny2,nz2,la2,la12;
    int     kx,ky,kz,maxkx,maxky,maxkz,kystart=0,kyend=0;
    real    m2,mx,my,mz;
    real    factor=M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    real    ets2,struct2,vfactor,ets2vf;
    real    eterm,d1,d2,energy=0;
    real    denom;
    real    bx,by;
    real    mhx,mhy,mhz;
    real    virxx=0,virxy=0,virxz=0,viryy=0,viryz=0,virzz=0;
    real    rxx,ryx,ryy,rzx,rzy,rzz;
    
    unpack_fftgrid(grid,&nx,&ny,&nz,
                   &nx2,&ny2,&nz2,&la2,&la12,FALSE,(real **)&ptr);
    
    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];
    
    maxkx = (nx+1)/2;
    maxky = (ny+1)/2;
    maxkz = nz/2+1;
    
    if (pme->ndecompdim > 0) { 
        /* transpose X & Y and only sum local cells */
#ifdef GMX_MPI
        kystart = grid->pfft.local_y_start_after_transpose;
        kyend   = kystart+grid->pfft.local_ny_after_transpose;
        if (debug)
            fprintf(debug,"solve_pme: kystart = %d, kyend=%d\n",kystart,kyend);
#else
        gmx_fatal(FARGS,"Parallel PME attempted without MPI and FFTW");
#endif /* end of parallel case loop */
    }
    else {
        kystart = 0;
        kyend   = ny;
    }
    
    for(ky=kystart; (ky<kyend); ky++) {  /* our local cells */
        
        if(ky<maxky) {
            my = ky;
        } else {
            my = (ky-ny);
        }
        by = M_PI*vol*pme->bsp_mod[YY][ky];
        
        for(kx=0; (kx<nx); kx++) {    
            if(kx < maxkx) {
                mx = kx;
            } else {
                mx = (kx-nx);
            }

            mhx = mx * rxx;
            mhy = mx * ryx + my * ryy;
            
            bx = pme->bsp_mod[XX][kx];
            
            if (pme->nnodes > 1) {
                p0 = ptr + INDEX(ky,kx,0); /* Pointer Arithmetic */
            } else {
                p0 = ptr + INDEX(kx,ky,0); /* Pointer Arithmetic */
            }
            
            for(kz=0; (kz<maxkz); kz++,p0++)  {
                if ((kx==0) && (ky==0) && (kz==0)) {
                    continue;
                }
                d1      = p0->re;
                d2      = p0->im;
                mz      = kz;
                
                mhz = mx * rzx + my * rzy + mz * rzz;
                
                m2      = mhx*mhx+mhy*mhy+mhz*mhz;
                denom   = m2*bx*by*pme->bsp_mod[ZZ][kz];
                eterm   = ONE_4PI_EPS0*exp(-factor*m2)/(pme->epsilon_r*denom);
                p0->re  = d1*eterm;
                p0->im  = d2*eterm;
                
                struct2 = d1*d1+d2*d2;
                if ((kz > 0) && (kz < (nz+1)/2)) {
                    struct2 *= 2;
                }
                ets2     = eterm*struct2;
                vfactor  = (factor*m2+1)*2.0/m2;
                energy  += ets2;
                
                ets2vf   = ets2*vfactor;
                virxx   += ets2vf*mhx*mhx-ets2;
                virxy   += ets2vf*mhx*mhy;   
                virxz   += ets2vf*mhx*mhz;  
                viryy   += ets2vf*mhy*mhy-ets2;
                viryz   += ets2vf*mhy*mhz;
                virzz   += ets2vf*mhz*mhz-ets2;
            }
        }
    }
    
    /* Update virial with local values. The virial is symmetric by definition.
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
    return(0.5*energy);
}

void gather_f_bsplines(gmx_pme_t pme,t_fftgrid *grid,
                       bool bClearF,pme_atomcomm_t *atc,real scale)
{
    /* sum forces for local particles */  
    int     nn,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
    int     nx,ny,nz,nx2,ny2,nz2,la2,la12,index_x,index_xy;
    real *  ptr;
    real    tx,ty,dx,dy,qn;
    real    fx,fy,fz,gval;
    real    fxy1,fz1;
    real    *thx,*thy,*thz,*dthx,*dthy,*dthz;
    int     norder,*idxptr;
    real    rxx,ryx,ryy,rzx,rzy,rzz;
    int     order;
    
    unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
    order = pme->pme_order;
    thx   = atc->theta[XX];
    thy   = atc->theta[YY];
    thz   = atc->theta[ZZ];
    dthx  = atc->dtheta[XX];
    dthy  = atc->dtheta[YY];
    dthz  = atc->dtheta[ZZ];
    ii0   = pme->nnx + nx2 + 1 - order/2;
    jj0   = pme->nny + ny2 + 1 - order/2;
    kk0   = pme->nnz + nz2 + 1 - order/2;
    
    rxx   = pme->recipbox[XX][XX];
    ryx   = pme->recipbox[YY][XX];
    ryy   = pme->recipbox[YY][YY];
    rzx   = pme->recipbox[ZZ][XX];
    rzy   = pme->recipbox[ZZ][YY];
    rzz   = pme->recipbox[ZZ][ZZ];

    for(nn=0; (nn<atc->n); nn++) {
        n = nn;
        qn      = scale*atc->q[n];
        
        if (bClearF) {
            atc->f[n][XX] = 0;
            atc->f[n][YY] = 0;
            atc->f[n][ZZ] = 0;
        }
        if (qn != 0) {
            fx     = 0;
            fy     = 0;
            fz     = 0;
            idxptr = atc->idx[n];
            norder = n*order;
            
            /* Pointer arithmetic alert, next nine statements */
            i0   = ii0 + idxptr[XX]; 
            j0   = jj0 + idxptr[YY];
            k0   = kk0 + idxptr[ZZ];
            thx  = atc->theta[XX] + norder;
            thy  = atc->theta[YY] + norder;
            thz  = atc->theta[ZZ] + norder;
            dthx = atc->dtheta[XX] + norder;
            dthy = atc->dtheta[YY] + norder;
            dthz = atc->dtheta[ZZ] + norder;
            
            for(ithx=0; (ithx<order); ithx++) {
                index_x = la12*i0[ithx];
                tx      = thx[ithx];
                dx      = dthx[ithx];
                
                for(ithy=0; (ithy<order); ithy++) {
                    index_xy = index_x+la2*j0[ithy];
                    ty       = thy[ithy];
                    dy       = dthy[ithy];
                    fxy1     = fz1 = 0;
                    
                    for(ithz=0; (ithz<order); ithz++) {
                        gval  = ptr[index_xy+k0[ithz]];
                        fxy1 += thz[ithz]*gval;
                        fz1  += dthz[ithz]*gval;
                    }
                    fx += dx*ty*fxy1;
                    fy += tx*dy*fxy1;
                    fz += tx*ty*fz1; 
                } 
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


void make_bsplines(splinevec theta,splinevec dtheta,int order,int nx,int ny,
		   int nz,rvec fractx[],int nr,real charge[],
		   bool bFreeEnergy)
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
  fprintf(log,"Destroying PME data structures.\n");
  sfree((*pmedata)->nnx);
  sfree((*pmedata)->nny);
  sfree((*pmedata)->nnz);
  done_fftgrid((*pmedata)->gridA);
  if((*pmedata)->gridB)
    done_fftgrid((*pmedata)->gridB);
    
  sfree(*pmedata);
  *pmedata = NULL;
  
  return 0;
}

int pme_inconvenient_nnodes(int nkx,int nky,int nnodes)
{
  int   nnx,nny;
  float imbal;
  int   ret;

  ret = 0;
  if (nnodes > nkx && nnodes > nky) {
    /* This is probably always bad */
    ret = 2;
  } else if (2*nnodes > nkx && nnodes != nkx) {
    /* This is inconvenient for the grid overlap communication */
    ret = 1;
  } 

  /* Determine the maximum number of grid slabs per PME node */
  nnx = (nkx + nnodes - 1)/nnodes;
  nny = (nky + nnodes - 1)/nnodes;
  /* Estimate the FFT + solve_pme load imbalance.
   * Imbalance in x for 2D FFT.
   * Imbalance in y for 1D FFT + solve_pme.
   * x and y imbalance affect the performance roughly equally.
   */
  imbal = (nnx*nnodes/(float)nkx + nny*nnodes/(float)nky)*0.5 - 1;
  if (debug)
    fprintf(debug,"PME load imbalance estimate for npme=%d: %f\n",
	    nnodes,imbal);

  /* The cost of charge spreading and force gathering (which is always
   * load balanced) is usually 1-2 times more than FFT+solve_pme.
   * So we compare the imbalance to (a rough guess of) the performance gain
   * in spreading and gathering with respect to one node less.
   */
  if (imbal > 2.0/nnodes) {
    ret = max(ret,2);
  } else if (imbal > 1.0/nnodes) {
    ret = max(ret,1);
  }

  return ret;
}

static void init_atomcomm(gmx_pme_t pme,pme_atomcomm_t *atc,bool bSpread)
{
    int lbnd,rbnd,maxlr,b,i;
    int nn,nk;
    pme_grid_comm_t *pgc;

#ifdef GMX_MPI
    atc->mpi_comm = pme->mpi_comm;
#endif

    atc->nslab   = pme->nnodes;
    atc->nodeid  = pme->nodeid;

    atc->bSpread   = bSpread;
    atc->pme_order = pme->pme_order;
}

static void init_overlap_comm(gmx_pme_t pme,pme_overlap_t *ol)
{
    int lbnd,rbnd,maxlr,b,i;
    int nn,nk;
    pme_grid_comm_t *pgc;

#ifdef GMX_MPI
    ol->mpi_comm = pme->mpi_comm;
#endif

    ol->nodeid = pme->nodeid;
    
    nn = pme->nnodes;
    nk = pme->nkx;

    /* Determine the grid boundary communication sizes and nodes */
    if (nk % nn == 0) {
        lbnd = pme->pme_order/2 - 1;
    } else {
        lbnd = pme->pme_order - (pme->pme_order/2 - 1) - 1;
    }
    rbnd     = pme->pme_order - (pme->pme_order/2 - 1) - 1;
    /* Round up */
    ol->nleftbnd  = (lbnd*nn + nk - 1)/nk;
    ol->nrightbnd = (rbnd*nn + nk - 1)/nk;
    maxlr = max(ol->nleftbnd,ol->nrightbnd);
    snew(ol->leftid, maxlr);
    snew(ol->rightid,maxlr);
    for(b=0; b<maxlr; b++) {
        ol->leftid[b]  = (ol->nodeid - (b + 1) + nn) % nn;
        ol->rightid[b] = (ol->nodeid + (b + 1)) % nn;
    }
    snew(ol->leftc, ol->nleftbnd);
    snew(ol->rightc,ol->nrightbnd);
    snew(ol->s2g,nn+1);
    for(i=0; i<nn+1; i++) {
        /* The definition of the grid position requires rounding up here */
        ol->s2g[i] = (i*nk + nn - 1)/nn;
    }
    /* The left boundary */
    for(b=0; b<ol->nleftbnd; b++) {
        pgc = &ol->leftc[b];
        /* Send */
        i = ol->s2g[ol->nodeid];
        if (ol->leftid[b] > ol->nodeid) {
            i += pme->nkx;
        }
        pgc->snd0 = max(i - lbnd,ol->s2g[ol->leftid[b]]);
        pgc->snds = min(i       ,ol->s2g[ol->leftid[b]+1]) - pgc->snd0;
        pgc->snds = max(pgc->snds,0);
        /* Receive */
        i = ol->s2g[ol->rightid[b]];
        if (ol->rightid[b] < ol->nodeid)
            i += pme->nkx;
        pgc->rcv0 = max(i - lbnd,ol->s2g[ol->nodeid]);
        pgc->rcvs = min(i       ,ol->s2g[ol->nodeid+1]) - pgc->rcv0;
        pgc->rcvs = max(pgc->rcvs,0);
    }
    /* The right boundary */
    for(b=0; b<ol->nrightbnd; b++) {
        pgc = &ol->rightc[b];
        /* Send */
        i = ol->s2g[ol->nodeid+1];
        if (ol->rightid[b] < ol->nodeid) {
            i -= nk;
        }
        pgc->snd0 = max(i       ,ol->s2g[ol->rightid[b]]);
        pgc->snds = min(i + rbnd,ol->s2g[ol->rightid[b]+1]) - pgc->snd0;
        pgc->snds = max(pgc->snds,0);
        /* Receive */
        i = ol->s2g[ol->leftid[b]+1];
        if (ol->leftid[b] > ol->nodeid)
            i -= pme->nkx;
        pgc->rcv0 = max(i       ,ol->s2g[ol->nodeid]);
        pgc->rcvs = min(i + rbnd,ol->s2g[ol->nodeid+1]) - pgc->rcv0;
        pgc->rcvs = max(pgc->rcvs,0);
    }
}

int gmx_pme_init(gmx_pme_t *pmedata,t_commrec *cr,
                 t_inputrec *ir,int homenr,
                 bool bFreeEnergy,
                 bool bReproducible)
{
    gmx_pme_t pme=NULL;
    
    pme_atomcomm_t *atc;
    int b,d,i,lbnd,rbnd,maxlr;
    
    if (debug)
        fprintf(debug,"Creating PME data structures.\n");
    snew(pme,1);
    
    if (PAR(cr)) {
        pme->ndecompdim = 1;
#ifdef GMX_MPI
        pme->mpi_comm = cr->mpi_comm_mygroup;
        MPI_Comm_rank(pme->mpi_comm,&pme->nodeid);
        MPI_Comm_size(pme->mpi_comm,&pme->nnodes);
#endif
        pme->bPPnode = (cr->duty & DUTY_PP);
    } else {
        pme->ndecompdim = 0;
        pme->nnodes = 1;
        pme->bPPnode = TRUE;
    }
    
    if (ir->ePBC == epbcSCREW) {
        gmx_fatal(FARGS,"pme does not (yet) work with pbc = screw");
    }
    
    pme->bFEP = ((ir->efep != efepNO) && bFreeEnergy);
    pme->nkx  = ir->nkx;
    pme->nky  = ir->nky;
    pme->nkz  = ir->nkz;
    pme->pme_order   = ir->pme_order;
    pme->epsilon_r   = ir->epsilon_r;
    
    /* Use atc[0] for spreading */
    init_atomcomm(pme,&pme->atc[0],TRUE);
    
    if (pme->nkx <= pme->pme_order*(pme->nnodes > 1 ? 2 : 1) ||
        pme->nky <= pme->pme_order ||
        pme->nkz <= pme->pme_order)
        gmx_fatal(FARGS,"The pme grid dimensions need to be larger than pme_order (%d) and in parallel larger than 2*pme_order for x",pme->pme_order);
    
    if (pme->nnodes > 1) {
#ifdef GMX_MPI
        MPI_Type_contiguous(DIM, mpi_type, &rvec_mpi);
        MPI_Type_commit(&rvec_mpi);
#endif
        
        /* Note that the charge spreading and force gathering, which usually
         * takes about the same amount of time as FFT+solve_pme,
         * is always fully load balanced
         * (unless the charge distribution is inhomogeneous).
         */
        
        if (pme_inconvenient_nnodes(pme->nkx,pme->nky,pme->nnodes) &&
            pme->nodeid == 0) {
            fprintf(stderr,
                    "\n"
                    "NOTE: For optimal PME load balancing at high parallelization\n"
                    "      PME grid_x (%d) and grid_y (%d) should be divisible by #PME_nodes (%d)\n"
                    "\n",
                    pme->nkx,pme->nky,pme->nnodes);
        }
        
        if (debug) {
            fprintf(debug,"Parallelized PME sum used. nkx=%d, npme=%d\n",
                    ir->nkx,pme->nnodes);
            if ((ir->nkx % pme->nnodes) != 0)
                fprintf(debug,"Warning: For load balance, fourier_nx should be divisible by the number of PME nodes\n");
        }
        
        atc = &pme->atc[0];
        if (DOMAINDECOMP(cr)) {
            snew(atc->node_dest,pme->nnodes);
            snew(atc->node_src,pme->nnodes);
            setup_coordinate_communication(atc);
        }
        snew(atc->count,pme->nnodes);
        snew(atc->rcount,pme->nnodes);
        snew(atc->buf_index,pme->nnodes);
        
        init_overlap_comm(pme,&pme->overlap[0]);
    } else {
        pme->overlap[0].s2g = NULL;
    }
    
    /* With domain decomposition we need nnx on the PP only nodes */
    snew(pme->nnx,5*pme->nkx);
    snew(pme->nny,5*pme->nky);
    snew(pme->nnz,5*pme->nkz);
    for(i=0; (i<5*pme->nkx); i++) {
        pme->nnx[i] = i % pme->nkx;
    }
    for(i=0; (i<5*pme->nky); i++) {
        pme->nny[i] = i % pme->nky;
    }
    for(i=0; (i<5*pme->nkz); i++) {
        pme->nnz[i] = i % pme->nkz;
    }
    
    snew(pme->bsp_mod[XX],pme->nkx);
    snew(pme->bsp_mod[YY],pme->nky);
    snew(pme->bsp_mod[ZZ],pme->nkz);
    
    pme->gridA = mk_fftgrid(pme->nkx,pme->nky,pme->nkz,
                            NULL,pme->overlap[0].s2g,cr,
                            bReproducible);
    if (bFreeEnergy) {
        pme->gridB = mk_fftgrid(pme->nkx,pme->nky,pme->nkz,
                                NULL,pme->overlap[0].s2g,cr,
                                bReproducible);
    } else {
        pme->gridB = NULL;
    }
    
    make_bspline_moduli(pme->bsp_mod,pme->nkx,pme->nky,pme->nkz,pme->pme_order);
    
    if (pme->nnodes == 1) {
        pme->atc[0].n = homenr;
        pme_realloc_atomcomm_things(&pme->atc[0]);
    }
    
    *pmedata = pme;
    
    return 0;
}

static void spread_on_grid(gmx_pme_t pme, pme_atomcomm_t *atc,
                           t_fftgrid *grid, t_commrec *cr,    
                           matrix box,
                           bool bGatherOnly, bool bHaveSplines)
{ 
    int nx,ny,nz,nx2,ny2,nz2,la2,la12;
    real *ptr;
    
    /* Unpack structure */
    unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
    
    if (!bHaveSplines) {
        /* Inverse box */
        calc_recipbox(box,pme->recipbox); 
    }
  
    if (!bGatherOnly) {
        if (!bHaveSplines) {
            /* Compute fftgrid index for all atoms,
             * with help of some extra variables.
             */
            calc_idx(pme,atc->x);
            
            /* make local bsplines  */
            make_bsplines(atc->theta,atc->dtheta,pme->pme_order,nx,ny,nz,
                          atc->fractx,atc->n,atc->q,pme->bFEP);
        }    
        
        /* put local atoms on grid. */
        spread_q_bsplines(pme,atc,grid);
        /*    pr_grid_dist(logfile,"spread",grid); */
    }
}

int gmx_pmeonly(gmx_pme_t pme,
                t_commrec *cr,    t_nrnb *nrnb,
                gmx_wallcycle_t wcycle,
                real ewaldcoeff,  bool bGatherOnly)
{
    gmx_pme_pp_t pme_pp;
    int  natoms;
    matrix box;
    rvec *x_pp=NULL,*f_pp=NULL;
    real *chargeA=NULL,*chargeB=NULL;
    real lambda=0;
    int  maxshift=0;
    real energy,dvdlambda;
    matrix vir;
    float cycles;
    int  count;
    
    pme_pp = gmx_pme_pp_init(cr);
    
    init_nrnb(nrnb);
    
    count = 0;
    do /****** this is a quasi-loop over time steps! */
    {
        /* Domain decomposition */
        natoms = gmx_pme_recv_q_x(pme_pp,
                                  &chargeA,&chargeB,box,&x_pp,&f_pp,
                                  &maxshift,&pme->bFEP,&lambda);
        
        if (natoms == -1) {
            /* We should stop: break out of the loop */
            break;
        }
        
        if (count == 0)
            wallcycle_start(wcycle,ewcRUN);
        
        wallcycle_start(wcycle,ewcPMEMESH_SEP);
        
        dvdlambda = 0;
        clear_mat(vir);
        gmx_pme_do(pme,0,natoms,x_pp,f_pp,chargeA,chargeB,box,
                   cr,maxshift,nrnb,vir,ewaldcoeff,
                   &energy,lambda,&dvdlambda,
                   bGatherOnly);
        
        cycles = wallcycle_stop(wcycle,ewcPMEMESH_SEP);
        
        gmx_pme_send_force_vir_ener(pme_pp,
                                    f_pp,vir,energy,dvdlambda,
                                    cycles,bGotTermSignal,bGotUsr1Signal);
        
        count++;
        
        /* MPI_Barrier(cr->mpi_comm_mysim); */ /* 100 */
    } /***** end of quasi-loop, we stop with the break above */
    while (TRUE);
    
    return 0;
}

int gmx_pme_do(gmx_pme_t pme,
               int start,       int homenr,
               rvec x[],        rvec f[],
               real *chargeA,   real *chargeB,
               matrix box,	t_commrec *cr,
               int  maxshift,   t_nrnb *nrnb,    
               matrix vir,      real ewaldcoeff,
               real *energy,    real lambda, 
               real *dvdlambda, bool bGatherOnly)
{
    int     q,i,j,ntot,npme;
    int     nx,ny,nz,nx2,ny2,nz2,la12,la2;
    int     local_ny;
    pme_atomcomm_t *atc;
    t_fftgrid *grid=NULL;
    real    *ptr;
    real    *charge=NULL,vol;
    real    energy_AB[2];
    matrix  vir_AB[2];
    bool    bClearF;

    if (pme->nnodes > 1) {
        atc = &pme->atc[0];
        atc->npd = homenr;
        if (atc->npd > atc->pd_nalloc) {
            atc->pd_nalloc = over_alloc_dd(atc->npd);
            srenew(atc->pd,atc->pd_nalloc);
        }
    }
    
    for(q=0; q<(pme->bFEP ? 2 : 1); q++) {
        if (q == 0) {
            grid = pme->gridA;
            charge = chargeA+start;
        } else {
            grid = pme->gridB;
            charge = chargeB+start;
        }
        /* Unpack structure */
        if (debug) {
            fprintf(debug,"PME: nnodes = %d, nodeid = %d\n",
                    cr->nnodes,cr->nodeid);
            fprintf(debug,"Grid = %p\n",grid);
            if (grid == NULL)
                gmx_fatal(FARGS,"No grid!");
        }
        where();
        unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
#ifdef GMX_MPI
        if (pme->nnodes > 1) {
            local_ny = grid->pfft.local_ny_after_transpose;
        } else {
            local_ny = ny;
        }
#else
        local_ny = ny;
#endif
        where();
        
        atc = &pme->atc[0];
        if (pme->nnodes == 1) {
            if (DOMAINDECOMP(cr)) {
                atc->n = homenr;
                pme_realloc_atomcomm_things(atc);
            }
            atc->x = x;
            atc->q = charge;
            atc->f = f;
        } else {
            pme_calc_pidx(pme->nnodes,homenr,box,x+start,atc);
            where();
            
            GMX_BARRIER(cr->mpi_comm_mysim);
            /* Redistribute x (only once) and qA or qB */
            if (DOMAINDECOMP(cr)) {
                dd_pmeredist_x_q(pme, maxshift, homenr, q==0, x, charge, atc);
            } else {
                pmeredist(pme, TRUE, homenr, q==0, x+start, charge, atc);
            }
            where();
        }
        
        if (debug)
            fprintf(debug,"Node= %6d, pme local particles=%6d\n",
                    cr->nodeid,atc->n);
        
        /* Spread the charges on a grid */
        GMX_MPE_LOG(ev_spread_on_grid_start);
        
        /* Spread the charges on a grid */
        spread_on_grid(pme,&pme->atc[0],grid,cr,box,bGatherOnly,
                       q==0 ? FALSE : TRUE);
        GMX_MPE_LOG(ev_spread_on_grid_finish);
        
        if (!bGatherOnly) {
            if (q == 0) {
                inc_nrnb(nrnb,eNR_WEIGHTS,DIM*atc->n);
            }
            inc_nrnb(nrnb,eNR_SPREADQBSP,
                     pme->pme_order*pme->pme_order*pme->pme_order*atc->n);
            
            /* sum contributions to local grid from other nodes */
            if (pme->nnodes > 1) {
#ifdef DEBUG
                if (debug)
                    pr_fftgrid(debug,"qgrid before dd sum",grid);
#endif
                GMX_BARRIER(cr->mpi_comm_mysim);
                gmx_sum_qgrid_dd(&pme->overlap[0],grid,GMX_SUM_QGRID_FORWARD);
                where();
            }
#ifdef DEBUG
            if (debug)
                pr_fftgrid(debug,"qgrid",grid);
#endif
            where();
            
            /* do 3d-fft */ 
            GMX_BARRIER(cr->mpi_comm_mysim);
            GMX_MPE_LOG(ev_gmxfft3d_start);
            gmxfft3D(grid,GMX_FFT_REAL_TO_COMPLEX,cr);
            GMX_MPE_LOG(ev_gmxfft3d_finish);
            where();
            
            /* solve in k-space for our local cells */
            vol = det(box);
            GMX_BARRIER(cr->mpi_comm_mysim);
            GMX_MPE_LOG(ev_solve_pme_start);
            energy_AB[q]=solve_pme(pme,grid,ewaldcoeff,vol,vir_AB[q],cr);
            where();
            GMX_MPE_LOG(ev_solve_pme_finish);
            inc_nrnb(nrnb,eNR_SOLVEPME,nx*local_ny*(nz/2+1));
            
            /* do 3d-invfft */
            GMX_BARRIER(cr->mpi_comm_mysim);
            GMX_MPE_LOG(ev_gmxfft3d_start);
            where();
            gmxfft3D(grid,GMX_FFT_COMPLEX_TO_REAL,cr);
            where();
            GMX_MPE_LOG(ev_gmxfft3d_finish);
            
            /* distribute local grid to all nodes */
            if (pme->nnodes > 1) {
                GMX_BARRIER(cr->mpi_comm_mysim);
                gmx_sum_qgrid_dd(&pme->overlap[0],grid,GMX_SUM_QGRID_BACKWARD);
            }
            where();
#ifdef DEBUG
            if (debug)
                pr_fftgrid(debug,"potential",grid);
#endif
            
            ntot  = grid->nxyz;  
            npme  = ntot*log((real)ntot)/log(2.0);
            if (pme->nnodes > 1) {
                npme /= (cr->nnodes - cr->npmenodes);
            }
            inc_nrnb(nrnb,eNR_FFT,2*npme);
            where();
        }
        /* interpolate forces for our local atoms */
        GMX_BARRIER(cr->mpi_comm_mysim);
        GMX_MPE_LOG(ev_gather_f_bsplines_start);
        
        where();
        /* If are running without parallelization,
         * atc->f is the actual force array, not a buffer,
         * therefore we should not clear it.
         */
        bClearF = (q == 0 && pme->ndecompdim > 0);
        gather_f_bsplines(pme,grid,bClearF,&pme->atc[0],
                          pme->bFEP ? (q==0 ? 1.0-lambda : lambda) : 1.0);
        where();
        
        GMX_MPE_LOG(ev_gather_f_bsplines_finish);
        
        inc_nrnb(nrnb,eNR_GATHERFBSP,
                 pme->pme_order*pme->pme_order*pme->pme_order*pme->atc[0].n);
    } /* of q-loop */
    
    if (pme->nnodes > 1) {
        GMX_BARRIER(cr->mpi_comm_mysim);
        if (DOMAINDECOMP(cr)) {
            dd_pmeredist_f(pme,&pme->atc[0],maxshift,homenr,f,pme->bPPnode);
        } else {
            pmeredist(pme, FALSE, homenr, TRUE, f+start, NULL, &pme->atc[0]);
        }
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
