/*
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
#include "shift_util.h"
#include "ewald_util.h"
#include "fftgrid.h"
#include "gmx_fatal.h"
#include "ewald.h"
#include "pme.h"
#include "network.h"
#include "physics.h"
#include "nrnb.h"
#include "copyrite.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif

#ifdef USE_MPE
#include "mpe.h"
#include "mpelogging.h"
#endif

#define DFT_TOL 1e-7

#ifdef GMX_DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif

#ifdef GMX_MPI
MPI_Datatype  rvec_mpi;
#endif

/* Global variables */

static bool bPar;               /* parallel? */
static int leftnode, rightnode; /* rank of left and right MPI task */
static int leftbnd, rightbnd;   /* size of left and right boundary */
static int my_homenr;

void pr_grid_dist(FILE *fp,char *title,t_fftgrid *grid)
{
  int     i,j,k,l,ntoti,ntot=0;
  int     nx,ny,nz,nx2,ny2,nz2,la12,la2;
  real *ptr;

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

void calc_recipbox(matrix box,matrix recipbox)
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


void calc_idx(int natoms,matrix recipbox,
	      rvec x[],rvec fractx[],ivec idx[],int nx,int ny,int nz,
	      int nx2,int ny2,int nz2,
	      int nnx[],int nny[],int nnz[])
{
  int  i;
  int  *idxptr,tix,tiy,tiz;
  real *xptr,tx,ty,tz;
  real rxx,ryx,ryy,rzx,rzy,rzz;
#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined GMX_DOUBLE && defined GMX_X86TRUNC)
  int x86_cw,x86_cwsave;

  asm("fnstcw %0" : "=m" (*&x86_cwsave));
  x86_cw = x86_cwsave | 3072;
  asm("fldcw %0" : : "m" (*&x86_cw));
  #define x86trunc(a,b) asm("fld %1\nfistpl %0\n" : "=m" (*&b) : "f" (a));
#endif
 
  rxx = recipbox[XX][XX];
  ryx = recipbox[YY][XX];
  ryy = recipbox[YY][YY];
  rzx = recipbox[ZZ][XX];
  rzy = recipbox[ZZ][YY];
  rzz = recipbox[ZZ][ZZ];

  for(i=0; (i<natoms); i++) {
    xptr   = x[i];
    idxptr = idx[i];
    
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

    fractx[i][XX] = tx - tix;
    fractx[i][YY] = ty - tiy;
    fractx[i][ZZ] = tz - tiz;   
    
    idxptr[XX] = nnx[tix];
    idxptr[YY] = nny[tiy];
    idxptr[ZZ] = nnz[tiz];
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

void pme_calc_pidx(t_commrec *cr,
		   int natoms,matrix box, rvec x[],int nx, int nnx[],
                   int *pidx, int *gidx)
{
  int  i,j,jj;
  int nx2;
  int  tix;
  int  irange, taskid;
  real *xptr,tx;
  real rxx,ryx,rzx;
  matrix recipbox;

/* Calculate PME task index (pidx) for each grid index */
  irange=1+(nx-1)/cr->nnodes;

  taskid=0;
  for(i=0; i<nx; i+=irange) {
    jj=min((i+irange),nx);
    for(j=i; (j<jj); j++) pidx[j]=taskid;
    taskid++;
  }

  nx2=2*nx;
  calc_recipbox(box,recipbox);
  rxx = recipbox[XX][XX];
  ryx = recipbox[YY][XX];
  rzx = recipbox[ZZ][XX];
/* Calculate grid index in x-dimension */
  for(i=0; (i<natoms); i++) {
    xptr   = x[i];
    /* Fractional coordinates along box vectors */
    tx = nx2 + nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx );
    tix = (int)(tx);
    gidx[i] = pidx[nnx[tix]];
  }
}

void pmeredist(t_commrec *cr, bool forw,
	       int n, rvec *x_f, real *charge, int *idxa,
               int *total, rvec *pme_x, real *pme_q)
/* Redistribute particle data for PME calculation */
/* domain decomposition by x coordinate           */
{
  static int *scounts, *rcounts,*sdispls, *rdispls, *sidx;
  static real *buf;
  int i, ii;
  static bool bFirst=TRUE;
  int nnodes;

  nnodes = cr->nnodes;
  if(bFirst) {
    snew(scounts,nnodes);
    snew(rcounts,nnodes);
    snew(sdispls,nnodes);
    snew(rdispls,nnodes);
    snew(sidx,nnodes);
    snew(buf,n*DIM);
    bFirst=FALSE;
  }

#ifdef GMX_MPI
  if (forw) { /* forward, redistribution from pp to pme */ 

/* Calculate send counts and exchange them with other nodes */
    for(i=0; (i<nnodes); i++) scounts[i]=0;
    for(i=0; (i<n); i++) scounts[idxa[i]]++;
    MPI_Alltoall( scounts, 1, MPI_INT, rcounts, 1, MPI_INT, MPI_COMM_WORLD );

/* Calculate send and receive displacements and index into send buffer */
    sdispls[0]=0;
    rdispls[0]=0;
    sidx[0]=0;
    for(i=1; (i<nnodes); i++) {
      sdispls[i]=sdispls[i-1]+scounts[i-1];
      rdispls[i]=rdispls[i-1]+rcounts[i-1];
      sidx[i]=sdispls[i];
    }
/* Total # of particles to be received */
    *total=rdispls[nnodes-1]+rcounts[nnodes-1];

/* Copy particle coordinates into send buffer and exchange*/
    for(i=0; (i<n); i++) {
      ii=DIM*sidx[idxa[i]];
      sidx[idxa[i]]++;
      buf[ii+XX]=x_f[i][XX];
      buf[ii+YY]=x_f[i][YY];
      buf[ii+ZZ]=x_f[i][ZZ];
    }
    MPI_Alltoallv(buf, scounts, sdispls, rvec_mpi,
                  pme_x, rcounts, rdispls, rvec_mpi,
                  MPI_COMM_WORLD);

/* Copy charge into send buffer and exchange*/
    for(i=0; (i<nnodes); i++) sidx[i]=sdispls[i];
    for(i=0; (i<n); i++) {
      ii=sidx[idxa[i]];
      sidx[idxa[i]]++;
      buf[ii]=charge[i];
    }
    MPI_Alltoallv(buf, scounts, sdispls, mpi_type,
                  pme_q, rcounts, rdispls, mpi_type,
                  MPI_COMM_WORLD);

  }
  else { /* backward, redistribution from pme to pp */ 
    MPI_Alltoallv(pme_x, rcounts, rdispls, rvec_mpi,
                  buf, scounts, sdispls, rvec_mpi,
                  MPI_COMM_WORLD);

/* Copy data from receive buffer */
    for(i=0; (i<nnodes); i++) sidx[i]=sdispls[i];
    for(i=0; (i<n); i++) {
      ii=DIM*sidx[idxa[i]];
      x_f[i][XX]+=buf[ii+XX];
      x_f[i][YY]+=buf[ii+YY];
      x_f[i][ZZ]+=buf[ii+ZZ];
      sidx[idxa[i]]++;
    }
  }
#endif 
}

void sum_qgrid(t_commrec *cr,t_nsborder *nsb,t_fftgrid *grid,
               int pme_order, bool bForward)
{
  static bool bFirst=TRUE;
  static real *tmp;
  int i;
  static int localsize;
  static int maxproc;
  int ny, la2r;
  int bndsize;
  real *from, *to;
#ifdef GMX_MPI
  MPI_Status stat;
#endif
  int nodeid;

#if (defined GMX_MPI)
  nodeid = cr->nodeid;
  if(bFirst) {
    localsize=grid->la12r*grid->pfft.local_nx;
    if(!grid->workspace)
      snew(tmp,localsize);
    maxproc=grid->nx/grid->pfft.local_nx;
  }
  bFirst=FALSE;
  if(grid->workspace)
    tmp=grid->workspace;
  if(bForward) { /* sum contributions to local grid */
/* #define USE_REDUCE */
#ifdef USE_REDUCE
    /********************************/
    /* Loop with several MPI_Reduce */
    /********************************/

#ifdef USE_MPE
    MPE_Log_event( ev_reduce_start, 0, "" );
#endif
    for(i=0;i<maxproc;i++) {
      MPI_Reduce(grid->ptr+i*localsize, /*ptr arithm.     */
		 tmp,localsize,      
		 GMX_MPI_REAL,MPI_SUM,i,MPI_COMM_WORLD);
    }
#ifdef USE_MPE
    MPE_Log_event(ev_reduce_finish, 0, "" );
#endif

/*#define ALL_REDUCE */
#ifdef ALL_REDUCE
      MPI_Allreduce(grid->ptr, grid->ptr,grid->nptr,      
		 GMX_MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif
    if (nodeid < maxproc)
      memcpy(grid->ptr+nodeid*localsize,tmp,localsize*sizeof(real));
#endif
#define EXCHANGE_GRID_BOUNDARY1
#ifdef EXCHANGE_GRID_BOUNDARY1

    ny=grid->ny;
    la2r=grid->la2r;
/* Send left Boundary */
    bndsize = leftbnd*ny*la2r;
    from = grid->ptr + (leftnode + 1)*localsize - bndsize;
    to   = grid->ptr + (nodeid +   1)*localsize - bndsize;
    MPI_Sendrecv(from,bndsize,mpi_type,
               leftnode,nodeid,
               tmp,bndsize,mpi_type,rightnode,rightnode,
               MPI_COMM_WORLD,&stat);
    for(i=0; (i<bndsize); i++) {
      to[i] += tmp[i];
    }
/* Send right boundary */
    bndsize = rightbnd*ny*la2r;
    from = grid->ptr + rightnode*localsize;
    to   = grid->ptr + nodeid*localsize;
    MPI_Sendrecv(from,bndsize,mpi_type,
               rightnode,nodeid,
               tmp,bndsize,mpi_type,leftnode,leftnode,
               MPI_COMM_WORLD,&stat);
    for(i=0; (i<bndsize); i++) {
      to[i] += tmp[i];
    }
#endif
  }
  else { /* distribute local grid to all processors */
/*  #define USE_BCAST */
#ifdef USE_BCAST
#ifdef USE_MPE
    MPE_Log_event( ev_bcast_start, 0, "" );
#endif
    for(i=0;i<maxproc;i++)
      MPI_Bcast(grid->ptr+i*localsize, /* ptr arithm     */
		localsize,       
		GMX_MPI_REAL,i,MPI_COMM_WORLD);
#ifdef USE_MPE
    MPE_Log_event(ev_bcast_finish, 0, "" );
#endif
#endif
#define EXCHANGE_GRID_BOUNDARY2
#ifdef EXCHANGE_GRID_BOUNDARY2

    ny=grid->ny;
    la2r=grid->la2r;

/* Send left Boundary */
    bndsize = rightbnd*ny*la2r;
    from = grid->ptr + nodeid*localsize;
    to   = grid->ptr + rightnode*localsize;
    MPI_Sendrecv(from,bndsize,mpi_type,
               leftnode,nodeid,
               to,bndsize,mpi_type,rightnode,rightnode,
               MPI_COMM_WORLD,&stat);
/* Send right boundary */
    bndsize = leftbnd*ny*la2r;
    from = grid->ptr + (nodeid   + 1)*localsize - bndsize;
    to   = grid->ptr + (leftnode + 1)*localsize - bndsize;
    MPI_Sendrecv(from,bndsize,mpi_type,
               rightnode,nodeid,
               to,bndsize,mpi_type,leftnode,leftnode,
               MPI_COMM_WORLD,&stat);
#endif
  }
#else
  gmx_fatal(FARGS,"Parallel grid summation requires MPI.\n");    
#endif
}

void spread_q_bsplines(t_fftgrid *grid,ivec idx[],real charge[],
		       splinevec theta,int nr,int order,
		       int nnx[],int nny[],int nnz[])
{
  /* spread charges from home atoms to local grid */
  real *ptr;
  int      i,j,k,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
  int      nx,ny,nz,nx2,ny2,nz2,la2,la12,xidx,yidx,zidx;
  int      norder,norder1,*idxptr,ind0;
  real     valx,valxy,qn;
  real     *thx,*thy,*thz;
  int localsize, bndsize;
  real thx0,thx1,thx2,thx3;
  real thy0,thy1,thy2,thy3;
  real thz0,thz1,thz2,thz3;
  int     pidx;
  
  if (!bPar) {
    clear_fftgrid(grid); 
#if (defined GMX_MPI)
  } else {
    localsize = grid->la12r*grid->pfft.local_nx;
    ptr = grid->localptr;
    for (i=0; (i<localsize); i++)
      ptr[i] = 0;
/* clear left boundary area */
    bndsize = leftbnd*grid->la12r;
    ptr = grid->ptr + (leftnode + 1)*localsize - bndsize;
    for (i=0; (i<bndsize); i++)
      ptr[i] = 0;
/* clear right boundary area */
    bndsize = rightbnd*grid->la12r;
    ptr = grid->ptr + rightnode*localsize;
    for (i=0; (i<bndsize); i++)
      ptr[i] = 0;
#endif
  }
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  ii0   = nnx+nx2+1-order/2;
  jj0   = nny+ny2+1-order/2;
  kk0   = nnz+nz2+1-order/2;
  thx   = theta[XX];
  thy   = theta[YY];
  thz   = theta[ZZ];
  
  for(n=0; (n<nr); n++) {
    qn     = charge[n];
    idxptr = idx[n];
    
    if (qn != 0) {
      xidx    = idxptr[XX];
      yidx    = idxptr[YY];
      zidx    = idxptr[ZZ];
#ifdef DEBUG
      range_check(xidx,0,nx);
      range_check(yidx,0,ny);
      range_check(zidx,0,nz);
#endif
      i0      = ii0+xidx; /* Pointer arithmetic */
      norder  = n*order;
      norder1 = norder+order;

      i = ii0[xidx];
      j = jj0[yidx];
      k = kk0[zidx];
      
      if(order == 4 && i<nx-3 && j<ny-3 && k<nz-3)
      {
          thx0 = thx[norder];
          thx1 = thx[norder+1];
          thx2 = thx[norder+2];
          thx3 = thx[norder+3];
          thy0 = thy[norder];
          thy1 = thy[norder+1];
          thy2 = thy[norder+2];
          thy3 = thy[norder+3];
          thz0 = thz[norder];
          thz1 = thz[norder+1];
          thz2 = thz[norder+2];
          thz3 = thz[norder+3];

          pidx = INDEX(i,j,k);
          ptr[pidx]   += qn*thx0*thy0*thz0;
          ptr[pidx+1] += qn*thx0*thy0*thz1;
          ptr[pidx+2] += qn*thx0*thy0*thz2;
          ptr[pidx+3] += qn*thx0*thy0*thz3;
          pidx = INDEX(i,j+1,k);
          ptr[pidx]   += qn*thx0*thy1*thz0;
          ptr[pidx+1] += qn*thx0*thy1*thz1;
          ptr[pidx+2] += qn*thx0*thy1*thz2;
          ptr[pidx+3] += qn*thx0*thy1*thz3;
          pidx = INDEX(i,j+2,k);
          ptr[pidx]   += qn*thx0*thy2*thz0;
          ptr[pidx+1] += qn*thx0*thy2*thz1;
          ptr[pidx+2] += qn*thx0*thy2*thz2;
          ptr[pidx+3] += qn*thx0*thy2*thz3;
          pidx = INDEX(i,j+3,k);
          ptr[pidx]   += qn*thx0*thy3*thz0;
          ptr[pidx+1] += qn*thx0*thy3*thz1;
          ptr[pidx+2] += qn*thx0*thy3*thz2;
          ptr[pidx+3] += qn*thx0*thy3*thz3;
          pidx = INDEX(i+1,j,k);
          ptr[pidx]   += qn*thx1*thy0*thz0;
          ptr[pidx+1] += qn*thx1*thy0*thz1;
          ptr[pidx+2] += qn*thx1*thy0*thz2;
          ptr[pidx+3] += qn*thx1*thy0*thz3;
          pidx = INDEX(i+1,j+1,k);
          ptr[pidx]   += qn*thx1*thy1*thz0;
          ptr[pidx+1] += qn*thx1*thy1*thz1;
          ptr[pidx+2] += qn*thx1*thy1*thz2;
          ptr[pidx+3] += qn*thx1*thy1*thz3;
          pidx = INDEX(i+1,j+2,k);
          ptr[pidx]   += qn*thx1*thy2*thz0;
          ptr[pidx+1] += qn*thx1*thy2*thz1;
          ptr[pidx+2] += qn*thx1*thy2*thz2;
          ptr[pidx+3] += qn*thx1*thy2*thz3;
          pidx = INDEX(i+1,j+3,k);
          ptr[pidx]   += qn*thx1*thy3*thz0;
          ptr[pidx+1] += qn*thx1*thy3*thz1;
          ptr[pidx+2] += qn*thx1*thy3*thz2;
          ptr[pidx+3] += qn*thx1*thy3*thz3;
          pidx = INDEX(i+2,j,k);
          ptr[pidx]   += qn*thx2*thy0*thz0;
          ptr[pidx+1] += qn*thx2*thy0*thz1;
          ptr[pidx+2] += qn*thx2*thy0*thz2;
          ptr[pidx+3] += qn*thx2*thy0*thz3;
          pidx = INDEX(i+2,j+1,k);
          ptr[pidx]   += qn*thx2*thy1*thz0;
          ptr[pidx+1] += qn*thx2*thy1*thz1;
          ptr[pidx+2] += qn*thx2*thy1*thz2;
          ptr[pidx+3] += qn*thx2*thy1*thz3;
          pidx = INDEX(i+2,j+2,k);
          ptr[pidx]   += qn*thx2*thy2*thz0;
          ptr[pidx+1] += qn*thx2*thy2*thz1;
          ptr[pidx+2] += qn*thx2*thy2*thz2;
          ptr[pidx+3] += qn*thx2*thy2*thz3;
          pidx = INDEX(i+2,j+3,k);
          ptr[pidx]   += qn*thx2*thy3*thz0;
          ptr[pidx+1] += qn*thx2*thy3*thz1;
          ptr[pidx+2] += qn*thx2*thy3*thz2;
          ptr[pidx+3] += qn*thx2*thy3*thz3;
          pidx = INDEX(i+3,j,k);
          ptr[pidx]   += qn*thx3*thy0*thz0;
          ptr[pidx+1] += qn*thx3*thy0*thz1;
          ptr[pidx+2] += qn*thx3*thy0*thz2;
          ptr[pidx+3] += qn*thx3*thy0*thz3;
          pidx = INDEX(i+3,j+1,k);
          ptr[pidx]   += qn*thx3*thy1*thz0;
          ptr[pidx+1] += qn*thx3*thy1*thz1;
          ptr[pidx+2] += qn*thx3*thy1*thz2;
          ptr[pidx+3] += qn*thx3*thy1*thz3;
          pidx = INDEX(i+3,j+2,k);
          ptr[pidx]   += qn*thx3*thy2*thz0;
          ptr[pidx+1] += qn*thx3*thy2*thz1;
          ptr[pidx+2] += qn*thx3*thy2*thz2;
          ptr[pidx+3] += qn*thx3*thy2*thz3;
          pidx = INDEX(i+3,j+3,k);
          ptr[pidx]   += qn*thx3*thy3*thz0;
          ptr[pidx+1] += qn*thx3*thy3*thz1;
          ptr[pidx+2] += qn*thx3*thy3*thz2;
          ptr[pidx+3] += qn*thx3*thy3*thz3;
      }
    else
    {
      for(ithx=norder; (ithx<norder1); ithx++,i0++) {
	i    = *i0;
	j0   = jj0+yidx; /* Pointer arithmetic */
	valx = qn*thx[ithx];
	
	for(ithy=norder; (ithy<norder1); ithy++,j0++) {
	  j     = *j0;
	  k0    = kk0+zidx; /* Pointer arithmetic */
	  valxy = valx*thy[ithy];
	  ind0  = INDEX(i,j,0);

	  for(ithz=norder; (ithz<norder1); ithz++,k0++) {
	    k = *k0;
#ifdef DEBUG
	    range_check(i,0,nx);
	    range_check(j,0,ny);
	    range_check(k,0,nz);
	    range_check(ind0+k,0,grid->nptr);
#endif
	    ptr[ind0+k] += valxy*thz[ithz];
	  }
	}
      }
      }
    }
  }
}

real solve_pme(t_fftgrid *grid,real ewaldcoeff,real vol,real epsilon_r,
	       splinevec bsp_mod,matrix recipbox,
	       matrix vir,t_commrec *cr)
{
  /* do recip sum over local cells in grid */
  t_complex *ptr,*p0;
  int     nx,ny,nz,nx2,ny2,nz2,la2,la12;
  int     kx,ky,kz,idx,idx0,maxkx,maxky,maxkz,kystart=0,kyend=0;
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
   
  rxx = recipbox[XX][XX];
  ryx = recipbox[YY][XX];
  ryy = recipbox[YY][YY];
  rzx = recipbox[ZZ][XX];
  rzy = recipbox[ZZ][YY];
  rzz = recipbox[ZZ][ZZ];
 
  maxkx = (nx+1)/2;
  maxky = (ny+1)/2;
  maxkz = nz/2+1;
    
  if (bPar) { /* transpose X & Y and only sum local cells */
#if (defined GMX_MPI)
    kystart = grid->pfft.local_y_start_after_transpose;
    kyend   = kystart+grid->pfft.local_ny_after_transpose;
    if (debug)
      fprintf(debug,"solve_pme: kystart = %d, kyend=%d\n",kystart,kyend);
#else
    gmx_fatal(FARGS,"Parallel PME attempted without MPI.");
#endif /* end of parallel case loop */
  }
  else {
    kystart = 0;
    kyend   = ny;
  }
  
  for(ky=kystart; (ky<kyend); ky++) {  /* our local cells */
    
    if(ky<maxky)
      my = ky;
    else
      my = (ky-ny);
    by = M_PI*vol*bsp_mod[YY][ky];
    
    for(kx=0; (kx<nx); kx++) {    
      if(kx < maxkx)
	mx = kx;
      else
	mx = (kx-nx);

      mhx = mx * rxx;
      mhy = mx * ryx + my * ryy;

      bx = bsp_mod[XX][kx];
      
      if (bPar)
	p0 = ptr + INDEX(ky,kx,0); /* Pointer Arithmetic */
      else
	p0 = ptr + INDEX(kx,ky,0); /* Pointer Arithmetic */

      for(kz=0; (kz<maxkz); kz++,p0++)  {
	if ((kx==0) && (ky==0) && (kz==0))
	  continue;
	d1      = p0->re;
	d2      = p0->im;
	mz      = kz;

	mhz = mx * rzx + my * rzy + mz * rzz;

	m2      = mhx*mhx+mhy*mhy+mhz*mhz;
	denom   = m2*bx*by*bsp_mod[ZZ][kz];
	eterm   = ONE_4PI_EPS0*exp(-factor*m2)/(epsilon_r*denom);
	p0->re  = d1*eterm;
	p0->im  = d2*eterm;
	
	struct2 = d1*d1+d2*d2;
	if ((kz > 0) && (kz < (nz+1)/2))
	  struct2*=2;
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

void gather_f_bsplines(t_fftgrid *grid,matrix recipbox,
		       ivec idx[],rvec f[],real *charge,real scale,
		       splinevec theta,splinevec dtheta,
		       int nr,int order,int nnx[],int nny[],int nnz[])
{
  /* sum forces for local particles */  
  int     i,j,k,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
  int     nx,ny,nz,nx2,ny2,nz2,la2,la12;
  real *ptr;
  int     xidx,yidx,zidx;
  real    tx,ty,dx,dy,qn;
  real    fx,fy,fz,gval,tgz,dgz;
  real    gval1,gval2,gval3,gval4;
  real    fxy1,fz1;
  real    *thx,*thy,*thz,*dthx,*dthy,*dthz;
  int     sn,norder,norder1,*idxptr,ind0;
  real    rxx,ryx,ryy,rzx,rzy,rzz;

  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
 
  thx  = theta[XX];
  thy  = theta[YY];
  thz  = theta[ZZ];
  dthx = dtheta[XX];
  dthy = dtheta[YY];
  dthz = dtheta[ZZ];
  ii0  = nnx+nx2+1-order/2;
  jj0  = nny+ny2+1-order/2;
  kk0  = nnz+nz2+1-order/2;
  
  rxx = recipbox[XX][XX];
  ryx = recipbox[YY][XX];
  ryy = recipbox[YY][YY];
  rzx = recipbox[ZZ][XX];
  rzy = recipbox[ZZ][YY];
  rzz = recipbox[ZZ][ZZ];


  for(n=0; (n<nr); n++) {
    qn      = scale*charge[n];
    fx      = 0.0;
    fy      = 0.0;
    fz      = 0.0;

    if (qn != 0) {
      idxptr = idx[n];
      xidx = idxptr[XX];
      yidx = idxptr[YY];
      zidx = idxptr[ZZ];
#ifdef DEBUG
      range_check(xidx,0,nx);
      range_check(yidx,0,ny);
      range_check(zidx,0,nz);
#endif
      
      i0      = ii0+xidx;   /* Pointer arithemtic */
      norder  = n*order;
      norder1 = norder+order;
      for(ithx=norder; (ithx<norder1); ithx++,i0++) {
        i     = *i0;
        tx    = thx[ithx];
        dx    = dthx[ithx];
        j0    = jj0+yidx;   /* Pointer arithemtic */

        if (order == 4) {
          for(ithy=norder; (ithy<norder1); ithy++,j0++) {
            j     = *j0;
            ty    = thy[ithy];
            dy    = dthy[ithy];
            k0    = kk0+zidx;     /* Pointer arithemtic */
            ind0  = INDEX(i,j,0);
            gval1 = ptr[ind0+k0[0]];
            gval2 = ptr[ind0+k0[1]];
            gval3 = ptr[ind0+k0[2]];
            gval4 = ptr[ind0+k0[3]];
            
            ithz  = norder;
            
            /* First iteration */
            fxy1  = thz[ithz]*gval1;
            fz1   = dthz[ithz]*gval1;
            ithz++;
            
            /* Second iteration */
            fxy1 += thz[ithz]*gval2;
            fz1  += dthz[ithz]*gval2;
            ithz++;
            
            /* Third iteration */
            fxy1 += thz[ithz]*gval3;
            fz1  += dthz[ithz]*gval3;
            ithz++;
            
            /* Fourth iteration */
            fxy1 += thz[ithz]*gval4;
            fz1  += dthz[ithz]*gval4;
            fx    = fx+dx*ty*fxy1;
            fy    = fy+tx*dy*fxy1;
            fz    = fz+tx*ty*fz1;    
          } 
        }
	   else {
        	  for(ithy=norder; (ithy<norder1); ithy++,j0++) {
            j     = *j0;
            ty    = thy[ithy];
            dy    = dthy[ithy];
            k0    = kk0+zidx; /* Pointer arithemtic */
            ind0  = INDEX(i,j,0);
            fxy1 = fz1 = 0;
            for(ithz=norder; (ithz<norder1); ithz++,k0++) {
              k     = *k0;
#ifdef DEBUG
	      range_check(i,0,nx);
	      range_check(j,0,ny);
	      range_check(k,0,nz);
	      range_check(ind0+k,0,grid->nptr);
#endif            
              gval  = ptr[ind0+k];
              fxy1 += thz[ithz]*gval;
              fz1  += dthz[ithz]*gval;
            }
            fx += dx*ty*fxy1;
            fy += tx*dy*fxy1;
            fz += tx*ty*fz1; 
          } 
        } 
      }
      f[n][XX] -= qn*( fx*nx*rxx );
      f[n][YY] -= qn*( fx*nx*ryx + fy*ny*ryy );
      f[n][ZZ] -= qn*( fx*nx*rzx + fy*ny*rzy + fz*nz*rzz );
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
		   int nz,rvec fractx[],ivec idx[],real charge[],int nr)
{
  /* construct splines for local atoms */
  int  i,j,k,l;
  real drXX,drYY,drZZ;
  real dr,div,rcons;
  real *dataXX,*dataYY,*dataZZ;
  real *ddataXX,*ddataYY,*ddataZZ;
  real tmpX,tmpY,tmpZ,lastX,lastY,lastZ;
  real *data,*ddata,*xptr;

  if( order == 4)
  {
    for(i=0; (i<nr); i++)
	{

      if (charge[i] != 0.0) 
	  {

	    drXX = fractx[i][XX];
	    drYY = fractx[i][YY];
	    drZZ = fractx[i][ZZ];

  	    /* dr is relative offset from lower cell limit */
	    dataXX=theta[XX]+i*4;
	    dataYY=theta[YY]+i*4;
	    dataZZ=theta[ZZ]+i*4;

	    dataXX[3]=0;
	    dataYY[3]=0;
	    dataZZ[3]=0;
	    dataXX[1]=drXX;
	    dataYY[1]=drYY;
	    dataZZ[1]=drZZ;
	    dataXX[0]=1.0-drXX;
	    dataYY[0]=1.0-drYY;
	    dataZZ[0]=1.0-drZZ;

        dataXX[2]=0.5*drXX*dataXX[1];
        dataYY[2]=0.5*drYY*dataYY[1];
        dataZZ[2]=0.5*drZZ*dataZZ[1];
    
        dataXX[1]=0.5*((drXX+1.0)*dataXX[0]+(2.0-drXX)*dataXX[1]);
        dataYY[1]=0.5*((drYY+1.0)*dataYY[0]+(2.0-drYY)*dataYY[1]);
        dataZZ[1]=0.5*((drZZ+1.0)*dataZZ[0]+(2.0-drZZ)*dataZZ[1]);
        
        dataXX[0]=0.5*(1.0-drXX)*dataXX[0];
        dataYY[0]=0.5*(1.0-drYY)*dataYY[0];
        dataZZ[0]=0.5*(1.0-drZZ)*dataZZ[0];
    
	    /* differentiate */
	    ddataXX  = dtheta[XX]+i*4;
	    ddataYY  = dtheta[YY]+i*4;
	    ddataZZ  = dtheta[ZZ]+i*4;
    
	    ddataXX[0] = -dataXX[0];
	    ddataYY[0] = -dataYY[0];
	    ddataZZ[0] = -dataZZ[0];
	    ddataXX[1] = dataXX[0]-dataXX[1];
	    ddataYY[1] = dataYY[0]-dataYY[1];
	    ddataZZ[1] = dataZZ[0]-dataZZ[1];
	    ddataXX[2] = dataXX[1]-dataXX[2];
	    ddataYY[2] = dataYY[1]-dataYY[2];
	    ddataZZ[2] = dataZZ[1]-dataZZ[2];
	    ddataXX[3] = dataXX[2]-dataXX[3];
	    ddataYY[3] = dataYY[2]-dataYY[3];
	    ddataZZ[3] = dataZZ[2]-dataZZ[3];
    
	    div=1.0/3.0;
	    dataXX[3]=div*drXX*dataXX[2];
	    dataYY[3]=div*drYY*dataYY[2];
	    dataZZ[3]=div*drZZ*dataZZ[2];

        dataXX[2]=div*((drXX+1.0)*dataXX[1]+(3.0-drXX)*dataXX[2]);
        dataYY[2]=div*((drYY+1.0)*dataYY[1]+(3.0-drYY)*dataYY[2]);
        dataZZ[2]=div*((drZZ+1.0)*dataZZ[1]+(3.0-drZZ)*dataZZ[2]);
    
        dataXX[1]=div*((drXX+2.0)*dataXX[0]+(2.0-drXX)*dataXX[1]);
        dataYY[1]=div*((drYY+2.0)*dataYY[0]+(2.0-drYY)*dataYY[1]);
        dataZZ[1]=div*((drZZ+2.0)*dataZZ[0]+(2.0-drZZ)*dataZZ[1]);
    
	    dataXX[0]=div*(1.0-drXX)*dataXX[0]; 
	    dataYY[0]=div*(1.0-drYY)*dataYY[0]; 
	    dataZZ[0]=div*(1.0-drZZ)*dataZZ[0]; 
      } 
    }
  }
  else
  {
    /* general case, order != 4 */
    for(i=0; (i<nr); i++) 
	{
      if (charge[i] != 0.0) 
	  {
        xptr = fractx[i];
        for(j=0; (j<DIM); j++) 
	    {
	      dr  = xptr[j];
	
	      /* dr is relative offset from lower cell limit */
	      data=&(theta[j][i*order]);
	      data[order-1]=0;
	      data[1]=dr;
	      data[0]=1-dr;
		
	      for(k=3; (k<order); k++) 
		  {
	        div=1.0/(k-1.0);    
	        data[k-1]=div*dr*data[k-2];
	        for(l=1; (l<(k-1)); l++)
	        {
			  data[k-l-1]=div*((dr+l)*data[k-l-2]+(k-l-dr)*data[k-l-1]);
		    }
  	        data[0]=div*(1-dr)*data[0];
	      }
	      /* differentiate */
	      ddata    = &(dtheta[j][i*order]);
	      ddata[0] = -data[0];
	      for(k=1; (k<order); k++)
          {
	        ddata[k]=data[k-1]-data[k];
		  }
	      div=1.0/(order-1);
	      data[order-1]=div*dr*data[order-2];
	      for(l=1; (l<(order-1)); l++)
	      {
            data[order-l-1]=div*((dr+l)*data[order-l-2]+(order-l-dr)*data[order-l-1]);
          }
	      data[0]=div*(1-dr)*data[0]; 
        }
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

/* Global variables! Yucky... */
static    t_fftgrid *gridA=NULL,*gridB=NULL;
/*static    int  nx,ny,nz;*/
static    int  *nnx,*nny,*nnz;
static    ivec *idx=NULL;
static    rvec *fractx; /* Fractional coordinate relative to the
			 * lower cell boundary 
			 */
static    matrix    recipbox;
static    splinevec theta;
static    splinevec dtheta;
static    splinevec bsp_mod;


t_fftgrid *init_pme(FILE *log,t_commrec *cr,
		    int nkx,int nky,int nkz,int pme_order,int homenr,
		    bool bFreeEnergy,bool bOptFFT,int ewald_geometry)
{
  int i;

  fprintf(log,"Will do PME sum in reciprocal space.\n");
  please_cite(log,"Essman95a");

  if(ewald_geometry==eewg3DC) {
    fprintf(log,"Using the Ewald3DC correction for systems with a slab geometry.\n");
    please_cite(log,"In-Chul99a");
  }

  bPar = cr && (cr->nnodes>1);
  if (bPar) {
    leftnode  = cr->nodeid - 1;
    rightnode = cr->nodeid + 1;
    if(leftnode == -1)
      leftnode = cr->nnodes - 1;
    if(rightnode == cr->nnodes)
      rightnode = 0;
    leftbnd  = pme_order/2 - 1;
    rightbnd = pme_order - leftbnd - 1;

    fprintf(log,"Parallelized PME sum used.\n");
    if ((nkx % cr->nnodes) != 0)
    {
        gmx_fatal(FARGS,"fourier_nx must be divisible by NNODES\n");
    }
  } 
 
  /* allocate space for things */
  snew(bsp_mod[XX],nkx);
  snew(bsp_mod[YY],nky);
  snew(bsp_mod[ZZ],nkz);
  for(i=0;i<DIM;i++) 
  {
    snew(theta[i],pme_order*homenr); 
    snew(dtheta[i],pme_order*homenr);
  }
  snew(fractx,homenr); 

  snew(idx,homenr);
  snew(nnx,5*nkx);
  snew(nny,5*nky);
  snew(nnz,5*nkz);
  for(i=0; (i<5*nkx); i++)
    nnx[i] = i % nkx;
  for(i=0; (i<5*nky); i++)
    nny[i] = i % nky;
  for(i=0; (i<5*nkz); i++)
    nnz[i] = i % nkz;

  gridA = mk_fftgrid(log,nkx,nky,nkz,cr);
  if (bFreeEnergy)
  {
    gridB = mk_fftgrid(log,nkx,nky,nkz,cr);
  }
  
  make_bspline_moduli(bsp_mod,nkx,nky,nkz,pme_order);   

  return gridA;
}

void spread_on_grid(FILE *logfile,   
		    t_fftgrid *grid,  int homenr,
		    int pme_order,    rvec x[],
		    real charge[],    matrix box,
		    bool bGatherOnly, bool bHaveSplines)
{ 
  int nx,ny,nz,nx2,ny2,nz2,la2,la12;
  real *ptr;
  
  /* Unpack structure */
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  
  if (!bHaveSplines)
    /* Inverse box */
    calc_recipbox(box,recipbox); 
  
  if (!bGatherOnly) {
    if (!bHaveSplines) {
      /* Compute fftgrid index for all atoms, with help of some extra variables */
      calc_idx(homenr,recipbox,x,fractx,idx,nx,ny,nz,nx2,ny2,nz2,nnx,nny,nnz);
      
      /* make local bsplines  */
      make_bsplines(theta,dtheta,pme_order,nx,ny,nz,fractx,idx,charge,homenr);
    }    

    /* put local atoms on grid. */
    spread_q_bsplines(grid,idx,charge,theta,homenr,pme_order,nnx,nny,nnz);
/*    pr_grid_dist(logfile,"spread",grid); */
  }
}

real do_pme(FILE *logfile,   bool bVerbose,
	    t_inputrec *ir,  rvec x[],
	    rvec f[],
	    real *chargeA,   real *chargeB,
	    matrix box,	     t_commrec *cr,
	    t_nsborder *nsb, t_nrnb *nrnb,    
	    matrix vir,      real ewaldcoeff,
	    bool bFreeEnergy,
	    real lambda,     real *dvdlambda,
	    bool bGatherOnly)
{ 
  static  real energy_AB[2] = {0,0};
  int     q,i,j,ntot,npme;
  int     nx,ny,nz,nx2,ny2,nz2,la12,la2;
  t_fftgrid *grid=NULL;
  real *ptr;
  real    *homecharge=NULL,vol,energy;
  matrix  vir_AB[2];
  static bool bFirst=TRUE;
  static int *pidx;
  static int *gidx=NULL; /* Grid index of particle, used for PME redist */
  static rvec *x_tmp;
  static rvec *f_tmp;
  static real *q_tmp;
  static FILE *fp_f;
  static char fn[8];

  if(bFirst) {
#if defined GMX_MPI
    if (bPar) {
      snew(gidx,nsb->natoms);
      snew(x_tmp,nsb->natoms);
      snew(f_tmp,nsb->natoms);
      snew(q_tmp,nsb->natoms);
      snew(pidx,nsb->natoms);
      MPI_Type_contiguous(DIM, mpi_type, &rvec_mpi);
      MPI_Type_commit(&rvec_mpi);
    }
#endif
    /*#define PRT_FORCE */
#ifdef PRT_FORCE
    sprintf(fn,"force%2.2d",cr->nodeid);
    fp_f=(FILE *)fopen(fn,"w");
#endif
    where();
    bFirst=FALSE;
  }
  for(q=0; q<(bFreeEnergy ? 2 : 1); q++) {
    if (q == 0) {
      grid = gridA;
      homecharge = chargeA+START(nsb);
    } else {
      grid = gridB;
      homecharge = chargeB+START(nsb);
    }
    /* Unpack structure */
    if (debug) {
      fprintf(debug,"PME: nnodes = %d, nodeid = %d\n",cr->nnodes,cr->nodeid);
      fprintf(debug,"Grid = %p\n",grid);
      if (grid == NULL)
	gmx_fatal(FARGS,"No grid!");
    }
    where();
    unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
    where();

    my_homenr=nsb->natoms;

    if (!bPar) {
      x_tmp=x;
      f_tmp=f;
      q_tmp=homecharge;
    } else {
      pme_calc_pidx(cr,HOMENR(nsb),box,x+START(nsb),ir->nkx,nnx,pidx,gidx);
      where();

      pmeredist(cr, TRUE, HOMENR(nsb), x+START(nsb), homecharge, gidx, 
		&my_homenr, x_tmp, q_tmp);
      where();
    }
    where();
    if (debug)
      fprintf(debug,"Node= %6d, pme local particles=%6d\n",
	      cr->nodeid,my_homenr);

    /* Spread the charges on a grid */
#ifdef USE_MPE
    MPE_Log_event( ev_spread_on_grid_start, 0, "" );
#endif

    /* Spread the charges on a grid */
    spread_on_grid(logfile,grid,my_homenr,ir->pme_order,
		   x_tmp,q_tmp,box,bGatherOnly,
		   q==0 ? FALSE : TRUE);
#ifdef USE_MPE
    MPE_Log_event( ev_spread_on_grid_finish, 0, "");
#endif
    if (!bGatherOnly) {
      inc_nrnb(nrnb,eNR_SPREADQBSP,
	       ir->pme_order*ir->pme_order*ir->pme_order*my_homenr);
      
      /* sum contributions to local grid from other nodes */
      if (bPar) {

	sum_qgrid(cr,nsb,grid,ir->pme_order,TRUE);
	where();
      }
#ifdef DEBUG
      if (debug)
	pr_fftgrid(debug,"qgrid",grid);
#endif
      where();

      /* do 3d-fft */ 
#ifdef USE_MPE
      MPE_Log_event( ev_gmxfft3d_start, 0, "" );
#endif
      gmxfft3D(grid,GMX_FFT_REAL_TO_COMPLEX,cr);
#ifdef USE_MPE
      MPE_Log_event( ev_gmxfft3d_finish, 0, "" );
#endif
      where();

      /* solve in k-space for our local cells */
      vol = det(box);
#ifdef USE_MPE
      MPE_Log_event( ev_solve_pme_start, 0, "" );
#endif
      energy_AB[q]=solve_pme(grid,ewaldcoeff,vol,ir->epsilon_r,
			     bsp_mod,recipbox,vir_AB[q],cr);
      where();
#ifdef USE_MPE
      MPE_Log_event( ev_solve_pme_finish, 0, "" );
#endif
      inc_nrnb(nrnb,eNR_SOLVEPME,nx*ny*nz*0.5);

      /* do 3d-invfft */
#ifdef USE_MPE
      MPE_Log_event( ev_gmxfft3d_start, 0, "" );
#endif
      where();
      gmxfft3D(grid,GMX_FFT_COMPLEX_TO_REAL,cr);
      where();
#ifdef USE_MPE
      MPE_Log_event( ev_gmxfft3d_finish, 0, "" );
#endif
      
      /* distribute local grid to all nodes */
      if (bPar) {
	sum_qgrid(cr,nsb,grid,ir->pme_order,FALSE);
      }
      where();
#ifdef DEBUG
      if (debug)
	pr_fftgrid(debug,"potential",grid);
#endif

      ntot  = grid->nxyz;  
      npme  = ntot*log((real)ntot)/log(2.0);
      if (bPar)
	npme /= cr->nnodes;
      inc_nrnb(nrnb,eNR_FFT,2*npme);
      where();
    }
    /* interpolate forces for our local atoms */
#ifdef USE_MPE
    MPE_Log_event( ev_gather_f_bsplines_start, 0, "" );
#endif
    if (bPar) {
      for(i=0; (i<my_homenr); i++) {
	f_tmp[i][XX]=0;
	f_tmp[i][YY]=0;
	f_tmp[i][ZZ]=0;
      }
    }
    where();
    gather_f_bsplines(grid,recipbox,idx,f_tmp,q_tmp,
		      bFreeEnergy ? (q==0 ? 1.0-lambda : lambda) : 1.0,
		      theta,dtheta,my_homenr,ir->pme_order,
		      nnx,nny,nnz);
    where();
    
#ifdef USE_MPE
    MPE_Log_event( ev_gather_f_bsplines_finish, 0, "" );
#endif
    
    if (bPar) {
      pmeredist(cr, FALSE,HOMENR(nsb), f+START(nsb), homecharge, gidx, 
		&my_homenr, f_tmp, q_tmp);
    }
    where();
    
#ifdef PRT_FORCE
    for(i=START(nsb); (i<START(nsb)+HOMENR(nsb)); i++) { 
      fprintf(fp_f,"force %5d  %12.4e %12.4e %12.4e %12.4e\n",
              i,f[i][XX],f[i][YY],f[i][ZZ],homecharge[i-START(nsb)]);
    }
#endif
    
    inc_nrnb(nrnb,eNR_GATHERFBSP,
	     ir->pme_order*ir->pme_order*ir->pme_order*my_homenr);
  }
  
  if (!bFreeEnergy) {
    energy = energy_AB[0];
    copy_mat(vir_AB[0],vir);
  } else {
    energy = (1.0-lambda)*energy_AB[0] + lambda*energy_AB[1];
    *dvdlambda += energy_AB[1] - energy_AB[0];
    for(i=0; i<DIM; i++)
      for(j=0; j<DIM; j++)
	vir[i][j] = (1.0-lambda)*vir_AB[0][i][j] + lambda*vir_AB[1][i][j];
  }
  /* fprintf(logfile,"PME mesh energy: %g\n",energy);*/
  if (debug)
    fprintf(debug,"PME mesh energy: %g\n",energy);

  return energy;
}
