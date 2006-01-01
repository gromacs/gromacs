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
#include "fatal.h"
#include "ewald.h"
#include "pme.h"
#include "network.h"
#include "physics.h"
#include "nrnb.h"
#include "copyrite.h"
#ifdef GMX_MPI
#include <mpi.h>
#endif
#include "mpelogging.h"

#define DFT_TOL 1e-7
/* #define PRT_FORCE */
/* conditions for on the fly time-measurement */
#define TAKETIME (step > 1 && timesteps < 10)

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
static int leftbnd, rightbnd;   /* size of left and right boundary */
static int my_homenr;

/* The following stuff is needed for signal handling on the PME nodes. 
 * signal_handler needs to be defined in md.c, the bGot..Signal variables
 * here */ 
extern RETSIGTYPE signal_handler(int n);

volatile bool bGotTermSignal = FALSE, bGotUsr1Signal = FALSE; 

/* #define SORTPME */


void pr_grid_dist(FILE *fp,char *title,t_fftgrid *grid)
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

static int comp_pmeatom(const void *a,const void *b)
{
  t_pmeatom *pa,*pb;
  pa=(t_pmeatom *)a;
  pb=(t_pmeatom *)b;
  
  return pa->index-pb->index;
}

void calc_idx(t_fftgrid *grid,int natoms,matrix recipbox,
	      rvec x[],rvec fractx[],ivec idx[],
	      //int nx,int ny,int nz,
	      //int nx2,int ny2,int nz2,
	      int nnx[],int nny[],int nnz[],
	      t_pmeatom pmeatom[])
{
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
 
  /* Unpack structure */
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  
  rxx = recipbox[XX][XX];
  ryx = recipbox[YY][XX];
  ryy = recipbox[YY][YY];
  rzx = recipbox[ZZ][XX];
  rzy = recipbox[ZZ][YY];
  rzz = recipbox[ZZ][ZZ];

  for(i=0; (i<natoms); i++) {
    xptr   = x[i];
    idxptr = idx[i];
    fptr   = fractx[i];
    
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
    
    idxptr[XX] = nnx[tix];
    idxptr[YY] = nny[tiy];
    idxptr[ZZ] = nnz[tiz];

#ifdef SORTPME    
    pmeatom[i].atom  = i;
    pmeatom[i].index = INDEX(idxptr[XX],idxptr[YY],idxptr[ZZ]);
#endif
    
#ifdef DEBUG
    range_check(idxptr[XX],0,nx);
    range_check(idxptr[YY],0,ny);
    range_check(idxptr[ZZ],0,nz);
#endif
  }  
#if (defined __GNUC__ && (defined i386 || defined __386__) && !defined GMX_DOUBLE && defined GMX_X86TRUNC)  
  asm("fldcw %0" : : "m" (*&x86_cwsave));
#endif

#ifdef SORTPME
  GMX_MPE_LOG(ev_sort_start);
  qsort(pmeatom,natoms,sizeof(pmeatom[0]),comp_pmeatom);
  GMX_MPE_LOG(ev_sort_finish);
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
  int  npme; /* how many nodes participate in PME? */

  if (cr->npmenodes)
    npme = cr->npmenodes;
  else 
    npme = cr->nnodes;
    
/* Calculate PME task index (pidx) for each grid index */
  irange=1+(nx-1)/npme;

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
  static int  npme; /* how many nodes participate in PME? */

  if(bFirst) {
    if (cr->npmenodes)
      npme = cr->npmenodes;
    else 
      npme = cr->nnodes;

    snew(scounts,npme);
    snew(rcounts,npme);
    snew(sdispls,npme);
    snew(rdispls,npme);
    snew(sidx,npme);
    snew(buf,n*DIM);
    bFirst=FALSE;
  }

#ifdef GMX_MPI
  if (forw) { /* forward, redistribution from pp to pme */ 

/* Calculate send counts and exchange them with other nodes */
    for(i=0; (i<npme); i++) scounts[i]=0;
    for(i=0; (i<n); i++) scounts[idxa[i]]++;
    MPI_Alltoall( scounts, 1, MPI_INT, rcounts, 1, MPI_INT, cr->mpi_comm_mygroup );

/* Calculate send and receive displacements and index into send buffer */
    sdispls[0]=0;
    rdispls[0]=0;
    sidx[0]=0;
    for(i=1; i<npme; i++) {
      sdispls[i]=sdispls[i-1]+scounts[i-1];
      rdispls[i]=rdispls[i-1]+rcounts[i-1];
      sidx[i]=sdispls[i];
    }
/* Total # of particles to be received */
    *total=rdispls[npme-1]+rcounts[npme-1];

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
                  cr->mpi_comm_mygroup);

/* Copy charge into send buffer and exchange*/
    for(i=0; i<npme; i++) sidx[i]=sdispls[i];
    for(i=0; (i<n); i++) {
      ii=sidx[idxa[i]];
      sidx[idxa[i]]++;
      buf[ii]=charge[i];
    }
    MPI_Alltoallv(buf, scounts, sdispls, mpi_type,
                  pme_q, rcounts, rdispls, mpi_type,
                  cr->mpi_comm_mygroup);

  }
  else { /* backward, redistribution from pme to pp */ 
    MPI_Alltoallv(pme_x, rcounts, rdispls, rvec_mpi,
                  buf, scounts, sdispls, rvec_mpi,
                  cr->mpi_comm_mygroup);

/* Copy data from receive buffer */
    for(i=0; i<npme; i++) sidx[i]=sdispls[i];
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

void sum_qgrid_dd(t_commrec *cr,t_nsborder *nsb,t_fftgrid *grid,
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
  static int nodeshift=0;

  GMX_MPE_LOG(ev_sum_qgrid_start);
  
#ifdef GMX_MPI
  if(bFirst) {
    bFirst=FALSE;
    localsize=grid->la12r*grid->pfft.local_nx;
    if(!grid->workspace)
      snew(tmp,localsize);
    maxproc=grid->nx/grid->pfft.local_nx;
    /* Nodeshift: same meaning as in spread_q_bsplines */
    if (pmeduty(cr)==epmePMEONLY)
      nodeshift=cr->nnodes-cr->npmenodes;
  }
  /* NOTE: FFTW doesnt necessarily use all processors for the fft;
     * above I assume that the ones that do have equal amounts of data.
     * this is bad since its not guaranteed by fftw, but works for now...
     * This will be fixed in the next release.
     */
  if(grid->workspace)
    tmp=grid->workspace;
  if(bForward) { /* sum contributions to local grid */
    ny=grid->ny;
    la2r=grid->la2r;
/* Send left Boundary */
    bndsize = leftbnd*ny*la2r;
    from = grid->ptr + (cr->left   - nodeshift + 1)*localsize - bndsize;
    to   = grid->ptr + (cr->nodeid - nodeshift + 1)*localsize - bndsize;
    MPI_Sendrecv(from,bndsize,mpi_type,
               cr->left,cr->nodeid,
               tmp,bndsize,mpi_type,cr->right,cr->right,
               MPI_COMM_WORLD,&stat);
    GMX_MPE_LOG(ev_test_start); 
    for(i=0; (i<bndsize); i++) {
      to[i] += tmp[i];
    }
    GMX_MPE_LOG(ev_test_finish);

/* Send right boundary */
    bndsize = rightbnd*ny*la2r;
    from = grid->ptr + (cr->right  - nodeshift)*localsize;
    to   = grid->ptr + (cr->nodeid - nodeshift)*localsize;
    MPI_Sendrecv(from,bndsize,mpi_type,
               cr->right,cr->nodeid,
               tmp,bndsize,mpi_type,cr->left,cr->left,
               MPI_COMM_WORLD,&stat); 
    GMX_MPE_LOG(ev_test_start);
    for(i=0; (i<bndsize); i++) {
      to[i] += tmp[i];
    }
    GMX_MPE_LOG(ev_test_finish);
  }
  else { /* distribute local grid to all processors */
    ny=grid->ny;
    la2r=grid->la2r;

/* Send left Boundary */
    bndsize = rightbnd*ny*la2r;
    from = grid->ptr + (cr->nodeid - nodeshift)*localsize;
    to   = grid->ptr + (cr->right  - nodeshift)*localsize;
    MPI_Sendrecv(from,bndsize,mpi_type,
               cr->left,cr->nodeid,
               to,bndsize,mpi_type,cr->right,cr->right,
               MPI_COMM_WORLD,&stat);
/* Send right boundary */
    bndsize = leftbnd*ny*la2r;
    from = grid->ptr + (cr->nodeid - nodeshift + 1)*localsize - bndsize;
    to   = grid->ptr + (cr->left   - nodeshift + 1)*localsize - bndsize;
    MPI_Sendrecv(from,bndsize,mpi_type,
               cr->right,cr->nodeid,
               to,bndsize,mpi_type,cr->left,cr->left,
               MPI_COMM_WORLD,&stat);
  }
#else
  gmx_fatal(FARGS,"Parallel grid summation requires MPI and FFTW.\n");    
#endif
  GMX_MPE_LOG(ev_sum_qgrid_finish);
}

void sum_qgrid(t_commrec *cr,t_nsborder *nsb,t_fftgrid *grid,
               int pme_order, bool bForward)
{
  static bool bFirst=TRUE;
  static real *tmp;
  int i;
  static int localsize;
  static int maxproc;

#ifdef GMX_MPI
  if(bFirst) {
    localsize=grid->la12r*grid->pfft.local_nx;
    if(!grid->workspace)
      snew(tmp,localsize);
    maxproc=grid->nx/grid->pfft.local_nx;
  }
  /* NOTE: FFTW doesnt necessarily use all processors for the fft;
     * above I assume that the ones that do have equal amounts of data.
     * this is bad since its not guaranteed by fftw, but works for now...
     * This will be fixed in the next release.
     */
  bFirst=FALSE;
  if(grid->workspace)
    tmp=grid->workspace;
  if(bForward) { /* sum contributions to local grid */

    GMX_BARRIER(cr->mpi_comm_mygroup);
    GMX_MPE_LOG(ev_reduce_start);
    for(i=0;i<maxproc;i++) {
      MPI_Reduce(grid->ptr+i*localsize, /*ptr arithm.     */
		 tmp,localsize,      
		 GMX_MPI_REAL,MPI_SUM,i,cr->mpi_comm_mygroup);
    }
    GMX_MPE_LOG(ev_reduce_finish);

    if(cr->nodeid<maxproc)
      memcpy(grid->ptr+cr->nodeid*localsize,tmp,localsize*sizeof(real));
  }
  else { /* distribute local grid to all processors */
    for(i=0;i<maxproc;i++)
      MPI_Bcast(grid->ptr+i*localsize, /* ptr arithm     */
		localsize,       
		GMX_MPI_REAL,i,cr->mpi_comm_mygroup);
  }
#else
  gmx_fatal(FARGS,"Parallel grid summation requires MPI and FFTW.\n");    
#endif
}

void spread_q_bsplines(t_fftgrid *grid,ivec idx[],real charge[],
		       splinevec theta,int nr,int order,
		       int nnx[],int nny[],int nnz[], t_commrec *cr,
		       t_pmeatom pmeatom[])
{
  /* spread charges from home atoms to local grid */
  real     *ptr;
  int      i,nn,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
  int      nx,ny,nz,nx2,ny2,nz2,la2,la12;
  int      norder,*idxptr,index_x,index_xy,index_xyz;
  real     valx,valxy,qn;
  real     *thx,*thy,*thz;
  int localsize, bndsize;
  static bool bFirst=TRUE;
  static int nodeshift=0;

  if(bFirst) {
    bFirst=FALSE;
    /* Variable nodeshift accounts for correct calculation of ptr
       if the first PME node has a nodeid of >0. */
    if (pmeduty(cr)==epmePMEONLY)
      nodeshift = cr->nnodes-cr->npmenodes;
  }
  
  if (!bPar) {
    clear_fftgrid(grid); 
#ifdef GMX_MPI
  } else {
    localsize = grid->la12r*grid->pfft.local_nx;
    ptr = grid->localptr;
    for (i=0; (i<localsize); i++)
      ptr[i] = 0;
/* clear left boundary area */
    bndsize = leftbnd*grid->la12r;
    ptr = grid->ptr + (cr->left - nodeshift + 1)*localsize - bndsize;
    for (i=0; (i<bndsize); i++)
      ptr[i] = 0;
/* clear right boundary area */
    bndsize = rightbnd*grid->la12r;
    ptr = grid->ptr + (cr->right - nodeshift)*localsize;
    for (i=0; (i<bndsize); i++)
      ptr[i] = 0;
#endif
  }
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  ii0   = nnx+nx2+1-order/2;
  jj0   = nny+ny2+1-order/2;
  kk0   = nnz+nz2+1-order/2;
  
  for(nn=0; (nn<nr); nn++) {
#ifdef SORTPME
    n = pmeatom[nn].atom;
#else
    n = nn;
#endif
    qn     = charge[n];
    idxptr = idx[n];
    
    if (qn != 0) {
      norder  = n*order;
      
      /* Pointer arithmetic alert, next six statements */
      i0  = ii0 + idxptr[XX]; 
      j0  = jj0 + idxptr[YY];
      k0  = kk0 + idxptr[ZZ];
      thx = theta[XX] + norder;
      thy = theta[YY] + norder;
      thz = theta[ZZ] + norder;
      
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
		       int nr,int order,int nnx[],int nny[],int nnz[],
		       t_pmeatom pmeatom[])
{
  /* sum forces for local particles */  
  int     nn,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
  int     nx,ny,nz,nx2,ny2,nz2,la2,la12,index_x,index_xy;
  real *  ptr;
  real    tx,ty,dx,dy,qn;
  real    fx,fy,fz,fnx,fny,fnz,gval,tgz,dgz;
  real    gval1,gval2,gval3,gval4;
  real    fxy1,fz1;
  real    *thx,*thy,*thz,*dthx,*dthy,*dthz;
  int     sn,norder,*idxptr,ind0;
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

  for(nn=0; (nn<nr); nn++) {
#ifdef SORTPME
    n = pmeatom[nn].atom;
#else
    n = nn;
#endif
    qn      = scale*charge[n];

    if (qn != 0) {
      fx     = 0.0;
      fy     = 0.0;
      fz     = 0.0;
      fnx    = f[n][XX];
      fny    = f[n][YY];
      fnz    = f[n][ZZ];
      idxptr = idx[n];
      norder = n*order;
      
      /* Pointer arithmetic alert, next nine statements */
      i0   = ii0 + idxptr[XX]; 
      j0   = jj0 + idxptr[YY];
      k0   = kk0 + idxptr[ZZ];
      thx  = theta[XX] + norder;
      thy  = theta[YY] + norder;
      thz  = theta[ZZ] + norder;
      dthx = dtheta[XX] + norder;
      dthy = dtheta[YY] + norder;
      dthz = dtheta[ZZ] + norder;
      
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
      f[n][XX] = fnx - qn*( fx*nx*rxx );
      f[n][YY] = fny - qn*( fx*nx*ryx + fy*ny*ryy );
      f[n][ZZ] = fnz - qn*( fx*nx*rzx + fy*ny*rzy + fz*nz*rzz );
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
  real dr,div,rcons;
  real *data,*ddata,*xptr;

  for(i=0; (i<nr); i++) {
    if (charge[i] != 0.0) {
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
	  for(l=1; (l<(k-1)); l++)
	    data[k-l-1]=div*((dr+l)*data[k-l-2]+(k-l-dr)*
			     data[k-l-1]);
	  data[0]=div*(1-dr)*data[0];
	}
	/* differentiate */
	ddata    = &(dtheta[j][i*order]);
	ddata[0] = -data[0];
	for(k=1; (k<order); k++)
	  ddata[k]=data[k-1]-data[k];
		
	div=1.0/(order-1);
	data[order-1]=div*dr*data[order-2];
	for(l=1; (l<(order-1)); l++)
	  data[order-l-1]=div*((dr+l)*data[order-l-2]+
			       (order-l-dr)*data[order-l-1]);
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

/* Global variables! Yucky... */
static    t_fftgrid *gridA=NULL,*gridB=NULL;
static    int  *nnx,*nny,*nnz;
static    ivec *idx=NULL;
static    rvec *fractx; /* Fractional coordinate relative to the
			 * lower cell boundary 
			 */
static    matrix    recipbox;
static    splinevec theta;
static    splinevec dtheta;
static    splinevec bsp_mod;
static    t_pmeatom *pmeatom;

t_fftgrid *init_pme(FILE *log,t_commrec *cr,
		    int nkx,int nky,int nkz,int pme_order,int homenr,
		    bool bFreeEnergy,bool bOptFFT,int ewald_geometry)
{
  int i,npme; /* how many nodes participate in PME? */


  if (cr->npmenodes)
    npme = cr->npmenodes;
  else
    npme = cr->nnodes;
 
  fprintf(log,"Will do PME sum in reciprocal space.\n");
  please_cite(log,"Essman95a");

  if(ewald_geometry==eewg3DC) {
    fprintf(log,"Using the Ewald3DC correction for systems with a slab geometry.\n");
    please_cite(log,"In-Chul99a");
  }

  bPar = cr && (npme>1); /* bPar is a global variable! */
  if (bPar) {
    leftbnd  = pme_order/2 - 1;
    rightbnd = pme_order - leftbnd - 1;

    fprintf(log,"Parallelized PME sum used. nkx=%d, npme=%d\n",nkx,npme);
    if ((nkx % npme) != 0)
      fprintf(log,"Warning: For load balance, "
	      "fourier_nx should be divisible by the number of PME nodes\n");
  } 
 
  /* allocate space for things */
  snew(bsp_mod[XX],nkx);
  snew(bsp_mod[YY],nky);
  snew(bsp_mod[ZZ],nkz);
  for(i=0;i<DIM;i++) {
    snew(theta[i],pme_order*homenr); 
    snew(dtheta[i],pme_order*homenr);
  }
  snew(fractx,homenr); 
  snew(pmeatom,homenr);
  
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
    gridB = mk_fftgrid(log,nkx,nky,nkz,cr);

  make_bspline_moduli(bsp_mod,nkx,nky,nkz,pme_order);   

  return gridA;
}

void spread_on_grid(FILE *logfile,    t_commrec *cr, 
		    t_fftgrid *grid,  int homenr,
		    int pme_order,    rvec x[],
		    real charge[],    matrix box,
		    bool bGatherOnly, bool bHaveSplines)
/*		    t_pmeatom pmeatom[]) */
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
//      calc_idx(grid,homenr,recipbox,x,fractx,idx,nx,ny,nz,nx2,ny2,nz2,
      calc_idx(grid,homenr,recipbox,x,fractx,idx,
	       nnx,nny,nnz,pmeatom);
      
      /* make local bsplines  */
      make_bsplines(theta,dtheta,pme_order,nx,ny,nz,fractx,idx,charge,homenr);
    }    

    /* put local atoms on grid. */
    spread_q_bsplines(grid,idx,charge,theta,homenr,pme_order,nnx,nny,nnz,cr,
		      pmeatom);
/*    pr_grid_dist(logfile,"spread",grid); */
  }
}

/*************/
#ifdef GMX_MPI  
/*************/
/* The following three routines are for split-off pme nodes */
void send_coordinates(t_nsborder *nsb, rvec x[], 
                      real chargeA[],  real chargeB[], 
                      bool bLastTime,  bool bFreeEnergy, 
		      t_commrec *cr)
{
  int receiver, i;
  int firstpmenode=nsb->nnodes-nsb->npmenodes;
  int lastpmenode =nsb->nnodes-1;
  int firstindex,lastindex;
  static bool bFirst=TRUE;
  static int nsend;
  static MPI_Request req[2];
  static MPI_Status  stat[2];
  FILE *fp;
  char fn[8];

#ifdef PRT_COORDS
  if (bFirst) {
    sprintf(fn,"xsend%2.2d",cr->nodeid);
    fp=fopen(fn,"w");
    for (i=START(nsb); (i<START(nsb)+HOMENR(nsb)); i++) { 
      fprintf(fp,"x %5d  %20.11e %20.11e %20.11e\n",
              i,x[i][XX],x[i][YY],x[i][ZZ]);
      }
    fclose(fp);
  }
#endif


  /* Communicate coordinates and charges to the PME nodes */
  /* find out which of our local data has to be sent do pme node i!
   * local data has indices from START(nsb) to START(nsb)+HOMENR(nsb)-1 
   */
  nsend=0;
  for (receiver=firstpmenode; receiver<=lastpmenode; receiver++) 
  {
    firstindex=max(START(nsb),nsb->pmeindex[receiver]);
    lastindex =min(START(nsb)+HOMENR(nsb)-1,nsb->pmeindex[receiver]+nsb->pmehomenr[receiver]-1);
    if (lastindex >= firstindex)
    {
/*    fprintf(stderr,"SEND: Node %d -> node %d (Indices %d to %d) - %f %f %f\n",
 *            nsb->nodeid,receiver,firstindex,lastindex,x[firstindex][0],x[firstindex][1],x[firstindex][2]);
 */
      if (MPI_Isend(&x[firstindex][0],(lastindex-firstindex+1)*DIM,
                    GMX_MPI_REAL,receiver,bLastTime,MPI_COMM_WORLD,&req[nsend]) !=0 )
        gmx_comm("MPI_Isend failed in send_coordinates");
      nsend++;
      
      if (bFirst) /* we need to send the charges in the first step 
                   * (this can be done using blocking sends) */
      {
        MPI_Send(&chargeA[firstindex],(lastindex-firstindex+1),
	         GMX_MPI_REAL,receiver,bLastTime,MPI_COMM_WORLD);
        if (bFreeEnergy) 
          MPI_Send(&chargeB[firstindex],(lastindex-firstindex+1),
	           GMX_MPI_REAL,receiver,bLastTime,MPI_COMM_WORLD);
      }
    }
  }
  bFirst=FALSE;
  /* wait for the communication to finish! */
  MPI_Waitall(nsend, req, stat);
}


void receive_lrforces(t_commrec *cr, t_nsborder *nsb, 
                      rvec f[]     , matrix vir, 
		      real *energy, int step)
{
  int  i,j,sender,elements;
  static int  firstpmenode, lastpmenode;
  int  firstindex, lastindex;
  real  vbuf[DIM*DIM];
  MPI_Status  status;
  static rvec  *f_rec;  /* force receive buffer, don't lose address across calls! */
  static bool  bFirst=TRUE;
  MPI_Request req[2];
  MPI_Status  stat[2];
  int nrecv;
  double t0, t1;
  static double t2=0, t3;
  static double t_receive_f=0.0, tstep_sum=0.0;
  static int timesteps=0;

  if (bFirst)
  {
    firstpmenode=nsb->nnodes-nsb->npmenodes;
    lastpmenode =nsb->nnodes-1;
    snew(f_rec,HOMENR(nsb));
    bFirst=FALSE;
  }

  *energy=0.0;
  
#ifdef PRT_FORCELR
    static FILE *fp_f1;
    static char fn1[10];

    /* long-range forces from the PME nodes */
    sprintf(fn1,"f_recv%2.2d",cr->nodeid);
    fp_f1=(FILE *)fopen(fn1,"w");
#endif

  /* this is to help the user balance the ratio
   * of PME versus PP nodes */
  t3=MPI_Wtime()-t2;
  t2=MPI_Wtime();
  if (TAKETIME)
  {
    t0=MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    t1=MPI_Wtime();
    t_receive_f += t1-t0;
    tstep_sum   += t3;
    timesteps++;
/*    fprintf(stderr," PP node %d, this time step %f, sum %f, steps %d, length %f, sum %f, ratio %f\n",
 *            cr->nodeid, t1-t0, t_receive_f, timesteps, t3, tstep_sum, t_receive_f/tstep_sum);
 */
    if (timesteps == 10)
      fprintf(stderr, "PP node %d waits %3.0f%% of the time in do_force.\n", cr->nodeid, 100*t_receive_f/tstep_sum);
  }

  /* now receive the forces */
  nrecv=0;
  for (sender=firstpmenode; sender<=lastpmenode; sender++)
  {
    firstindex=max(START(nsb), nsb->pmeindex[sender]);
    lastindex =min(START(nsb)+HOMENR(nsb)-1, nsb->pmeindex[sender]+nsb->pmehomenr[sender]-1);
    if (lastindex >= firstindex)
    {
      if (MPI_Irecv(f_rec[firstindex-START(nsb)],(lastindex-firstindex+1)*DIM,
                    GMX_MPI_REAL,sender,0,MPI_COMM_WORLD,&req[nrecv]) != 0)
        gmx_comm("MPI_Irecv failed in receive_lrforces");
      nrecv++;
/*      fprintf(stderr,"   RECEIVE F: Node %d from node %d (Indices %d to %d) - %12.4e %12.4e %12.4e\n",
 *              nsb->nodeid,sender,firstindex,lastindex,
 *  	      f_rec[firstindex-START(nsb)][0],f_rec[firstindex-START(nsb)][1],f_rec[firstindex-START(nsb)][2]);
 */
    }
  }
  MPI_Waitall(nrecv, req, stat);

  /* add PME forces to f (that already contains corrections 
   * and must not be overridden at this point!) */
  for (i=START(nsb); i<START(nsb)+HOMENR(nsb); i++)
    rvec_inc(f[i], f_rec[i-START(nsb)]);
  
  /* receive energy, virial */
  sender=cr->nodeid+(cr->nnodes-cr->npmenodes);
  if (sender < cr->nnodes) 
  {
    /* receive energy 
     * (The tag of the energy message is used to transport a
     * TERM/USR1 signal if necessary)
     */
    MPI_Recv(energy,1,GMX_MPI_REAL,sender,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    if (status.MPI_TAG == 1)
      bGotTermSignal=TRUE;
    else if (status.MPI_TAG == 2)
      bGotUsr1Signal=TRUE;
	                
    /* receive virial buffer */
    MPI_Recv(&vbuf,DIM*DIM,GMX_MPI_REAL,sender,2,MPI_COMM_WORLD,&status);

    /* add receive buffer to vir variable */
    elements=0;
    for (i=0; i<DIM; i++)
      for (j=0; j<DIM; j++)
        {
	  vir[i][j] += vbuf[elements];
	  elements++;
	}
  }
    
#ifdef PRT_FORCELR
  for(i=0; i<HOMENR(nsb); i++) 
  { 
    fprintf(fp_f1,"force %5d  %20.12e %20.12e %20.12e\n",
            i+START(nsb),f_rec[i][XX],f_rec[i][YY],f_rec[i][ZZ]);

  }
#endif

}


void do_pmeonly(FILE *logfile, 
                t_inputrec *ir,  t_commrec *cr,
		matrix box, 
		t_nsborder *nsb, t_nrnb *mynrnb,
		                 real ewaldcoeff,
		real lambda,
		bool bGatherOnly)
		
{
  MPI_Status  status;
  rvec  *xpme,*fpme;
  real  *chargeA,*chargeB=NULL,*chargepme;
  int  sender,receiver,elements,tag;
  int  firstindex, lastindex;
  int  nx,ny,nz,nx2,ny2,nz2,la12,la2;
  int  i,j,q,ntot,npme;
  real *ptr;
  t_fftgrid  *grid;
  int  *pidx=NULL;
  int  *gidx=NULL; /* Grid index of particle, used for PME redist */
  rvec  *x_tmp=NULL;
  rvec  *f_tmp;
  real  *q_tmp=NULL;
  FILE  *fp_f3,*fp_f4;
  char  fn[8];
  real  energy_AB[2] = {0,0};
  matrix  vir,vir_AB[2];
  real  vbuf[DIM*DIM];
  bool  bFreeEnergy, bFirst=TRUE;
  real  vol,energy;
  real  *dvdlambda=0;
  int messages; /* count the Isends or Ireceives */    
  MPI_Request *req;
  MPI_Status  *stat;
  bool bDone=FALSE;
  int count=0;
  double t0, t1, t2=0, t3, t_send_f=0.0, tstep_sum=0.0;
  bool bNS;
  unsigned long int step=0, timesteps=0;

  bFreeEnergy = (ir->efep != efepNO);
 
#define ME (nsb->nodeid) 

  /* Allocate space for the request handles, a lonely PME node may
     receive (and later send back) as many messages as there are PP nodes */
  snew(req ,nsb->nnodes-nsb->npmenodes);
  snew(stat,nsb->nnodes-nsb->npmenodes);
          
  snew(xpme,PMEHOMENR(nsb));
  snew(fpme,PMEHOMENR(nsb));
  snew(chargeA,PMEHOMENR(nsb));
  if (bFreeEnergy)
    snew(chargeB,PMEHOMENR(nsb));
  snew(f_tmp,nsb->natoms);
  if (cr->npmenodes > 1) {
    snew(x_tmp,nsb->natoms);
    snew(q_tmp,nsb->natoms);
    snew(pidx,nsb->natoms);
    snew(gidx,nsb->natoms);
  }
  MPI_Type_contiguous(DIM, mpi_type, &rvec_mpi);
  MPI_Type_commit(&rvec_mpi);

  
#ifdef PRT_FORCE
  sprintf(fn,"f_pme%2.2d",cr->nodeid);
  fp_f3=(FILE *)fopen(fn,"w");
#endif

  do /****** this is a quasi-loop over time steps! */
  {  
    /* I am a PME node, I need the data with indices starting from                
     * pmeindex[nsb->nodeid] up to pmeindex[nsb->nodeid]+pmehomenr[nsb->nodeid]-1
     * Find out from which ppnodes I have to receive this data.
     */
    messages=0;
    for (sender=0; sender < nsb->nnodes-nsb->npmenodes; sender++) /* index of sender pp node */
    {
      /* First data index to be received from PP node i: */
      firstindex=max(nsb->index[sender], nsb->pmeindex[ME]); 
      /* Last data index to be received from PP node i: */
      lastindex =min(nsb->index[sender]+nsb->homenr[sender]-1, nsb->pmeindex[ME]+nsb->pmehomenr[ME]-1);
      if (lastindex >= firstindex)
      {
        if (MPI_Irecv(&xpme[firstindex-nsb->pmeindex[ME]][0],
	              (lastindex-firstindex+1)*DIM,GMX_MPI_REAL,sender,
                      MPI_ANY_TAG,MPI_COMM_WORLD,&req[messages]) != 0)
          gmx_comm("MPI_Irecv failed in do_pmeonly");
        messages++;
		       
	if (bFirst) /* receive charges at first time-step */
	{
          MPI_Recv(&chargeA[firstindex-nsb->pmeindex[ME]],
	           (lastindex-firstindex+1),GMX_MPI_REAL,sender,
                   MPI_ANY_TAG,MPI_COMM_WORLD,&status);
          if (bFreeEnergy)
            MPI_Recv(&chargeB[firstindex-nsb->pmeindex[ME]],
	             (lastindex-firstindex+1),GMX_MPI_REAL,sender,
                     MPI_ANY_TAG,MPI_COMM_WORLD,&status);		      
        }
/*        fprintf(stderr, "    RECEIVE: Node %d from node %d (Indices %d to %d) - %f %f %f\n",
 *                ME,sender,firstindex,lastindex,xpme[firstindex-nsb->pmeindex[ME]][0],
 *  	                                       xpme[firstindex-nsb->pmeindex[ME]][1],
 *	    		  		       xpme[firstindex-nsb->pmeindex[ME]][2]);
 */
      }
    }
    /* now we got the data! */
    bFirst=FALSE;
    MPI_Waitall(messages, req, stat);

    /* Check MPI_TAG if this was the last time step: */
    if (stat[0].MPI_TAG == 1) bDone=TRUE;
     

    for(q=0; q<(bFreeEnergy ? 2 : 1); q++) /* loop over q */
    {
      if (q == 0) {
        grid = gridA;
	chargepme = chargeA;
      } else {
        grid = gridB;
        chargepme = chargeB;
      }

      /* Unpack structure */
      unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
    
      my_homenr=nsb->natoms;

#ifdef PRT_FORCE_X
      if (count==0) {
        for(i=0; i<PMEHOMENR(nsb); i++) { 
          fprintf(fp_f4,"x %5d  %20.11e %20.11e %20.11e\n",
                 (i+PMESTART(nsb)),xpme[i][XX],xpme[i][YY],xpme[i][ZZ]);
	}       
      }
#endif
 
      if(cr->npmenodes == 1) {
        x_tmp=xpme;
        f_tmp=fpme;
        q_tmp=chargepme;
      } else {
/*        fprintf(stderr,"node %d -- PMEHOMENR(nsb)=%d -- xpme[0][0]=%f -- ir->nkx=%d\n",
 *                cr->nodeid,PMEHOMENR(nsb),xpme[0][0],ir->nkx);
 */
        GMX_MPE_LOG(ev_pmeredist_start);
        pme_calc_pidx(cr,PMEHOMENR(nsb),box,xpme,ir->nkx,nnx,pidx,gidx);
        pmeredist(cr,TRUE,PMEHOMENR(nsb), xpme, chargepme, gidx, 
                  &my_homenr, x_tmp, q_tmp);
        GMX_MPE_LOG(ev_pmeredist_finish);
      }
      /* fprintf(logfile,"Node= %6d, pme local particles=%6d\n", cr->nodeid,my_homenr); */

      /* Spread the charges on a grid */
      GMX_MPE_LOG(ev_spread_on_grid_start);
      spread_on_grid(logfile,cr,grid,my_homenr,ir->pme_order,
  		     x_tmp,q_tmp,box,bGatherOnly,
		     q==0 ? FALSE : TRUE); 
      GMX_MPE_LOG(ev_spread_on_grid_finish);

      if (!bGatherOnly) 
      {
        inc_nrnb(mynrnb,eNR_SPREADQBSP,
	         ir->pme_order*ir->pme_order*ir->pme_order*my_homenr);
      
        /* sum contributions to local grid from other nodes */
        if (cr->npmenodes > 1) 
 	  sum_qgrid_dd(cr,nsb,grid,ir->pme_order,TRUE);
        if (debug)
	  pr_fftgrid(debug,"qgrid",grid);

        /* do 3d-fft */ 
        GMX_MPE_LOG(ev_gmxfft3d_start);
        gmxfft3D(grid,GMX_FFT_REAL_TO_COMPLEX,cr);
        GMX_MPE_LOG(ev_gmxfft3d_finish);

        /* solve in k-space for our local cells */
        vol = det(box);
        GMX_MPE_LOG(ev_solve_pme_start);
/*      energy_AB[q]=solve_pme(grid,ewaldcoeff,vol,bsp_mod,recipbox,
	  		       vir_AB[q],cr); */			       
        energy_AB[q]=solve_pme(grid,ewaldcoeff,vol,ir->epsilon_r,
			       bsp_mod,recipbox,vir_AB[q],cr);
        GMX_MPE_LOG(ev_solve_pme_finish);

        inc_nrnb(mynrnb,eNR_SOLVEPME,nx*ny*nz*0.5);

        /* do 3d-invfft */
        GMX_MPE_LOG(ev_gmxfft3d_start);
        gmxfft3D(grid,GMX_FFT_COMPLEX_TO_REAL,cr);
        GMX_MPE_LOG(ev_gmxfft3d_finish);

        /* distribute local grid to all nodes */
        if (cr->npmenodes > 1) 
	  sum_qgrid_dd(cr,nsb,grid,ir->pme_order,FALSE); 
        if (debug)
	  pr_fftgrid(debug,"potential",grid);
	
        ntot  = grid->nxyz;  
        npme  = ntot*log((real)ntot)/(cr->npmenodes*log(2.0));
        inc_nrnb(mynrnb,eNR_FFT,2*npme);
      }
    
      /* interpolate forces for our local atoms */
      GMX_MPE_LOG(ev_gather_f_bsplines_start);
      gather_f_bsplines(grid,recipbox,idx,f_tmp,q_tmp,
  		        bFreeEnergy ? (q==0 ? 1.0-lambda : lambda) : 1.0,
		        theta,dtheta,my_homenr,ir->pme_order,
		        nnx,nny,nnz,
	                pmeatom);
      GMX_MPE_LOG(ev_gather_f_bsplines_finish);

      /* redistribute forces */
      if(cr->npmenodes > 1) 
      {
        GMX_MPE_LOG(ev_pmeredist_start);
        pmeredist(cr,FALSE,PMEHOMENR(nsb), fpme, chargepme, gidx, 
                  &my_homenr, f_tmp, q_tmp);
        GMX_MPE_LOG(ev_pmeredist_finish);
      }

#ifdef PRT_FORCE
      if (count==10) {
      	fprintf(stderr, "### Printing forces in do_pmeonly ###\n");

        for(i=0; i<PMEHOMENR(nsb); i++) 
	{ 
          fprintf(fp_f3,"force %5d  %20.12e %20.12e %20.12e\n",
                 (i+PMESTART(nsb)),fpme[i][XX],fpme[i][YY],fpme[i][ZZ]);
        }	
      }
#endif
  
      inc_nrnb(mynrnb,eNR_GATHERFBSP,
  	       ir->pme_order*ir->pme_order*ir->pme_order*my_homenr);
    } /* of q-loop */


    if (!bFreeEnergy) 
    {
      energy = energy_AB[0];
      copy_mat(vir_AB[0],vir);
    } else 
    {
      energy = (1.0-lambda)*energy_AB[0] + lambda*energy_AB[1];
      *dvdlambda += energy_AB[1] - energy_AB[0];
      for(i=0; i<DIM; i++)
        for(j=0; j<DIM; j++)
	  vir[i][j] = (1.0-lambda)*vir_AB[0][i][j] + lambda*vir_AB[1][i][j];
    }
    if (debug)
      fprintf(debug,"PME mesh energy: %g\n",energy);


    /* this is to help the user balance the ratio
     * of PME versus PP nodes */

    /* Determine whether or not Neighbour Searching is done on PP nodes 
     * bNS = ((ir->nstlist && (step % ir->nstlist==0))
     *        || step==0);
     */
    t3=MPI_Wtime()-t2;
    t2=MPI_Wtime();
    if TAKETIME
    {
      t0=MPI_Wtime();
      MPI_Barrier(MPI_COMM_WORLD);
      t1=MPI_Wtime();
      t_send_f += t1-t0;
      tstep_sum   += t3;
      timesteps++;
/*      fprintf(stderr," PME node %d, this time step %f, sum %f, steps %d, length t3=%f, tstep_sum=%f, ratio %f\n",
 *              cr->nodeid, t1-t0, t_send_f, timesteps, t3, tstep_sum, t_send_f/tstep_sum);
 */
      if (timesteps % 10 ==0)
        fprintf(stderr, "PME node %d waits %3.0f%% of the time step.\n", cr->nodeid, 100*t_send_f/tstep_sum);
    }

    /* Now the evaluated forces have to be transferred to the PP nodes */
    messages=0;
    for (receiver=0; receiver < nsb->nnodes-nsb->npmenodes; receiver++)
    {
      firstindex=max(nsb->index[receiver], nsb->pmeindex[ME]);
      lastindex =min(nsb->index[receiver]+nsb->homenr[receiver]-1,
                     nsb->pmeindex[ME]+nsb->pmehomenr[ME]-1);
      if (lastindex >= firstindex)
      {
/*        fprintf(stderr, "SEND BACK f: Node %d to node %d (Indices %d to %d) - %12.4e %12.4e %12.4e\n",
 *	        ME,receiver,firstindex,lastindex,
 *       	        fpme[firstindex-nsb->pmeindex[ME]][XX],
 *		fpme[firstindex-nsb->pmeindex[ME]][YY],
 *		fpme[firstindex-nsb->pmeindex[ME]][ZZ]);
 */
	if (MPI_Isend(&fpme[firstindex-nsb->pmeindex[ME]][0],
	              (lastindex-firstindex+1)*DIM,GMX_MPI_REAL,
		      receiver,0,MPI_COMM_WORLD,&req[messages]) != 0)
	  gmx_comm("MPI_Isend failed in do_pmeonly");
        messages++;
      }
    }
    MPI_Waitall(messages, req, stat);

    /* send energy, virial */
    /* task i out of the PME nodes group sends to task i out of the PP nodes group.
     * since we normally have more PP nodes than PME nodes some of the PP nodes will 
     * get no values. This is no problem since all these values are summed anyway in 
     * global_stat 
     */
    receiver=cr->nodeid-(cr->nnodes-cr->npmenodes);

    /* send energy 
     *
     * The tag of the energy message is used to carry some extra information
     * upon wether the user has send a TERM or USR1 signal to one of the PME 
     * nodes */
    if (bGotTermSignal) {
      tag = 1;
      bGotTermSignal=FALSE;
    }
    else if (bGotUsr1Signal) {
      tag = 2;
      bGotUsr1Signal=FALSE;
    }
    else
      tag = 0;
 
    MPI_Send(&energy,1,GMX_MPI_REAL,receiver,tag,MPI_COMM_WORLD); 

    /* for virial first fill send buffer with matrix values */
    elements=0;
    for (i=0; i<DIM; i++)
      for (j=0; j<DIM; j++)
        {
	  vbuf[elements]=vir[i][j];
	  elements++;
	}
	
    /* send virial buffer */
    MPI_Send(&vbuf,elements,GMX_MPI_REAL,receiver,2,MPI_COMM_WORLD);	

    count++;

    /* These variables are needed clean in the next time step: 
     *
     * zero out f_tmp vector, since gather_f_bsplines
     * only subsequently adds to f_tmp:
     */
    clear_rvecs(my_homenr,f_tmp);
    /* zero out fpme vector */
    clear_rvecs(PMEHOMENR(nsb),fpme);
    
    /* In case of dynamic pressure coupling, here the box for the
     * next time step is received: */
    if (DYNAMIC_BOX(*ir))
    {
      /* receive buffer with box data from last PP node */
      sender = nsb->nnodes-nsb->npmenodes-1;
      /* I receive in case I am the first PME node */
      receiver = sender+1;
      if (nsb->nodeid == receiver)
        MPI_Recv(&vbuf,DIM*DIM,GMX_MPI_REAL,sender,0,MPI_COMM_WORLD,&status);
      /* Broadcast box to all other nodes in the PME group */
      /* attention: receiver has rank=0 in PME group mpi_comm_mygroup! */
      MPI_Bcast(&vbuf, DIM*DIM, GMX_MPI_REAL, 0, cr->mpi_comm_mygroup);
      /* copy box buffer into box variable */
      elements=0;
      for (i=0; i<DIM; i++)
        for (j=0; j<DIM; j++)
          {
	    box[i][j] = vbuf[elements];
	    elements++;
	  }
    }
    /* Keep track of time step */
    step++;
    MPI_Barrier(MPI_COMM_WORLD); /* 100 */
  } /***** end of quasi-loop */
  while (!bDone);

//  fprintf(stderr,"=== node %d received TAG=1 message!\n", ME);
     
#undef ME
}

/******************/
#endif /* GMX_MPI */  
/******************/


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
  static FILE *fp_f, *fp_f4;
  static char fn[11];
  static int count=0;   
  
  
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

#ifdef PRT_FORCE
    sprintf(fn,"f_do_pme%2.2d",cr->nodeid);
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

#ifdef PRT_FORCE_X
    if (count==0)
      for (i=START(nsb); (i<START(nsb)+HOMENR(nsb)); i++) { 
        fprintf(fp_f4,"x %5d  %20.11e %20.11e %20.11e\n",
                i,x[i][XX],x[i][YY],x[i][ZZ]);
      }
#endif

    if (!bPar) {
      x_tmp=x;
      f_tmp=f;
      q_tmp=homecharge;
    } else {
      pme_calc_pidx(cr,HOMENR(nsb),box,x+START(nsb),ir->nkx,nnx,pidx,gidx);
      where();

      GMX_BARRIER(MPI_COMM_WORLD);
      pmeredist(cr, TRUE, HOMENR(nsb), x+START(nsb), homecharge, gidx, 
		&my_homenr, x_tmp, q_tmp);
      where();
    }
    where();
    if (debug)
      fprintf(debug,"Node= %6d, pme local particles=%6d\n",
	      cr->nodeid,my_homenr);

    /* Spread the charges on a grid */
    GMX_MPE_LOG(ev_spread_on_grid_start);

    /* Spread the charges on a grid */
    spread_on_grid(logfile,cr,grid,my_homenr,ir->pme_order,
		   x_tmp,q_tmp,box,bGatherOnly,
		   q==0 ? FALSE : TRUE);
    GMX_MPE_LOG(ev_spread_on_grid_finish);

    if (!bGatherOnly) {
      inc_nrnb(nrnb,eNR_SPREADQBSP,
	       ir->pme_order*ir->pme_order*ir->pme_order*my_homenr);
      
      /* sum contributions to local grid from other nodes */
      if (bPar) {
        GMX_BARRIER(MPI_COMM_WORLD);
	sum_qgrid_dd(cr,nsb,grid,ir->pme_order,TRUE);
	where();
      }
#ifdef DEBUG
      if (debug)
	pr_fftgrid(debug,"qgrid",grid);
#endif
      where();

      /* do 3d-fft */ 
      GMX_BARRIER(MPI_COMM_WORLD);
      GMX_MPE_LOG(ev_gmxfft3d_start);
      gmxfft3D(grid,GMX_FFT_REAL_TO_COMPLEX,cr);
      GMX_MPE_LOG(ev_gmxfft3d_finish);

      where();

      /* solve in k-space for our local cells */
      vol = det(box);
      GMX_BARRIER(MPI_COMM_WORLD);
      GMX_MPE_LOG(ev_solve_pme_start);
      energy_AB[q]=solve_pme(grid,ewaldcoeff,vol,ir->epsilon_r,
			     bsp_mod,recipbox,vir_AB[q],cr);
      where();
      GMX_MPE_LOG(ev_solve_pme_finish);
      inc_nrnb(nrnb,eNR_SOLVEPME,nx*ny*nz*0.5);

      /* do 3d-invfft */
      GMX_BARRIER(MPI_COMM_WORLD);
      GMX_MPE_LOG(ev_gmxfft3d_start);
      where();
      gmxfft3D(grid,GMX_FFT_COMPLEX_TO_REAL,cr);
      where();
      GMX_MPE_LOG(ev_gmxfft3d_finish);
      
      /* distribute local grid to all nodes */
      if (bPar) {
        GMX_BARRIER(MPI_COMM_WORLD);
	sum_qgrid_dd(cr,nsb,grid,ir->pme_order,FALSE);
      }
      where();
#ifdef DEBUG
      if (debug)
	pr_fftgrid(debug,"potential",grid);
#endif

      ntot  = grid->nxyz;  
      npme  = ntot*log((real)ntot)/log(2.0);
      if (bPar)
	npme /= (cr->nnodes - cr->npmenodes);
      inc_nrnb(nrnb,eNR_FFT,2*npme);
      where();
    }
    /* interpolate forces for our local atoms */
    GMX_BARRIER(MPI_COMM_WORLD);
    GMX_MPE_LOG(ev_gather_f_bsplines_start);
    
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
		      nnx,nny,nnz,pmeatom);
    where();
    
    GMX_MPE_LOG(ev_gather_f_bsplines_finish);
    
    if (bPar) {
      GMX_BARRIER(MPI_COMM_WORLD);
      pmeredist(cr, FALSE,HOMENR(nsb), f+START(nsb), homecharge, gidx, 
		&my_homenr, f_tmp, q_tmp);
    }
    where();
    
#ifdef PRT_FORCE
    if (count==10) 
    {
      fprintf(stderr, "*** Printing forces in do_pme ***\n");  

      for(i=START(nsb); (i<START(nsb)+HOMENR(nsb)); i++) { 
/*      fprintf(fp_f,"force %5d  %20.12e %20.12e %20.12e %20.12e\n",
                      i,f[i][XX],f[i][YY],f[i][ZZ],homecharge[i-START(nsb)]); */
        fprintf(fp_f,"force %5d  %20.12e %20.12e %20.12e\n",
                i,f[i][XX],f[i][YY],f[i][ZZ]);
      }
    }
#endif
    count++;
    inc_nrnb(nrnb,eNR_GATHERFBSP,
	     ir->pme_order*ir->pme_order*ir->pme_order*my_homenr);
  } /* of q-loop */
  
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
