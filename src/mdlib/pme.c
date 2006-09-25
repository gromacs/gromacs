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
#include "coulomb.h"
#include "fftgrid.h"
#include "gmx_fatal.h"
#include "pme.h"
#include "network.h"
#include "physics.h"
#include "nrnb.h"
#include "copyrite.h"
#include "domdec.h"
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
  int atom,index;
} t_pmeatom;

typedef struct gmx_pme {
  int  nodeid;             /* Our nodeid in mpi->mpi_comm */
  int  nnodes;             /* The number of nodes doing PME */
#ifdef GMX_MPI
  MPI_Comm mpi_comm;
#endif

  bool bFEP;               /* Compute Free energy contribution */
  bool bVerbose;           /* Should we be talkative? */
  int nleftbnd,nrightbnd;  /* The number of nodes to communicate with */
  int *leftbnd,*rightbnd;  /* Size of left and right boundary */
  int *leftid,*rightid;
  int my_homenr;
  int nkx,nky,nkz;         /* Grid dimensions */
  int pme_order;
  real epsilon_r;           
  t_fftgrid *gridA,*gridB;
  int  *nnx,*nny,*nnz;
  ivec *idx;
  rvec *fractx;            /* Fractional coordinate relative to the
			    * lower cell boundary 
			    */
  rvec      *x_home;
  real      *q_home;
  rvec      *f_home;
  matrix    recipbox;
  splinevec theta,dtheta,bsp_mod;
  t_pmeatom *pmeatom;
  int  homenr_nalloc;

  int  nnodes_comm;       /* The number of nodes to communicate with with DD */
  int  *node_dest;        /* The nodes to send x and q to with DD */
  int  *node_src;         /* The nodes to receive x and q from with DD */
  int  *count;            /* The number of atoms to send to each node */
  int  *rcount;           /* The number of atoms to receive */
  int  *buf_index;        /* Index for commnode into the buffers */
  rvec *bufv;             /* Communication buffer */
  real *bufr;             /* Communication buffer */
  int  buf_nalloc;        /* The communication buffer size */
} t_gmx_pme;
  
typedef struct {
  int    natoms;
  matrix box;
  real   lambda;
  bool   bLastStep;
} gmx_pme_comm_n_box_t;

typedef struct {
  matrix vir;
  real   energy;
  real   dvdlambda;
  int    flags;
} gmx_pme_comm_vir_ene_t;

/* The following stuff is needed for signal handling on the PME nodes. 
 * signal_handler needs to be defined in md.c, the bGot..Signal variables
 * here */ 
extern RETSIGTYPE signal_handler(int n);

volatile bool bGotTermSignal = FALSE, bGotUsr1Signal = FALSE; 

#define PME_TERM (1<<0)
#define PME_USR1 (1<<1)

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

static int comp_pmeatom(const void *a,const void *b)
{
  t_pmeatom *pa,*pb;
  pa=(t_pmeatom *)a;
  pb=(t_pmeatom *)b;
  
  return pa->index-pb->index;
}

static void calc_idx(gmx_pme_t pme,rvec x[])
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
  unpack_fftgrid(pme->gridA,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  
  rxx = pme->recipbox[XX][XX];
  ryx = pme->recipbox[YY][XX];
  ryy = pme->recipbox[YY][YY];
  rzx = pme->recipbox[ZZ][XX];
  rzy = pme->recipbox[ZZ][YY];
  rzz = pme->recipbox[ZZ][ZZ];

  for(i=0; (i<pme->my_homenr); i++) {
    xptr   = x[i];
    idxptr = pme->idx[i];
    fptr   = pme->fractx[i];
    
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

#ifdef SORTPME    
    pme->pmeatom[i].atom  = i;
    pme->pmeatom[i].index = INDEX(idxptr[XX],idxptr[YY],idxptr[ZZ]);
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
  qsort(pme->pmeatom,pme->my_homenr,sizeof(pmeatom[0]),comp_pmeatom);
  GMX_MPE_LOG(ev_sort_finish);
#endif
}

static void pme_calc_pidx(int npmenodes,
			  int natoms,matrix box, rvec x[],int nx, int nnx[],
			  int *pidx, int *gidx, int *count)
{
  int  i,j,jj;
  int nx2;
  int  tix;
  int  irange, taskid;
  real *xptr,tx;
  real rxx,ryx,rzx;
  matrix recipbox;
    
/* Calculate PME task index (pidx) for each grid index */
  irange=1+(nx-1)/npmenodes;

  taskid=0;
  for(i=0; i<nx; i+=irange) {
    jj=min((i+irange),nx);
    for(j=i; (j<jj); j++) pidx[j]=taskid;
    taskid++;
  }

  /* Reset the count */
  for(i=0; i<npmenodes; i++)
    count[i] = 0;

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
    count[gidx[i]]++;
  }
}

static void pme_realloc_homenr_things(gmx_pme_t pme)
{
  int i;

  if (pme->my_homenr > pme->homenr_nalloc) {
    pme->homenr_nalloc = over_alloc(pme->my_homenr);

    if (pme->nnodes > 1) {
      srenew(pme->x_home,pme->homenr_nalloc);
      srenew(pme->q_home,pme->homenr_nalloc);
      srenew(pme->f_home,pme->homenr_nalloc);
    }
    for(i=0;i<DIM;i++) {
      srenew(pme->theta[i],pme->pme_order*pme->homenr_nalloc); 
      srenew(pme->dtheta[i],pme->pme_order*pme->homenr_nalloc);
    }
    srenew(pme->fractx,pme->homenr_nalloc); 
    srenew(pme->pmeatom,pme->homenr_nalloc);
    srenew(pme->idx,pme->homenr_nalloc);
  }
}

static void pmeredist(gmx_pme_t pme, bool forw,
		      int n, rvec *x_f, real *charge,int *idxa)
/* Redistribute particle data for PME calculation */
/* domain decomposition by x coordinate           */
{
  static int *scounts, *rcounts,*sdispls, *rdispls, *sidx;
  static real *buf=NULL;
  static int buf_nalloc=0;
  int i, ii;
  static bool bFirst=TRUE;
  static int  npme; /* how many nodes participate in PME? */

  if(bFirst) {
    npme = pme->nnodes;

    snew(scounts,npme);
    snew(rcounts,npme);
    snew(sdispls,npme);
    snew(rdispls,npme);
    snew(sidx,npme);
    bFirst=FALSE;
  }
  if (n > buf_nalloc) {
    buf_nalloc = over_alloc(n);
    srenew(buf,buf_nalloc*DIM);
  }

#ifdef GMX_MPI
  if (forw && x_f) { /* forward, redistribution from pp to pme */ 

/* Calculate send counts and exchange them with other nodes */
    for(i=0; (i<npme); i++) scounts[i]=0;
    for(i=0; (i<n); i++) scounts[idxa[i]]++;
    MPI_Alltoall( scounts, 1, MPI_INT, rcounts, 1, MPI_INT, pme->mpi_comm);

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
    pme->my_homenr = rdispls[npme-1] + rcounts[npme-1];

    pme_realloc_homenr_things(pme);

/* Copy particle coordinates into send buffer and exchange*/
    for(i=0; (i<n); i++) {
      ii=DIM*sidx[idxa[i]];
      sidx[idxa[i]]++;
      buf[ii+XX]=x_f[i][XX];
      buf[ii+YY]=x_f[i][YY];
      buf[ii+ZZ]=x_f[i][ZZ];
    }
    MPI_Alltoallv(buf, scounts, sdispls, rvec_mpi,
                  pme->x_home, rcounts, rdispls, rvec_mpi,
                  pme->mpi_comm);
  }
  if (forw) {
/* Copy charge into send buffer and exchange*/
    for(i=0; i<npme; i++) sidx[i]=sdispls[i];
    for(i=0; (i<n); i++) {
      ii=sidx[idxa[i]];
      sidx[idxa[i]]++;
      buf[ii]=charge[i];
    }
    MPI_Alltoallv(buf, scounts, sdispls, mpi_type,
                  pme->q_home, rcounts, rdispls, mpi_type,
                  pme->mpi_comm);
  }
  else { /* backward, redistribution from pme to pp */ 
    MPI_Alltoallv(pme->f_home, rcounts, rdispls, rvec_mpi,
                  buf, scounts, sdispls, rvec_mpi,
                  pme->mpi_comm);

    /* Copy data from receive buffer */
    for(i=0; i<npme; i++)
      sidx[i] = sdispls[i];
    for(i=0; (i<n); i++) {
      ii = DIM*sidx[idxa[i]];
      x_f[i][XX] = buf[ii+XX];
      x_f[i][YY] = buf[ii+YY];
      x_f[i][ZZ] = buf[ii+ZZ];
      sidx[idxa[i]]++;
    }
  }
#endif 
}

static void pme_dd_sendrecv(gmx_pme_t pme,bool bBackward,int shift,
			    void *buf_s,int nbyte_s,
			    void *buf_r,int nbyte_r)
{
#ifdef GMX_MPI
  int dest,src;
  MPI_Status stat;
  
  if (bBackward == FALSE) {
    dest = pme->node_dest[shift];
    src  = pme->node_src[shift];
  } else {
    dest = pme->node_src[shift];
    src  = pme->node_dest[shift];
  }

  if (nbyte_s > 0 && nbyte_r > 0) {
    MPI_Sendrecv(buf_s,nbyte_s,MPI_BYTE,
		 dest,shift,
		 buf_r,nbyte_r,MPI_BYTE,
		 src,shift,
		 pme->mpi_comm,&stat);
  } else if (nbyte_s > 0) {
    MPI_Send(buf_s,nbyte_s,MPI_BYTE,
	     dest,shift,
	     pme->mpi_comm);
  } else if (nbyte_r > 0) {
    MPI_Recv(buf_r,nbyte_r,MPI_BYTE,
	     src,shift,
	     pme->mpi_comm,&stat);
  }
#endif
}
static void dd_pmeredist_x_q(gmx_pme_t pme,
			     int n, rvec *x, real *charge, int *idxa)
{
  int *commnode,*buf_index;
  int i,nsend,local_pos,buf_pos,node,scount;

  commnode  = pme->node_dest;
  buf_index = pme->buf_index;
  
  nsend = 0;
  for(i=0; i<pme->nnodes_comm; i++) {
    buf_index[commnode[i]] = nsend;
    nsend += pme->count[commnode[i]];
  }
  if (x) {
    if (pme->count[pme->nodeid] + nsend != n)
      gmx_fatal(FARGS,"%d particles communicated to PME node %d are more than a cell length out of the domain decomposition cell of their charge group",
		n - (pme->count[pme->nodeid] + nsend),pme->nodeid);
    
    if (nsend > pme->buf_nalloc) {
      pme->buf_nalloc = over_alloc(nsend);
      srenew(pme->bufv,pme->buf_nalloc);
      srenew(pme->bufr,pme->buf_nalloc);
    }
    
    pme->my_homenr = pme->count[pme->nodeid];
    for(i=0; i<pme->nnodes_comm; i++) {
      scount = pme->count[commnode[i]];
      /* Communicate the count */
      if (debug)
	fprintf(debug,"PME node %d send to node %d: %d\n",
		pme->nodeid,commnode[i],scount);
      pme_dd_sendrecv(pme,FALSE,i,
		      &scount,sizeof(int),
		      &pme->rcount[i],sizeof(int));
      pme->my_homenr += pme->rcount[i];
    }

    pme_realloc_homenr_things(pme);
  }
  
  local_pos = 0;
  for(i=0; i<n; i++) {
    node = idxa[i];
    if (node == pme->nodeid) {
      /* Copy direct to the receive buffer */
      if (x)
	copy_rvec(x[i],pme->x_home[local_pos]);
      pme->q_home[local_pos] = charge[i];
      local_pos++;
    } else {
      /* Copy to the send buffer */
      if (x)
	copy_rvec(x[i],pme->bufv[buf_index[node]]);
      pme->bufr[buf_index[node]] = charge[i];
      buf_index[node]++;
    }
  }

  buf_pos = 0;
  for(i=0; i<pme->nnodes_comm; i++) {
    scount = pme->count[commnode[i]];
    if (x) {
      /* Communicate the coordinates */
      pme_dd_sendrecv(pme,FALSE,i,
		      pme->bufv[buf_pos],scount*sizeof(rvec),
		      pme->x_home[local_pos],pme->rcount[i]*sizeof(rvec));
    }
    /* Communicate the charges */
    pme_dd_sendrecv(pme,FALSE,i,
		    pme->bufr+buf_pos,scount*sizeof(real),
		    pme->q_home+local_pos,pme->rcount[i]*sizeof(real));
    buf_pos   += scount;
    local_pos += pme->rcount[i];
  }
}

static void dd_pmeredist_f(gmx_pme_t pme, int n, rvec *f,int *idxa)
{
  int *commnode,*buf_index;
  int local_pos,buf_pos,i,rcount,node;

  commnode  = pme->node_dest;
  buf_index = pme->buf_index;

  local_pos = pme->count[pme->nodeid];
  buf_pos = 0;
  for(i=0; i<pme->nnodes_comm; i++) {
    rcount = pme->count[commnode[i]];
    /* Communicate the forces */
    pme_dd_sendrecv(pme,TRUE,i,
		    pme->f_home[local_pos],pme->rcount[i]*sizeof(rvec),
		    pme->bufv[buf_pos],rcount*sizeof(rvec));
    local_pos += pme->rcount[i];
    buf_index[commnode[i]] = buf_pos;
    buf_pos   += rcount;
  }

  local_pos = 0;
  for(i=0; i<n; i++) {
    node = idxa[i];
    if (node == pme->nodeid) {
      /* Copy from the local force array */
      copy_rvec(pme->f_home[local_pos],f[i]);
      local_pos++;
    } else {
      /* Copy from the receive buffer */
      copy_rvec(pme->bufv[buf_index[node]],f[i]);
      buf_index[node]++;
    }
  }
}

static void gmx_sum_qgrid_dd(gmx_pme_t pme,t_fftgrid *grid,
			     int direction)
{
  static bool bFirst=TRUE;
  static real *tmp;
  int b,i;
  static int localsize;
  static int maxproc;
  int ny, la2r;
  int bndsize;
  real *from, *to;
#ifdef GMX_MPI
  MPI_Status stat;
#endif

  GMX_MPE_LOG(ev_sum_qgrid_start);
  
#ifdef GMX_MPI
  if(bFirst) {
    bFirst=FALSE;
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
  if (grid->workspace)
    tmp=grid->workspace;
  if (direction  == GMX_SUM_QGRID_FORWARD) { 
    /* sum contributions to local grid */
    ny=grid->ny;
    la2r=grid->la2r;
    /* Send left boundaries */
    for(b=0; b<pme->nleftbnd; b++) {
      bndsize = (pme->leftbnd[b+1] - pme->leftbnd[b])*ny*la2r;
      from = grid->ptr + (pme->leftid[b] + 1)*localsize - bndsize;
      to   = grid->ptr + (pme->nodeid    + 1)*localsize - bndsize;
      MPI_Sendrecv(from,bndsize,mpi_type,
		   pme->leftid[b],pme->nodeid,
		   tmp,bndsize,mpi_type,pme->rightid[b],pme->rightid[b],
		   pme->mpi_comm,&stat);
      GMX_MPE_LOG(ev_test_start); 
      for(i=0; (i<bndsize); i++) {
	to[i] += tmp[i];
      }
    }
    GMX_MPE_LOG(ev_test_finish);
    /* Send right boundaries */
    for(b=0; b<pme->nrightbnd; b++) {
      bndsize = (pme->rightbnd[b+1] - pme->rightbnd[b])*ny*la2r;
      from = grid->ptr + (pme->rightid[b])*localsize;
      to   = grid->ptr + (pme->nodeid    )*localsize;
      MPI_Sendrecv(from,bndsize,mpi_type,
		   pme->rightid[b],pme->nodeid,
		   tmp,bndsize,mpi_type,pme->leftid[b],pme->leftid[b],
		   pme->mpi_comm,&stat); 
      GMX_MPE_LOG(ev_test_start);
      for(i=0; (i<bndsize); i++) {
	to[i] += tmp[i];
      }
    }
    GMX_MPE_LOG(ev_test_finish);
  }
  else if (direction  == GMX_SUM_QGRID_BACKWARD) { 
    /* distribute local grid to all processors */
    ny=grid->ny;
    la2r=grid->la2r;
    /* Send right boundaries */
    for(b=0; b<pme->nrightbnd; b++) {
      bndsize = (pme->rightbnd[b+1] - pme->rightbnd[b])*ny*la2r;
      from = grid->ptr + (pme->nodeid    )*localsize;
      to   = grid->ptr + (pme->rightid[b])*localsize;
      MPI_Sendrecv(from,bndsize,mpi_type,
		   pme->leftid[b],pme->nodeid,
		   to,bndsize,mpi_type,pme->rightid[b],pme->rightid[b],
		   pme->mpi_comm,&stat);
    }
    /* Send left boundaries */
    for(b=0; b<pme->nleftbnd; b++) {
      bndsize = (pme->leftbnd[b+1] - pme->leftbnd[b])*ny*la2r;
      from = grid->ptr + (pme->nodeid    + 1)*localsize - bndsize;
      to   = grid->ptr + (pme->leftid[b] + 1)*localsize - bndsize;
      MPI_Sendrecv(from,bndsize,mpi_type,
		   pme->rightid[b],pme->nodeid,
		   to,bndsize,mpi_type,pme->leftid[b],pme->leftid[b],
		   pme->mpi_comm,&stat);
    }
  }
  else
    gmx_fatal(FARGS,"Invalid direction %d for summing qgrid",direction);
    
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
  if (grid->workspace)
    tmp=grid->workspace;
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

    if(cr->nodeid<maxproc)
      memcpy(grid->ptr+cr->nodeid*localsize,tmp,localsize*sizeof(real));
  }
  else if (direction == GMX_SUM_QGRID_BACKWARD) { 
    /* distribute local grid to all processors */
    for(i=0;i<maxproc;i++)
      MPI_Bcast(grid->ptr+i*localsize, /* ptr arithm     */
		localsize,       
		GMX_MPI_REAL,i,cr->mpi_comm_mygroup);
  }
  else
    gmx_fatal(FARGS,"Invalid direction %d for summing qgrid",direction);
    
#else
  gmx_fatal(FARGS,"Parallel grid summation requires MPI and FFTW.\n");    
#endif
}

static void spread_q_bsplines(gmx_pme_t pme,t_fftgrid *grid,real charge[],
			      int nr,int order)
{
  /* spread charges from home atoms to local grid */
  real     *ptr;
  int      b,i,nn,n,*i0,*j0,*k0,*ii0,*jj0,*kk0,ithx,ithy,ithz;
  int      nx,ny,nz,nx2,ny2,nz2,la2,la12;
  int      norder,*idxptr,index_x,index_xy,index_xyz;
  real     valx,valxy,qn;
  real     *thx,*thy,*thz;
  int localsize, bndsize;
  
  if (pme->nnodes <= 1) {
    clear_fftgrid(grid); 
#ifdef GMX_MPI
  } else {
    localsize = grid->la12r*grid->pfft.local_nx;
    ptr = grid->localptr;
    for (i=0; (i<localsize); i++)
      ptr[i] = 0;
    /* clear left boundary area */
    for(b=0; b<pme->nleftbnd; b++) {
      bndsize = (pme->leftbnd[b+1] - pme->leftbnd[b])*grid->la12r;
      ptr = grid->ptr + (pme->leftid[b] + 1)*localsize - bndsize;
      for (i=0; (i<bndsize); i++)
	ptr[i] = 0;
    }
    /* clear right boundary area */
    for(b=0; b<pme->nrightbnd; b++) {
      bndsize = (pme->rightbnd[b+1] - pme->rightbnd[b])*grid->la12r;
      ptr = grid->ptr + (pme->rightid[b])*localsize;
      for (i=0; (i<bndsize); i++)
	ptr[i] = 0;
    }
#endif
  }
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  ii0   = pme->nnx+nx2+1-order/2;
  jj0   = pme->nny+ny2+1-order/2;
  kk0   = pme->nnz+nz2+1-order/2;
  
  for(nn=0; (nn<nr); nn++) {
#ifdef SORTPME
    n = pme->pmeatom[nn].atom;
#else
    n = nn;
#endif
    qn     = charge[n];
    idxptr = pme->idx[n];
    
    if (qn != 0) {
      norder  = n*order;
      
      /* Pointer arithmetic alert, next six statements */
      i0  = ii0 + idxptr[XX]; 
      j0  = jj0 + idxptr[YY];
      k0  = kk0 + idxptr[ZZ];
      thx = pme->theta[XX] + norder;
      thy = pme->theta[YY] + norder;
      thz = pme->theta[ZZ] + norder;
      
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
    
  if (pme->nnodes > 1) { 
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
    
    if(ky<maxky)
      my = ky;
    else
      my = (ky-ny);
    by = M_PI*vol*pme->bsp_mod[YY][ky];
    
    for(kx=0; (kx<nx); kx++) {    
      if(kx < maxkx)
	mx = kx;
      else
	mx = (kx-nx);

      mhx = mx * rxx;
      mhy = mx * ryx + my * ryy;

      bx = pme->bsp_mod[XX][kx];
      
      if (pme->nnodes > 1)
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
	denom   = m2*bx*by*pme->bsp_mod[ZZ][kz];
	eterm   = ONE_4PI_EPS0*exp(-factor*m2)/(pme->epsilon_r*denom);
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

void gather_f_bsplines(gmx_pme_t pme,t_fftgrid *grid,
		       bool bClear,rvec f[],real *charge,real scale)
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
  thx   = pme->theta[XX];
  thy   = pme->theta[YY];
  thz   = pme->theta[ZZ];
  dthx  = pme->dtheta[XX];
  dthy  = pme->dtheta[YY];
  dthz  = pme->dtheta[ZZ];
  ii0   = pme->nnx+nx2+1-order/2;
  jj0   = pme->nny+ny2+1-order/2;
  kk0   = pme->nnz+nz2+1-order/2;
  
  rxx   = pme->recipbox[XX][XX];
  ryx   = pme->recipbox[YY][XX];
  ryy   = pme->recipbox[YY][YY];
  rzx   = pme->recipbox[ZZ][XX];
  rzy   = pme->recipbox[ZZ][YY];
  rzz   = pme->recipbox[ZZ][ZZ];

  for(nn=0; (nn<pme->my_homenr); nn++) {
#ifdef SORTPME
    n = pme->pmeatom[nn].atom;
#else
    n = nn;
#endif
    qn      = scale*charge[n];

    if (bClear) {
      f[n][XX] = 0;
      f[n][YY] = 0;
      f[n][ZZ] = 0;
    }
    if (qn != 0) {
      fx     = 0;
      fy     = 0;
      fz     = 0;
      idxptr = pme->idx[n];
      norder = n*order;
      
      /* Pointer arithmetic alert, next nine statements */
      i0   = ii0 + idxptr[XX]; 
      j0   = jj0 + idxptr[YY];
      k0   = kk0 + idxptr[ZZ];
      thx  = pme->theta[XX] + norder;
      thy  = pme->theta[YY] + norder;
      thz  = pme->theta[ZZ] + norder;
      dthx = pme->dtheta[XX] + norder;
      dthy = pme->dtheta[YY] + norder;
      dthz = pme->dtheta[ZZ] + norder;
      
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
      f[n][XX] += -qn*( fx*nx*rxx );
      f[n][YY] += -qn*( fx*nx*ryx + fy*ny*ryy );
      f[n][ZZ] += -qn*( fx*nx*rzx + fy*ny*rzy + fz*nz*rzz );
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
  real dr,div;
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

static void get_my_ddnodes(FILE *logfile,t_commrec *cr,int pmenodeid,
			   int *nmy_ddnodes,int **my_ddnodes)
{
  int i;

  snew(*my_ddnodes,(cr->dd->nnodes+cr->npmenodes-1)/cr->npmenodes);
  
  *nmy_ddnodes = 0;
  for(i=0; i<cr->dd->nnodes; i++) {
    if (gmx_ddindex2pmeslab(cr,i) == pmenodeid)
      (*my_ddnodes)[(*nmy_ddnodes)++] = gmx_ddindex2nodeid(cr,i);
  }

  fprintf(logfile,"PME node %d, receive coordinates from %d PP nodes\n",
	  cr->nodeid,*nmy_ddnodes);
}

static void setup_coordinate_communication(FILE *log, t_commrec *cr,
					   gmx_pme_t pme)
{
  static bool bFirst=TRUE;
  int  npme,shmax,i,slab,dd_cx0,dd_cx1,fw,bw;
  ivec coords;
  bool bPPnode;

  npme = pme->nnodes;

  /* This code assumes a uniform DD cell size in x direction
   * and that particles are no further than one DD cell size
   * out of their charge group home cell due to charge group size
   * and diffusion between neighbor list updates.
   */
  shmax = 1;
  for(i=0; i<cr->dd->nnodes; i++) {
    gmx_ddindex2xyz(cr->dd->nc,i,coords);
    slab = gmx_ddindex2pmeslab(cr,i);
    /* Initial (ns step) charge group center x in range 0 - dd_nx */
    dd_cx0 = coords[XX];
    dd_cx1 = coords[XX] + 1;
    /* Add one DD cell size */
    dd_cx0--;
    dd_cx1++;
    /* Now we multiply with npme, so the x range is 0 - dd_nx*npme */
    dd_cx0 *= npme;
    dd_cx1 *= npme;
    /* Check if we need to increase the maximum shift */
    while (dd_cx1 > (slab + 1 + shmax)*cr->dd->nc[XX])
      shmax++;
    while (dd_cx0 < (slab - shmax)*cr->dd->nc[XX])
      shmax++;
}

  if (bFirst) {
    fprintf(log,"PME maximum node shift for coordinate communication: %d\n",
	    shmax);
    bFirst = FALSE;
  }

  /* Set the communication ids.
   * We interleave forward and backward directions,
   * as for large shifts often only one of the two will
   * occur on each node and we can therefore avoid nodes waiting
   * while other nodes are communicating.
   */
  pme->nnodes_comm = 0;
  for(i=1; i<=shmax; i++) {
    fw = (pme->nodeid + i) % npme;
    bw = (pme->nodeid - i + npme) % npme;
    if (pme->nnodes_comm < npme - 1) {
      pme->node_dest[pme->nnodes_comm] = fw;
      pme->node_src[pme->nnodes_comm] = bw;
      pme->nnodes_comm++;
    } 
    if (pme->nnodes_comm < npme - 1) {
      pme->node_dest[pme->nnodes_comm] = bw;
      pme->node_src[pme->nnodes_comm] = fw;
      pme->nnodes_comm++;
    }
  }
}

int gmx_pme_destroy(FILE *log,gmx_pme_t *pmedata)
{
  fprintf(log,"Destroying PME data structures.\n");
  sfree((*pmedata)->nnx);
  sfree((*pmedata)->nny);
  sfree((*pmedata)->nnz);
  sfree((*pmedata)->idx);
  sfree((*pmedata)->fractx);
  sfree((*pmedata)->pmeatom);
  done_fftgrid((*pmedata)->gridA);
  if((*pmedata)->gridB)
    done_fftgrid((*pmedata)->gridB);
    
  sfree(*pmedata);
  *pmedata = NULL;
  
  return 0;
}

int gmx_pme_init(FILE *log,gmx_pme_t *pmedata,t_commrec *cr,
		 t_inputrec *ir,int homenr,
		 bool bFreeEnergy,bool bVerbose)
{
  gmx_pme_t pme=NULL;
  
  int b,d,i,totbnd;
  
  fprintf(log,"Creating PME data structures.\n");
  snew(pme,1);

  if (PAR(cr)) {
#ifdef GMX_MPI
      pme->mpi_comm = cr->mpi_comm_mygroup;
      MPI_Comm_rank(pme->mpi_comm,&pme->nodeid);
      MPI_Comm_size(pme->mpi_comm,&pme->nnodes);
#endif
  } else {
    pme->nnodes = 1;
  }

  fprintf(log,"Will do PME sum in reciprocal space.\n");
  please_cite(log,"Essman95a");

  if (ir->ewald_geometry == eewg3DC) {
    fprintf(log,"Using the Ewald3DC correction for systems with a slab geometry.\n");
    please_cite(log,"In-Chul99a");
  }

  pme->bFEP = ((ir->efep != efepNO) && bFreeEnergy);
  pme->nkx  = ir->nkx;
  pme->nky  = ir->nky;
  pme->nkz  = ir->nkz;
  pme->pme_order   = ir->pme_order;
  pme->epsilon_r   = ir->epsilon_r;
  
  if (pme->nkx <= pme->pme_order ||
      pme->nky <= pme->pme_order ||
      pme->nkz <= pme->pme_order)
    gmx_fatal(FARGS,
	      "The pme grid dimensions need to be larger than pme_order (%d)",
	      pme->pme_order);

  if (pme->nnodes > 1) {
#ifdef GMX_MPI
    MPI_Type_contiguous(DIM, mpi_type, &rvec_mpi);
    MPI_Type_commit(&rvec_mpi);
#endif

    d = pme->nkx/pme->nnodes;
    /* The left boundary */
    totbnd = ir->pme_order/2 - 1;
    pme->nleftbnd = (totbnd + (d - 1))/d;
    snew(pme->leftbnd,pme->nleftbnd+1);
    pme->leftbnd[0] = 0;
    for(b=0; b<pme->nleftbnd; b++) {
      pme->leftbnd[b+1] = min(totbnd,pme->leftbnd[b] + d);
    }
    /* The right boundary */
    totbnd = ir->pme_order - (ir->pme_order/2 - 1) - 1;
    pme->nrightbnd = (totbnd + (d - 1))/d;
    snew(pme->rightbnd,pme->nrightbnd+1);
    pme->rightbnd[0] = 0;
    for(b=0; b<pme->nrightbnd; b++) {
      pme->rightbnd[b+1] = min(totbnd,pme->rightbnd[b] + d);
    }
    totbnd = max(pme->nleftbnd,pme->nrightbnd);
    snew(pme->leftid,totbnd);
    snew(pme->rightid,totbnd);
    for(b=0; b<totbnd; b++) {
      pme->leftid[b]  = (pme->nodeid - (b + 1) + pme->nnodes) % pme->nnodes;
      pme->rightid[b] = (pme->nodeid + (b + 1)) % pme->nnodes;
    }

    fprintf(log,"Parallelized PME sum used. nkx=%d, npme=%d\n",
	    ir->nkx,pme->nnodes);
    if ((ir->nkx % pme->nnodes) != 0)
      fprintf(log,"Warning: For load balance, "
	      "fourier_nx should be divisible by the number of PME nodes\n");

    if (DOMAINDECOMP(cr)) {
      snew(pme->node_dest,pme->nnodes);
      snew(pme->node_src,pme->nnodes);
      setup_coordinate_communication(log,cr,pme);
    }
    snew(pme->count,pme->nnodes);
    snew(pme->rcount,pme->nnodes);
    snew(pme->buf_index,pme->nnodes);
  }

  /* With domain decomposition we need nnx on the PP only nodes */
  snew(pme->nnx,5*pme->nkx);
  snew(pme->nny,5*pme->nky);
  snew(pme->nnz,5*pme->nkz);
  for(i=0; (i<5*pme->nkx); i++)
    pme->nnx[i] = i % pme->nkx;
  for(i=0; (i<5*pme->nky); i++)
    pme->nny[i] = i % pme->nky;
  for(i=0; (i<5*pme->nkz); i++)
    pme->nnz[i] = i % pme->nkz;

  snew(pme->bsp_mod[XX],pme->nkx);
  snew(pme->bsp_mod[YY],pme->nky);
  snew(pme->bsp_mod[ZZ],pme->nkz);

  pme->gridA = mk_fftgrid(log,pme->nkx,pme->nky,pme->nkz,NULL,cr);
  if (bFreeEnergy)
    pme->gridB = mk_fftgrid(log,pme->nkx,pme->nky,pme->nkz,NULL,cr);
  else
    pme->gridB = NULL;
  
  make_bspline_moduli(pme->bsp_mod,pme->nkx,pme->nky,pme->nkz,pme->pme_order);
  
  if (pme->nnodes == 1) {
    pme->my_homenr = homenr;
    pme_realloc_homenr_things(pme);
  }

  *pmedata = pme;
  
  return 0;
}

static void spread_on_grid(FILE *logfile,    gmx_pme_t pme,
			   t_fftgrid *grid,  t_commrec *cr,    
			   rvec x[],
			   real charge[],    matrix box,
			   bool bGatherOnly, bool bHaveSplines)
{ 
  int nx,ny,nz,nx2,ny2,nz2,la2,la12;
  real *ptr;
  
  /* Unpack structure */
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  
  if (!bHaveSplines)
    /* Inverse box */
    calc_recipbox(box,pme->recipbox); 
  
  if (!bGatherOnly) {
    if (!bHaveSplines) {
      /* Compute fftgrid index for all atoms, with help of some extra variables */
      calc_idx(pme,x);
      
      /* make local bsplines  */
      make_bsplines(pme->theta,pme->dtheta,pme->pme_order,nx,ny,nz,
		    pme->fractx,pme->idx,charge,pme->my_homenr);
    }    

    /* put local atoms on grid. */
    spread_q_bsplines(pme,grid,charge,pme->my_homenr,pme->pme_order);
/*    pr_grid_dist(logfile,"spread",grid); */
  }
}

void gmx_pme_send_x_q(t_commrec *cr, matrix box,
		      rvec x[], real chargeA[], real chargeB[],
		      bool bFreeEnergy, real lambda, bool bLastStep)
{
  int  n,nreq;
  gmx_pme_comm_n_box_t cnb;
#ifdef GMX_MPI
  MPI_Request req[4];
#endif

  n = cr->dd->nat_home;
  cnb.natoms = n;
  copy_mat(box,cnb.box);
  cnb.lambda = lambda;
  cnb.bLastStep = bLastStep;

#ifdef GMX_MPI
  nreq = 0;
  /* Communicate bLastTime (0/1) via the MPI tag */
  MPI_Isend(&cnb,sizeof(cnb),MPI_BYTE,
	    cr->dd->pme_nodeid,0,cr->mpi_comm_mysim,&req[nreq++]);

  MPI_Isend(x[0],n*sizeof(rvec),MPI_BYTE,
	    cr->dd->pme_nodeid,1,cr->mpi_comm_mysim,&req[nreq++]);
  MPI_Isend(chargeA,n*sizeof(real),MPI_BYTE,
	    cr->dd->pme_nodeid,2,cr->mpi_comm_mysim,&req[nreq++]);
  if (bFreeEnergy)
    MPI_Isend(chargeB,n*sizeof(real),MPI_BYTE,
	      cr->dd->pme_nodeid,3,cr->mpi_comm_mysim,&req[nreq++]);

  /* We would like to postpone this wait until dd_pme_receive_f,
   * but currently the coordinates can be (un)shifted.
   * We will have to wait with this until we have separate
   * arrays for the shifted and unshifted coordinates.
   */
  MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);
#endif
}

static void receive_virial_energy(t_commrec *cr,
				  matrix vir,real *energy,real *dvdlambda) 
{
  gmx_pme_comm_vir_ene_t cve;

  if (cr->dd->pme_receive_vir_ener) {
    /* receive virial and energy */
#ifdef GMX_MPI
    MPI_Recv(&cve,sizeof(cve),MPI_BYTE,cr->dd->pme_nodeid,1,cr->mpi_comm_mysim,
	     MPI_STATUS_IGNORE);
#else
    memset(&cve,0,sizeof(cve));
#endif
	
    m_add(vir,cve.vir,vir);
    *energy = cve.energy;
    *dvdlambda += cve.dvdlambda;
 
    bGotTermSignal = (cve.flags & PME_TERM);
    bGotUsr1Signal = (cve.flags & PME_USR1);
  } else {
    *energy = 0;
  }
}

void gmx_pme_receive_f(t_commrec *cr,
		       rvec f[], matrix vir, 
		       real *energy, real *dvdlambda)
{
  static rvec *f_rec=NULL;
  static int  nalloc=0;
  int natoms,i;

  natoms = cr->dd->nat_home;

  if (natoms > nalloc) {
    nalloc = over_alloc(natoms);
    srenew(f_rec,nalloc);
  }

#ifdef GMX_MPI  
  MPI_Recv(f_rec[0],natoms*sizeof(rvec),MPI_BYTE,
	   cr->dd->pme_nodeid,0,cr->mpi_comm_mysim,
	   MPI_STATUS_IGNORE);
#endif

  for(i=0; i<natoms; i++)
    rvec_inc(f[i],f_rec[i]);
  
  receive_virial_energy(cr,vir,energy,dvdlambda);
}

int gmx_pmeonly(FILE *logfile,    gmx_pme_t pme,
		t_commrec *cr,    t_nrnb *nrnb,
		real ewaldcoeff,  bool bGatherOnly)
     
{
#ifdef GMX_MPI
  rvec  *x_pp=NULL,*f_pp=NULL;
  real  *chargeA=NULL,*chargeB=NULL,*charge_pp;
  int  nppnodes,sender,receiver;
  int  nmy_ddnodes,*my_ddnodes;
  int  ind_start, ind_end;
  gmx_pme_comm_n_box_t *cnb;
  int  nalloc_natpp=0;
  int  n,i,q;
  gmx_pme_comm_vir_ene_t cve;
  int messages; /* count the Isends or Ireceives */    
  MPI_Request *req;
  MPI_Status  *stat;
  bool bDone=FALSE;
  int count=0;
  double t0, t1, t2=0, t3, t_send_f=0.0, tstep_sum=0.0;
  unsigned long int step=0, timesteps=0;

  nppnodes = cr->nnodes - cr->npmenodes;

  /* Allocate space for the request handles, a lonely PME node may
     receive (and later send back) as many messages as there are PP nodes */
  snew(req ,3*nppnodes);
  snew(stat,3*nppnodes);

  /* Determine the DD nodes we need to talk with */
  get_my_ddnodes(logfile,cr,pme->nodeid,&nmy_ddnodes,&my_ddnodes);
  snew(cnb,nmy_ddnodes);

  init_nrnb(nrnb);

  do /****** this is a quasi-loop over time steps! */
  {  
    /* Domain decomposition */
    /* Receive the send counts and box */
    messages=0;
    for(sender=0; sender<nmy_ddnodes; sender++) {
      MPI_Irecv(cnb+sender,sizeof(cnb[0]),MPI_BYTE,
		my_ddnodes[sender],0,
		cr->mpi_comm_mysim,&req[messages++]);
    }
    MPI_Waitall(messages, req, stat);
    messages = 0;
    bDone = cnb[0].bLastStep;
    
    n = 0;
    for(sender=0; sender<nmy_ddnodes; sender++) {
      n += cnb[sender].natoms;
    }
    if (n > nalloc_natpp) {
      nalloc_natpp = over_alloc(n);
      
      srenew(x_pp,nalloc_natpp);
      srenew(chargeA,nalloc_natpp);
      if (pme->bFEP)
	srenew(chargeB,nalloc_natpp);
      srenew(f_pp,nalloc_natpp);
    }
    
    /* Receive the coordinates in place */
    i = 0;
    for(sender=0; sender<nmy_ddnodes; sender++) {
      if (cnb[sender].natoms > 0) {
	MPI_Irecv(x_pp[i],cnb[sender].natoms*sizeof(rvec),MPI_BYTE,
		  my_ddnodes[sender],1,
		  cr->mpi_comm_mysim,&req[messages++]);
	i += cnb[sender].natoms;
      }
    }
    
    /* Receive the charges in place */
    for(q=0; q<(pme->bFEP ? 2 : 1); q++) {
      if (q == 0)
	charge_pp = chargeA;
      else
	charge_pp = chargeB;
      i = 0;
      for(sender=0; sender<nmy_ddnodes; sender++) {
	if (cnb[sender].natoms > 0) {
	  MPI_Irecv(&charge_pp[i],cnb[sender].natoms*sizeof(real),MPI_BYTE,
		    my_ddnodes[sender],2+q,
		    cr->mpi_comm_mysim,&req[messages++]);
	  i += cnb[sender].natoms;
	}
      }
    }

    /* Wait for the coordinates and charges to arrive */
    MPI_Waitall(messages, req, stat);
    messages = 0;

    cve.dvdlambda = 0;
    gmx_pme_do(logfile,pme,0,n,x_pp,f_pp,chargeA,chargeB,cnb[0].box,
	       cr,nrnb,cve.vir,ewaldcoeff,
	       &cve.energy,cnb[0].lambda,&cve.dvdlambda,
	       bGatherOnly);

    t3=MPI_Wtime()-t2;
    t2=MPI_Wtime();
    if (TAKETIME)
    {
      t0=MPI_Wtime();
      MPI_Barrier(cr->mpi_comm_mysim);
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
    messages = 0;
    ind_end = 0;
    for (receiver=0; receiver<nmy_ddnodes; receiver++) {
      ind_start = ind_end;
      ind_end   = ind_start + cnb[receiver].natoms;
      if (MPI_Isend(f_pp[ind_start],(ind_end-ind_start)*sizeof(rvec),MPI_BYTE,
		    my_ddnodes[receiver],0,
		    cr->mpi_comm_mysim,&req[messages++]) != 0)
	gmx_comm("MPI_Isend failed in do_pmeonly");
    }

    /* send virial and energy to our last PP node */
    cve.flags    = 0;
    if (bGotTermSignal)
      cve.flags |= PME_TERM;
    if (bGotUsr1Signal)
      cve.flags |= PME_USR1;

    MPI_Isend(&cve,sizeof(cve),MPI_BYTE,
	      my_ddnodes[nmy_ddnodes-1],1,cr->mpi_comm_mysim,&req[messages++]);

    /* Wait for the forces to arrive */
    MPI_Waitall(messages, req, stat);

    count++;

    /* Keep track of time step */
    step++;
    /* MPI_Barrier(cr->mpi_comm_mysim); */ /* 100 */
  } /***** end of quasi-loop */
  while (!bDone);
     
#endif
  return 0;
}

int gmx_pme_do(FILE *logfile,   gmx_pme_t pme,
	       int start,       int homenr,
	       rvec x[],        rvec f[],
	       real *chargeA,   real *chargeB,
	       matrix box,	t_commrec *cr,
	       t_nrnb *nrnb,    
	       matrix vir,      real ewaldcoeff,
	       real *energy,    real lambda, 
	       real *dvdlambda, bool bGatherOnly)
{
  int     q,i,j,ntot,npme;
  int     nx,ny,nz,nx2,ny2,nz2,la12,la2;
  t_fftgrid *grid=NULL;
  real    *ptr;
  real    *charge=NULL,vol;
  real    energy_AB[2];
  matrix  vir_AB[2];
  static int  *pidx=NULL;
  static int  *gidx=NULL;  /* Grid index of particle, used for PME redist */
  static int  gidx_nalloc=0;

  if (pme->nnodes > 1) {
    /* Allocate indices for sorting the atoms */
    if (pidx == NULL)
      snew(pidx,pme->nkx);
    
    if (homenr > gidx_nalloc) {
      gidx_nalloc = over_alloc(homenr);
      srenew(gidx,gidx_nalloc);
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
      fprintf(debug,"PME: nnodes = %d, nodeid = %d\n",cr->nnodes,cr->nodeid);
      fprintf(debug,"Grid = %p\n",grid);
      if (grid == NULL)
	gmx_fatal(FARGS,"No grid!");
    }
    where();
    unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
    where();

    if (pme->nnodes == 1) {
      pme->x_home = x;
      pme->q_home = charge;
      pme->f_home = f;
    } else {
      pme_calc_pidx(pme->nnodes,homenr,box,x+start,
		    pme->nkx,pme->nnx,pidx,gidx,pme->count);
      where();

      GMX_BARRIER(cr->mpi_comm_mysim);
      /* Redistribute x (only once) and qA or qB */
      if (DOMAINDECOMP(cr)) {
	dd_pmeredist_x_q(pme, homenr, q==0 ? x : NULL, charge, gidx);
      } else {
	pmeredist(pme, TRUE, homenr, q==0 ? x+start : NULL, charge, gidx);
      }
      where();
    }

    if (debug)
      fprintf(debug,"Node= %6d, pme local particles=%6d\n",
	      cr->nodeid,pme->my_homenr);

    /* Spread the charges on a grid */
    GMX_MPE_LOG(ev_spread_on_grid_start);

    /* Spread the charges on a grid */
    spread_on_grid(logfile,pme,grid,cr,
		   pme->x_home,pme->q_home,box,bGatherOnly,
		   q==0 ? FALSE : TRUE);
    GMX_MPE_LOG(ev_spread_on_grid_finish);

    if (!bGatherOnly) {
      inc_nrnb(nrnb,eNR_SPREADQBSP,
	       pme->pme_order*pme->pme_order*pme->pme_order*pme->my_homenr);
      
      /* sum contributions to local grid from other nodes */
      if (pme->nnodes > 1) {
        GMX_BARRIER(cr->mpi_comm_mysim);
	gmx_sum_qgrid_dd(pme,grid,GMX_SUM_QGRID_FORWARD);
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
      inc_nrnb(nrnb,eNR_SOLVEPME,nx*ny*nz*0.5);

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
	gmx_sum_qgrid_dd(pme,grid,GMX_SUM_QGRID_BACKWARD);
      }
      where();
#ifdef DEBUG
      if (debug)
	pr_fftgrid(debug,"potential",grid);
#endif

      ntot  = grid->nxyz;  
      npme  = ntot*log((real)ntot)/log(2.0);
      if (pme->nnodes > 1)
	npme /= (cr->nnodes - cr->npmenodes);
      inc_nrnb(nrnb,eNR_FFT,2*npme);
      where();
    }
    /* interpolate forces for our local atoms */
    GMX_BARRIER(cr->mpi_comm_mysim);
    GMX_MPE_LOG(ev_gather_f_bsplines_start);
    
    where();
    gather_f_bsplines(pme,grid,q==0,pme->f_home,pme->q_home,
		      pme->bFEP ? (q==0 ? 1.0-lambda : lambda) : 1.0);
    where();
    
    GMX_MPE_LOG(ev_gather_f_bsplines_finish);
    
    inc_nrnb(nrnb,eNR_GATHERFBSP,
	     pme->pme_order*pme->pme_order*pme->pme_order*pme->my_homenr);
  } /* of q-loop */

  if (pme->nnodes > 1) {
    GMX_BARRIER(cr->mpi_comm_mysim);
    if (DOMAINDECOMP(cr)) {
      dd_pmeredist_f(pme, homenr, f, gidx);
    } else {
      pmeredist(pme, FALSE, homenr, f+start, NULL, gidx);
    }
  }
  where();

  if (!pme->bFEP) {
    *energy = energy_AB[0];
    copy_mat(vir_AB[0],vir);
  } else {
    *energy = (1.0-lambda)*energy_AB[0] + lambda*energy_AB[1];
    *dvdlambda += energy_AB[1] - energy_AB[0];
    for(i=0; i<DIM; i++)
      for(j=0; j<DIM; j++)
	vir[i][j] = (1.0-lambda)*vir_AB[0][i][j] + lambda*vir_AB[1][i][j];
  }
  if (debug)
    fprintf(debug,"PME mesh energy: %g\n",*energy);

  return 0;
}
