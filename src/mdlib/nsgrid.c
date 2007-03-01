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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "nsgrid.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "network.h"
#include "domdec.h"
#include "partdec.h"

/***********************************
 *         Grid Routines
 ***********************************/

static void init_range_check()
{
  sprintf(warn_buf,"Explanation: During neighborsearching, we assign each particle to a grid\n"
	  "based on its coordinates. If your system contains collisions or parameter\n"
	  "errors that give particles very high velocities you might end up with some\n"
	  "coordinates being +-Infinity or NaN (not-a-number). Obviously, we cannot\n"
	  "put these on a grid, so this is usually where we detect those errors.\n"
	  "Make sure your system is properly energy-minimized and that the potential\n"
	  "energy seems reasonable before trying again.\n");
}

static void set_grid_sizes(FILE *fplog,
			   int ePBC,matrix box,real rlist,int delta,
			   gmx_domdec_t *dd,
			   t_grid *grid)
{
  int  i,j;
  bool bDD,bDDRect;
  rvec dd_cell_size;
  real size,radd,add_tric;

  if (dd) {
    for(i=0; (i<DIM); i++)
      dd_cell_size[i] = dd->cell_ns_x1[i] - dd->cell_ns_x0[i];
  }

  clear_rvec(grid->cell_offset);
  for(i=0; (i<DIM); i++) {
    bDD = dd && (dd->nc[i] > 1);
    if (!bDD) {
      bDDRect = FALSE;
      size = box[i][i];
    } else {
      /* With DD grid cell jumps only the first decomposition
       * direction has uniform DD cell boundaries.
       */
      bDDRect = !(dd->tric_dir[i] || (dd->bGridJump && i != dd->dim[0]));

      if (bDDRect) {
	radd = rlist;
      } else {
	radd = rlist/dd->skew_fac[i];
      }
      /* With DD we only need a grid of one DD cell size + rlist */
      grid->cell_offset[i] = dd->cell_ns_x0[i];
      size = dd_cell_size[i] + radd;
      /* Check if the cell boundary in this direction is
       * perpendicular to the Cartesian axis.
       */
      for(j=i+1; j<DIM; j++) {
	if (box[j][i] != 0) {
	  /* Correct the offset for the home cell location */
	  grid->cell_offset[i] += dd->cell_ns_x0[j]*box[j][i]/box[j][j];
	  /* Correct the offset and size for the off-diagonal
	   * displacement of opposing DD cell corners.
	   */
	  /* Without rouding we would need to add box[j][i]*radd/box[j][j]); */
	  add_tric = dd_cell_size[j]*box[j][i]/box[j][j];
	  if (box[j][i] < 0) {
	    grid->cell_offset[i] += add_tric;
	    size -= add_tric;
	  } else {
	    size += add_tric;
	  }
	}
      }
    }
    if (!bDDRect) {
      /* No DD or the box is triclinic is this direction:
       * we will use the normal grid ns that checks all cells
       * that are within cut-off distance of the i-particle.
       */
      if (ePBC==epbcXY && i==ZZ)
	grid->n[i] = 1;
      else
	grid->n[i] = delta*size/rlist;
      grid->cell_size[i] = size/grid->n[i];
      grid->ncpddc[i] = grid->n[i];
    } else {
      /* We use grid->ncpddc[i] != nc[i] such that all particles
       * in one ns cell belong to one DD cell only.
       * We can then beforehand exclude certain ns grid cells
       * for non-home i-particles.
       */
      size = dd_cell_size[i];
      grid->ncpddc[i] = delta*size/rlist;
      grid->cell_size[i] = size/grid->ncpddc[i];
      grid->n[i] = grid->ncpddc[i] + (int)(rlist/grid->cell_size[i]) + 1;
    }
    if (debug)
      fprintf(debug,"grid dim %d size %d x %f: %f - %f\n",
	      i,grid->n[i],grid->cell_size[i],
	      grid->cell_offset[i],
	      grid->cell_offset[i]+grid->n[i]*grid->cell_size[i]);
  }
}

void init_grid(FILE *fplog,t_grid *grid,int delta,gmx_domdec_t *dd,
	       int ePBC,matrix box,real rlistlong,int ncg)
{
  int     d,m;
  
  set_grid_sizes(fplog,ePBC,box,rlistlong,delta,dd,grid);

  fprintf(fplog,"Grid: %d x %d x %d cells\n",
	  grid->n[XX],grid->n[YY],grid->n[ZZ]);

  grid->nr      = ncg;
  grid->ncells  = grid->n[XX]*grid->n[YY]*grid->n[ZZ];
  grid->maxcells= 2*grid->ncells;
  grid->delta	= delta;
  grid->nr_alloc= over_alloc(grid->nr)+1;
  snew(grid->cell_index,grid->nr_alloc);
  snew(grid->a,grid->nr_alloc);
  snew(grid->index,grid->maxcells);
  snew(grid->nra,grid->maxcells);
    
  if (debug) 
    fprintf(debug,"Succesfully allocated memory for grid pointers.");
}

void done_grid(t_grid *grid)
{
  grid->nr      = 0;
  clear_ivec(grid->n);
  grid->ncells  = 0;
  grid->maxcells= 0;
  grid->delta	= 0;
  sfree(grid->cell_index);
  sfree(grid->a);
  sfree(grid->index);
  sfree(grid->nra);
  
  if (debug) 
    fprintf(debug,"Succesfully freed memory for grid pointers.");
}

int xyz2ci_(int nry,int nrz,int x,int y,int z)
/* Return the cell index */
{
  return (nry*nrz*x+nrz*y+z);
}

void ci2xyz(t_grid *grid, int i, int *x, int *y, int *z)
/* Return x,y and z from the cell index */
{
  int ci;

  range_check(i,0,grid->nr);

  ci = grid->cell_index[i];
  *x  = ci / (grid->n[YY]*grid->n[ZZ]);
  ci -= (*x)*grid->n[YY]*grid->n[ZZ];
  *y  = ci / grid->n[ZZ];
  ci -= (*y)*grid->n[ZZ];
  *z  = ci;
}

static int ci_not_used(ivec n)
{
  /* Return an improbable value */
  return xyz2ci(n[YY],n[ZZ],3*n[XX],3*n[YY],3*n[ZZ]);
}

void grid_first(FILE *fplog,t_grid *grid,gmx_domdec_t *dd,
		int ePBC,matrix box,real rlistlong,int ncg)
{
  int    i,k,m,ncells;
  ivec   cx;

  /* Must do this every step because other routines may override it. */
  init_range_check();

  set_grid_sizes(fplog,ePBC,box,rlistlong,grid->delta,dd,grid);

  ncells = grid->n[XX]*grid->n[YY]*grid->n[ZZ];

  if (grid->ncells != ncells) {
    if (ncells > grid->maxcells) { 
      srenew(grid->nra,ncells);
      srenew(grid->index,ncells);
      for(i=grid->maxcells; (i<ncells); i++) {
	grid->nra[i] = 0;
	grid->index[i] = 0;
      }
      fprintf(fplog,"WARNING: your box is exploding! (ncells = %d)\n",ncells);
      grid->maxcells = ncells;
    }
    grid->ncells = ncells;
  }

  grid->nr = ncg;
  if (grid->nr+1 > grid->nr_alloc) {
    grid->nr_alloc = over_alloc(grid->nr) + 1;
    srenew(grid->cell_index,grid->nr_alloc);
    srenew(grid->a,grid->nr_alloc);
  }
  
  for(i=0; (i<ncells); i++)
    grid->nra[i]=0;
  for(i=0; i<grid->maxcells; i++)
    grid->index[i] = 10000*random();
}

static void calc_bor(int cg0,int cg1,int ncg,int CG0[2],int CG1[2])
{
  if (cg1 > ncg) {
    CG0[0]=cg0;
    CG1[0]=ncg;
    CG0[1]=0;
    CG1[1]=cg1-ncg;
  }
  else {
    CG0[0]=cg0;
    CG1[0]=cg1;
    CG0[1]=0;
    CG1[1]=0;
  }
  if (debug) {
    int m;
    
    fprintf(debug,"calc_bor: cg0=%d, cg1=%d, ncg=%d\n",cg0,cg1,ncg);
    for(m=0; (m<2); m++)
      fprintf(debug,"CG0[%d]=%d, CG1[%d]=%d\n",m,CG0[m],m,CG1[m]);
  }

}

void calc_elemnr(FILE *fplog,t_grid *grid,int cg0,int cg1,int ncg)
{
  int    CG0[2],CG1[2];
  int    *cell_index=grid->cell_index;
  int    *nra=grid->nra;
  int    i,m,ncells;
  int    ci,not_used;

  ncells=grid->ncells;
  if(ncells<=0) 
    gmx_fatal(FARGS,"Number of grid cells is zero. Probably the system and box collapsed.\n");

  not_used = ci_not_used(grid->n);

  calc_bor(cg0,cg1,ncg,CG0,CG1);
  for(m=0; (m<2); m++)
    for(i=CG0[m]; (i<CG1[m]); i++) {
      ci = cell_index[i];
      if (ci != not_used) {
	range_check(ci,0,ncells);
	nra[ci]++;
      }
    }
}

void calc_ptrs(t_grid *grid)
{
  int *index = grid->index;
  int *nra   = grid->nra;
  int ix,iy,iz,ci,nr;
  int nnra,ncells;

  ncells=grid->ncells;
  if(ncells<=0) 
    gmx_fatal(FARGS,"Number of grid cells is zero. Probably the system and box collapsed.\n");
  
  ci=nr=0;
  for(ix=0; (ix < grid->n[XX]); ix++)
    for(iy=0; (iy < grid->n[YY]); iy++) 
      for(iz=0; (iz < grid->n[ZZ]); iz++,ci++) {
	range_check(ci,0,ncells);
	index[ci] = nr;
	nnra      = nra[ci];
	nr       += nnra;
	nra[ci]   = 0;
      }
}

void grid_last(FILE *log,t_grid *grid,int cg0,int cg1,int ncg)
{
  int    CG0[2],CG1[2];
  int    i,m;
  int    ci,not_used,ind,ncells;
  int    *cell_index = grid->cell_index;
  int    *nra        = grid->nra;
  int    *index      = grid->index;
  int    *a          = grid->a;

  ncells=grid->ncells;
  if (ncells <= 0) 
    gmx_fatal(FARGS,"Number of grid cells is zero. Probably the system and box collapsed.\n");

  not_used = ci_not_used(grid->n);

  calc_bor(cg0,cg1,ncg,CG0,CG1);
  for(m=0; (m<2); m++)
    for(i=CG0[m]; (i<CG1[m]); i++) {
      ci     = cell_index[i];
      if (ci != not_used) {
	range_check(ci,0,ncells);
	ind    = index[ci]+nra[ci]++;
	range_check(ind,0,grid->nr);
	a[ind] = i;
      }
    }
}

void fill_grid(FILE *log,
	       gmx_domdec_t *dd,
	       t_grid *grid,matrix box,
	       int cg0,int cg1,rvec cg_cm[])
{
  int    *cell_index=grid->cell_index;
  int    nrx,nry,nrz;
  rvec   n_box,offset;
  int    cell,cg,d,not_used;
  ivec   shift0,b0,b1,ind;
  bool   bUse;
  
  /* Initiate cell borders */
  nrx = grid->n[XX];
  nry = grid->n[YY];
  nrz = grid->n[ZZ];
  for(d=0; d<DIM; d++) {
    if (grid->cell_size[d] > 0)
      n_box[d] = 1/grid->cell_size[d];
    else
      n_box[d] = 0;
  }

  if (debug)
    fprintf(debug,"Filling grid from %d to %d\n",cg0,cg1);

  /* We assume here that the charge group center of mass is always
   * 0 <= cgcm < box
   * If not this will generate errors (SEGV). If you suspect this, turn on
   * DEBUG_PBC
   */
  debug_gmx();
  if (dd == NULL) {
    for (cg=cg0; cg<cg1; cg++) {
      for(d=0; d<DIM; d++) {
	ind[d] = cg_cm[cg][d]*n_box[d];
	if (ind[d] == grid->n[d])
	  ind[d]--;
      }
#ifdef DEBUG_PBC
#define myrc(ixyz,n) if ((ixyz<0) || (ixyz>=n)) gmx_fatal(FARGS,"%s=%d(max=%d), index=%d, i=%d, cgcm=(%f,%f,%f)",#ixyz,ixyz,n,index,cg,cg_cm[index][XX],cg_cm[index][YY],cg_cm[index][ZZ])
      myrc(ix,nrx);
      myrc(iy,nry);
      myrc(iz,nrz);
#undef myrc
#endif
      cell_index[cg] = xyz2ci(nry,nrz,ind[XX],ind[YY],ind[ZZ]);
    }
  } else {
    copy_rvec(grid->cell_offset,offset);
    for(cell=0; cell<dd->ncell; cell++) {
      /* Determine the ns grid cell limits for this DD cell */
      for(d=0; d<DIM; d++) {
	shift0[d] = dd->shift[cell][d];
	if (grid->ncpddc[d] == grid->n[d]) {
	  b0[d] = 0;
	  b1[d] = grid->n[d];
	} else {
	  b0[d] = shift0[d];
	  b1[d] = b0[d] + 1;
	  b0[d] *= grid->ncpddc[d];
	  b1[d] *= grid->ncpddc[d];
	  if (b1[d] > grid->n[d])
	    b1[d] = grid->n[d];
	}
      }

      not_used = ci_not_used(grid->n);

      /* Put all the charge groups of this DD cell on the grid */
      cg0 = dd->ncg_cell[cell];
      cg1 = dd->ncg_cell[cell+1];
      for(cg=cg0; cg<cg1; cg++) {
	bUse = TRUE;
	for(d=0; d<DIM; d++) {
	  ind[d] = (cg_cm[cg][d] - offset[d])*n_box[d];
	  /* Here we have to correct for rounding problems,
	   * as this cg_cm to cell index operation is not necessarily
	   * binary identical to the operation for the DD cell assignment
	   * and therefore a cg could end up in an unused grid cell.
	   */
	  if (ind[d] < b0[d]) {
	    ind[d]++;
	  } else if (ind[d] >= b1[d]) {
	    if (shift0[d] == 0) {
	      ind[d]--;
	    } else {
	      /* Charge groups in this DD cell further away than the cut-off
	       * in direction do not participate in non-bonded interactions.
	       */
	      bUse = FALSE;
	    }
	  }
	}
	if (bUse)
	  cell_index[cg] = xyz2ci(nry,nrz,ind[XX],ind[YY],ind[ZZ]);
	else
	  cell_index[cg] = not_used;
      }
    }
  }
  debug_gmx();
}

void check_grid(FILE *log,t_grid *grid)
{
  int ix,iy,iz,ci,cci,nra;

  if(grid->ncells<=0) 
    gmx_fatal(FARGS,"Number of grid cells is zero. Probably the system and box collapsed.\n");
  
  ci=0;
  cci=0;
  for(ix=0; (ix<grid->n[XX]); ix++)
    for(iy=0; (iy<grid->n[YY]); iy++)
      for(iz=0; (iz<grid->n[ZZ]); iz++,ci++) {
	if (ci > 0) {
	  nra=grid->index[ci]-grid->index[cci];
	  if (nra != grid->nra[cci]) 
	    gmx_fatal(FARGS,"nra=%d, grid->nra=%d, cci=%d",
			nra,grid->nra[cci],cci);
	}
	cci=xyz2ci(grid->n[YY],grid->n[ZZ],ix,iy,iz);
	range_check(cci,0,grid->ncells);
	
	if (cci != ci) 
	  gmx_fatal(FARGS,"ci = %d, cci = %d",ci,cci);
      }
}

void print_grid(FILE *log,t_grid *grid)
{
  int i,nra,index;
  int ix,iy,iz,ci;

  fprintf(log,"nr:    %d\n",grid->nr);
  fprintf(log,"nrx:   %d\n",grid->n[XX]);
  fprintf(log,"nry:   %d\n",grid->n[YY]);
  fprintf(log,"nrz:   %d\n",grid->n[ZZ]);
  fprintf(log,"delta: %d\n",grid->delta);
  fprintf(log,"    i  cell_index\n");
  for(i=0; (i<grid->nr); i++)
    fprintf(log,"%5d  %5d\n",i,grid->cell_index[i]);
  fprintf(log,"cells\n");
  fprintf(log," ix iy iz   nr  index  cgs...\n");
  ci=0;
  for(ix=0; (ix<grid->n[XX]); ix++)
    for(iy=0; (iy<grid->n[YY]); iy++)
      for(iz=0; (iz<grid->n[ZZ]); iz++,ci++) {
	index=grid->index[ci];
	nra=grid->nra[ci];
	fprintf(log,"%3d%3d%3d%5d%5d",ix,iy,iz,nra,index);
	for(i=0; (i<nra); i++)
	  fprintf(log,"%5d",grid->a[index+i]);
	fprintf(log,"\n");
      }
  fflush(log);
}

void mv_grid(t_commrec *cr,t_grid *grid)
{
  int i,start,nr;
  int cur=cr->nodeid;
  int *ci,*cgindex;
#define next ((cur+1) % (cr->nnodes-cr->npmenodes))

  ci = grid->cell_index;
  cgindex = pd_cgindex(cr);
  for(i=0; (i<cr->nnodes-1); i++) {
    start = cgindex[cur];
    nr    = cgindex[cur+1] - start;
    gmx_tx(cr,cr->left,&(ci[start]),nr*sizeof(*ci));
    
    start = cgindex[next];
    nr    = cgindex[next+1] - start;
    gmx_rx(cr,cr->right,&(ci[start]),nr*sizeof(*ci));
    
    gmx_tx_wait(cr->left);
    gmx_rx_wait(cr->right);
    
    cur=next;
  }
}

