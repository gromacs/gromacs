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

#include <stdlib.h>
#include <stdio.h>

#include "typedefs.h"
#include "futil.h"
#include "smalloc.h"
#include "futil.h"
#include "macros.h"
#include "network.h"
#include "fftgrid.h"
#include "gmx_fft.h"
#include "gmx_parallel_3dfft.h"


#ifdef GMX_MPI
static void print_parfft(FILE *fp,char *title,t_parfft *pfft)
{
  fprintf(fp,"PARALLEL FFT DATA:\n"
	  "   local_nx:                 %3d  local_x_start:                 %3d\n"
	  "   local_ny_after_transpose: %3d  local_y_start_after_transpose  %3d\n",
	  pfft->local_nx,pfft->local_x_start,pfft->local_ny_after_transpose,
	  pfft->local_y_start_after_transpose);
}
#endif


void *
gmx_alloc_aligned(size_t size)
{
    void *p0,*p;
    
    p0 = malloc(size+32);
    
    if(p0 == NULL)
    {
        gmx_fatal(FARGS,"Failed to allocated %u bytes of aligned memory.",size+32);
    }
    
    p = (void *) (((size_t) p0 + 32) & (~((size_t) 31)));
    
    /* Yeah, yeah, we cannot free this pointer, but who cares... */
    return p;
}

t_fftgrid *mk_fftgrid(int          nx,
                      int          ny,
                      int          nz,
                      int          *node2slab,
		      int          *slab2grid_x,
                      t_commrec *  cr,
                      bool         bReproducible)
{
/* parallel runs with non-parallel ffts haven't been tested yet */
  int           nnodes;
  int           localsize;
  t_fftgrid *   grid;
  int           flags;
    
  nnodes = 1;
#ifdef GMX_MPI
  if (cr && cr->nnodes > 1) {
    MPI_Comm_size(cr->mpi_comm_mygroup,&nnodes);
  }
#endif
 
  snew(grid,1);
  grid->nx   = nx;
  grid->ny   = ny;
  grid->nz   = nz;
  grid->nxyz = nx*ny*nz;
  grid->bParallel = (nnodes > 1);
  
  if (grid->bParallel)
  {
      grid->la2r = (nz/2+1)*2;
  }
  else
  {
      grid->la2r = nz;  
  }
  
  grid->la2c = (nz/2+1);    

  grid->la12r = ny*grid->la2r;
  
  if (grid->bParallel)
  {
      grid->la12c = nx*grid->la2c;
  }
  else
  {
    grid->la12c = ny*grid->la2c;
  }
  
  grid->nptr = nx*ny*grid->la2c*2;

  if (grid->bParallel) 
  {
#ifdef GMX_MPI
      gmx_parallel_3dfft_init(&grid->mpi_fft_setup,nx,ny,nz,
			      node2slab,slab2grid_x,cr->mpi_comm_mygroup,
			      bReproducible);
          
      gmx_parallel_3dfft_limits(grid->mpi_fft_setup,
                                &(grid->pfft.local_x_start),                                
                                &(grid->pfft.local_nx),
                                &(grid->pfft.local_y_start_after_transpose),
                                &(grid->pfft.local_ny_after_transpose));
#else
    gmx_fatal(FARGS,"Parallel FFT supported with MPI only!");
#endif
  }
  else 
  {
      flags = bReproducible ? GMX_FFT_FLAG_CONSERVATIVE : 0;
      gmx_fft_init_3d_real(&grid->fft_setup,nx,ny,nz,flags);
  }
  grid->ptr = gmx_alloc_aligned(grid->nptr*sizeof(*(grid->ptr)));
  
#ifdef GMX_MPI
  if (grid->bParallel && debug) 
  {
    print_parfft(debug,"Plan", &grid->pfft);
  }
  if (grid->bParallel)
  {
      int maxlocalsize;
      maxlocalsize = max(grid->la2c *2*grid->pfft.local_nx*ny,
			 grid->la12c*2*grid->pfft.local_ny_after_transpose);
      grid->workspace = gmx_alloc_aligned(maxlocalsize*sizeof(*(grid->workspace)));
  }
  else
  {
      grid->workspace = gmx_alloc_aligned(grid->nptr*sizeof(*(grid->workspace)));
  }
#else /* no MPI */
  grid->workspace = gmx_alloc_aligned(grid->nptr*sizeof(*(grid->workspace)));
#endif

  return grid;
}

void 
pr_fftgrid(FILE *fp,char *title,t_fftgrid *grid)
{
  int     i,j,k,index_x,index_xy,ntot=0;
  int     nx,ny,nz,nx2,ny2,nz2,la12,la2;
  real *  ptr;

  /* Unpack structure */
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  for(i=0; (i<nx); i++) {
    index_x = la12*i;
    for(j=0; (j<ny); j++) {
      index_xy = index_x + la2*j;
      for(k=0; (k<nz); k++) {
	if (ptr[index_xy+k] != 0) {
	  fprintf(fp,"%-12s  %5d  %5d  %5d  %12.5e\n",
		  title,i,j,k,ptr[index_xy+k]);
	  ntot++;
	}
      }
    }
  }
  fprintf(fp,"%d non zero elements in %s\n",ntot,title);
}

void done_fftgrid(t_fftgrid *grid)
{
//memory can't be freed because it is allocated by gmx_alloc_aligned
  if (grid->ptr) {
    //sfree(grid->ptr);
    grid->ptr = NULL;
  }
 
  if (grid->workspace) {
      //sfree(grid->workspace);
      grid->workspace=NULL;
  }
}


void gmxfft3D(t_fftgrid *grid,int dir,t_commrec *cr)
{
  void *tmp;

  if (grid->bParallel) 
  {
#ifdef GMX_MPI
    if( dir == GMX_FFT_REAL_TO_COMPLEX || dir == GMX_FFT_COMPLEX_TO_REAL )
    {
      gmx_parallel_3dfft(grid->mpi_fft_setup,dir,grid->ptr,grid->ptr);
    }
    else
    {
        gmx_fatal(FARGS,"Invalid direction for FFT: %d",dir);
    }
#else
    gmx_fatal(FARGS,"Parallel FFT supported with MPI only!");   
#endif
  }
  else
  {
    if( dir == GMX_FFT_REAL_TO_COMPLEX || dir == GMX_FFT_COMPLEX_TO_REAL)
    {
        gmx_fft_3d_real(grid->fft_setup,dir,grid->ptr,grid->workspace);
        tmp = grid->ptr;
        grid->ptr = grid->workspace;
        grid->workspace = tmp;        
    }
    else
    {
      gmx_fatal(FARGS,"Invalid direction for FFT: %d",dir);
    }
  }
}

void clear_fftgrid(t_fftgrid *grid)
{
    /* clears the whole grid */
  int      i,ngrid;
  real *   ptr;
  
  ngrid = grid->nptr;
  ptr   = grid->ptr;
  
  for (i=0; (i<ngrid); i++) {
    ptr[i] = 0;
  }
}

void unpack_fftgrid(t_fftgrid *grid,int *nx,int *ny,int *nz,
		    int *nx2,int *ny2,int *nz2,
		    int *la2,int *la12,bool bReal,real **ptr)
{
  *nx  = grid->nx;
  *ny  = grid->ny;
  *nz  = grid->nz;
  *nx2 = 2*grid->nx;
  *ny2 = 2*grid->ny;
  *nz2 = 2*grid->nz;
  if(bReal) {
    *la2 = grid->la2r;
    *la12= grid->la12r;
  } 
  else {
    *la2 = grid->la2c;
    *la12= grid->la12c;
  }
  *ptr = grid->ptr;
}



/*****************************************************************
 * 
 * For backward compatibility (for testing the ewald code vs. PPPM etc)
 * some old grid routines are retained here.
 *
 ************************************************************************/

real ***mk_rgrid(int nx,int ny,int nz)
{
  real *ptr1;
  real **ptr2;
  real ***ptr3;
  int  i,j,n2,n3;
  
  snew(ptr1,nx*ny*nz);
  snew(ptr2,nx*ny);
  snew(ptr3,nx);
  
  n2=n3=0;
  for(i=0; (i<nx); i++) {
    ptr3[i]=&(ptr2[n2]);
    for(j=0; (j<ny); j++,n2++) { 
      ptr2[n2] = &(ptr1[n3]);
      n3 += nz;
    }
  }
  return ptr3;
}

void free_rgrid(real ***grid,int nx,int ny)
{
  int i;

  sfree(grid[0][0]);  
  for(i=0; (i<nx); i++) {
    sfree(grid[i]);
  }
  sfree(grid);
}

real print_rgrid(FILE *fp,char *title,int nx,int ny,int nz,real ***grid)
{
  int  ix,iy,iz;
  real g,gtot;
  
  gtot=0;
  if (fp)
    fprintf(fp,"Printing all non-zero real elements of %s\n",title);
  for(ix=0; (ix<nx); ix++)
    for(iy=0; (iy<ny); iy++)
      for(iz=0; (iz<nz); iz++) {
	g=grid[ix][iy][iz];
	if (fp && (g != 0))
	  fprintf(fp,"%s[%2d][%2d][%2d] = %12.5e\n",title,ix,iy,iz,g);
	gtot+=g;
      }
  return gtot;
}

void print_rgrid_pdb(char *fn,int nx,int ny,int nz,real ***grid)
{
  FILE *fp;
  int  ix,iy,iz,n,ig;
  real x,y,z,g;

  n=1;
  fp=ffopen(fn,"w");  
  for(ix=0; (ix<nx); ix++) {
    for(iy=0; (iy<ny); iy++) {
      for(iz=0; (iz<nz); iz++) {
	g=grid[ix][iy][iz];
	ig=g;
	if ((ig != 0) || (1)) {
	  x = 4*ix;
	  y = 4*iy;
	  z = 4*iz;
	  fprintf(fp,"ATOM  %5d  Na   Na     1    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
		  n++,x,y,z,0.0,g);
	}
      }
    }
  }
  fclose(fp);
}

void clear_rgrid(int nx,int ny,int nz,real ***grid)
{
  int i,j,k;
  
  for(i=0; (i<nx); i++)
    for(j=0; (j<ny); j++)
      for(k=0; (k<nz); k++)
	grid[i][j][k] = 0;
}

void clear_cgrid(int nx,int ny,int nz,t_complex ***grid)
{
  int i,j,k;
  
  for(i=0; (i<nx); i++)
    for(j=0; (j<ny); j++)
      for(k=0; (k<nz); k++)
	grid[i][j][k] = cnul;
}

t_complex ***mk_cgrid(int nx,int ny,int nz)
{
  t_complex *ptr1;
  t_complex **ptr2;
  t_complex ***ptr3;
  int  i,j,n2,n3;
  
  snew(ptr1,nx*ny*nz);
  snew(ptr2,nx*ny);
  snew(ptr3,nx);
  
  n2=n3=0;
  for(i=0; (i<nx); i++) {
    ptr3[i]=&(ptr2[n2]);
    for(j=0; (j<ny); j++,n2++) { 
      ptr2[n2] = &(ptr1[n3]);
      n3 += nz;
    }
  }
  return ptr3;
}

void free_cgrid(t_complex ***grid,int nx,int ny)
{
  int i;

  sfree(grid[0][0]);
  for(i=0; (i<nx); i++) 
    sfree(grid[i]);
  sfree(grid);
}

t_complex print_cgrid(FILE *fp,char *title,int nx,int ny,int nz,
		      t_complex ***grid)
{
  int     ix,iy,iz;
  t_complex g,gtot;
  
  gtot=cnul;
  if (fp)
    fprintf(fp,"Printing all non-zero complex elements of %s\n",title);
  for(ix=0; (ix<nx); ix++)
    for(iy=0; (iy<ny); iy++)
      for(iz=0; (iz<nz); iz++) {
	g=grid[ix][iy][iz];
	if (fp  && ((g.re != 0) || (g.im != 0)))
	  fprintf(fp,"%s[%2d][%2d][%2d] = %12.5e + i %12.5e\n",
		  title,ix,iy,iz,g.re,g.im);
	gtot = cadd(gtot,g);
      }
  return gtot;
}

void print_cgrid_pdb(char *fn,int nx,int ny,int nz,t_complex ***grid)
{
  FILE *fp;
  int  ix,iy,iz,n;
  real x,y,z,g;

  n=1;
  fp=ffopen(fn,"w");  
  for(ix=0; (ix<nx); ix++) {
    for(iy=0; (iy<ny); iy++) {
      for(iz=0; (iz<nz); iz++) {
	g=grid[ix][iy][iz].re;
	if (g != 0) {
	  x = 4*ix;
	  y = 4*iy;
	  z = 4*iz;
	  fprintf(fp,"ATOM  %5d  Na   Na     1    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
		  n++,x,y,z,0.0,g);
	}
      }
    }
  }
  fclose(fp);
}








