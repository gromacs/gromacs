/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#include <stdio.h>
#include "typedefs.h"
#include "futil.h"
#include "smalloc.h"
#include "grids.h"

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
  int i,j;

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
  int i,j;

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
  int  ix,iy,iz,n,ig;
  real x,y,z,g;

  n=1;
  fp=ffopen(fn,"w");  
  for(ix=0; (ix<nx); ix++) {
    for(iy=0; (iy<ny); iy++) {
      for(iz=0; (iz<nz); iz++) {
	g=grid[ix][iy][iz].re;
	ig=g;
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

