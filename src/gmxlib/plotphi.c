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
 * GROup of MAchos and Cynical Suckers
 */
#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "writeps.h"
#include "smalloc.h"
#include "plotphi.h"
#include "macros.h"
#include "futil.h"

static real rgbset(real col)
{
  int icol32;

  icol32=32.0*col;
  return icol32/32.0;
}

void plot_phi(char *fn,rvec box,int natoms,rvec x[],real phi[])
{
  FILE *eps;
  real phi_max,rr,gg,bb,fac,dx,x0,y0;
  real offset;
  int  i;
  
  phi_max=phi[0];
  rr=gg=bb=0.0;
  for(i=0; (i<natoms); i++) 
    phi_max=max(phi_max,fabs(phi[i]));
    
  if (phi_max==0.0) {
    fprintf(stderr,"All values zero, see .out file\n");
    return;
  }
  offset=20.0;
  fac=15.0;
#ifdef DEBUG
  fprintf(stderr,"Scaling box by %g\n",fac);
#endif
  eps=ps_open(fn,0,0,fac*box[XX]+2*offset,fac*box[YY]+2*offset);
  ps_translate(eps,offset,offset);
  ps_color(eps,0,0,0);
  ps_box(eps,1,1,fac*box[XX]-1,fac*box[YY]-1);
  dx=0.15*fac;
  for(i=0; (i<natoms); i++) {
    rr=gg=bb=1.0;
    if (phi[i] < 0)
      gg=bb=(1.0+(phi[i]/phi_max));
    else 
      rr=gg=(1.0-(phi[i]/phi_max));
    rr=rgbset(rr);
    gg=rgbset(gg);
    bb=rgbset(bb);
    ps_color(eps,rr,gg,bb);
    x0=fac*x[i][XX];
    y0=fac*x[i][YY];
    ps_fillbox(eps,x0-dx,y0-dx,x0+dx,y0+dx);
  }
  ps_close(eps);
}

void plot_qtab(char *fn,int nx,int ny,int nz,real ***qtab)
{
  rvec box;
  rvec *xx;
  real *phi;
  int  i,npt,ix,iy,iz;
  
  box[XX]=nx;
  box[YY]=ny;
  box[ZZ]=nz;

  npt=(box[XX]*box[YY]*box[ZZ]);
  snew(xx,npt);
  snew(phi,npt);
  nx/=2;
  ny/=2;
  nz/=2;
  i=0;
  for(ix=-nx; (ix<nx); ix++)
    for(iy=-ny; (iy<ny); iy++)
      for(iz=-nz; (iz<nz); iz++,i++) {
	xx[i][XX]=ix+nx+0.5;
	xx[i][YY]=iy+ny+0.5;
	xx[i][ZZ]=iz+nz+0.5; /* onzin */
	phi[i]=qtab[ix+nx][iy+ny][iz+nz];
      }
  
  plot_phi(fn,box,npt,xx,phi);
  
  sfree(xx);
  sfree(phi);
}

void print_phi(char *fn,int natoms,rvec x[],real phi[])
{
  FILE *fp;
  int  i;
  
  fp=ffopen(fn,"w");
  for(i=0; (i<natoms); i++)
    fprintf(fp,"%10d  %12.5e\n",i,phi[i]);
  fclose(fp);
}

void write_pqr(char *fn,t_atoms *atoms,rvec x[],real phi[],real dx)
{
  FILE *fp;
  int  i,rnr;
  
  fp=ffopen(fn,"w");
  for(i=0; (i<atoms->nr); i++) {
    rnr=atoms->atom[i].resnr;
    fprintf(fp,"%-6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
	    "ATOM",i+1,*atoms->atomname[i],*atoms->resname[rnr],' ',rnr+1,
	    10*(dx+x[i][XX]),10*x[i][YY],10*(x[i][ZZ]),0,phi[i]);
  }
  fclose(fp);
}

