/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_lrutil_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "lrutil.h"
#include "smalloc.h"
#include "physics.h"
#include "futil.h"
#include "fftgrid.h"
#include "writeps.h"
#include "macros.h"

#define p2(x) ((x)*(x))
#define p3(x) ((x)*(x)*(x)) 
#define p4(x) ((x)*(x)*(x)*(x)) 

static real A,A_3,B,B_4,C,c1,c2,c3,c4,c5,c6,One_4pi,FourPi_V,N0;

void set_LRconsts(FILE *log,real r1,real rc,rvec box,t_forcerec *fr)
{
  if (r1 < rc) {
    A   = (2*r1-5*rc)/(p3(rc)*p2(rc-r1));
    B   = (4*rc-2*r1)/(p3(rc)*p3(rc-r1));
    /*C   = (10*rc*rc-5*rc*r1+r1*r1)/(6*rc*rc); Hermans Eq. not correct */
  }
  else {
    fatal_error(0,"r1 (%f) >= rc (%f) in %s, line %d",
		r1,rc,__FILE__,__LINE__);
  }
  A_3 = A/3.0;
  B_4 = B/4.0;
  C   = 1/rc-A_3*p3(rc-r1)-B_4*p4(rc-r1);
  N0  = 2.0*M_PI*p3(rc)*p3(rc-r1);

  if (fr) {
    fr->A   = A;
    fr->B   = B;
    fr->A_3 = A_3;
    fr->B_4 = B_4;
    fr->C   = C;
  }
  
  FourPi_V=4.0*M_PI/(box[XX]*box[YY]*box[ZZ]);

  fprintf(log,"Constants for short-range and fourier stuff:\n"
	  "r1 = %10.3f,  rc = %10.3f\n"
	  "A  = %10.3e,  B  = %10.3e,  C  = %10.3e, FourPi_V = %10.3e\n",
	  r1,rc,A,B,C,FourPi_V);

  /* Constants derived by Mathematica */
  c1 = -40*rc*rc    + 50*rc*r1    - 16*r1*r1;
  c2 =  60*rc       - 30*r1;
  c3 = -10*rc*rc*rc + 20*rc*rc*r1 - 13*rc*r1*r1 + 3*r1*r1*r1;
  c4 = -20*rc*rc    + 40*rc*r1    - 14*r1*r1;
  c5 = -c2;
  c6 = -5*rc*rc*r1  +  7*rc*r1*r1 - 2*r1*r1*r1;

  fprintf(log,"c1 = %10.3e,  c2 = %10.3e,  c3 = %10.3e\n"
	  "c4 = %10.3e,  c5 = %10.3e,  c6 = %10.3e,  N0 = %10.3e\n",
	  c1,c2,c3,c4,c5,c6,N0);
    
  One_4pi = 1.0/(4.0*M_PI);
}

real gk(real k,real rc,real r1)
/* Spread function in Fourier space. Eqs. 56-64 */
{
  real N,gg;
  
  N  = N0*p4(k);

  /* c1 thru c6 consts are global variables! */  
  gg = (FourPi_V / (N*k)) * 
    ( c1*k*cos(k*rc) + (c2+c3*k*k)*sin(k*rc) + 
      c4*k*cos(k*r1) + (c5+c6*k*k)*sin(k*r1) );
  
  return gg;
}

real calc_dx2(rvec xi,rvec xj,rvec box)
{
  int  m;
  real dx,r2;
  
  r2=0;
  for(m=0; (m<DIM); m++) {
    dx=xj[m]-xi[m];
    if (dx < -0.5*box[m])
      dx+=box[m];
    else if (dx >= 0.5*box[m])
      dx-=box[m];
    r2+=dx*dx;
  }
  return r2;
}

void calc_dx(rvec xi,rvec xj,rvec box,rvec dx)
{
  int  m;
  real ddx;
  
  for(m=0; (m<DIM); m++) {
    ddx=xj[m]-xi[m];
    if (ddx < -0.5*box[m])
      ddx+=box[m];
    else if (ddx >= 0.5*box[m])
      ddx-=box[m];
    dx[m]=ddx;
  }
}

real phi_sr(int nj,rvec x[],real charge[],real rc,real r1,rvec box,
	    real phi[],t_block *excl)
{
  int  i,j,k,ni,i1,i2;
  real pp,r2,R,R_1,r12,rc2;
  real qi,qj,vsr,eps;
  
  vsr = 0.0;
  eps = ONE_4PI_EPS0;
  rc2 = rc*rc;
  r12 = r1*r1;
  ni=0;
  for(i=0; (i<nj-1); i++) {
    qi=charge[i];
    for(j=i+1; (j<nj); j++) {
      i1=excl->index[i];
      i2=excl->index[i+1];
      for(k=i1; (k<i2); k++)
	if (excl->a[k] == j)
	  break;
      if (k == i2) {
	r2=calc_dx2(x[i],x[j],box);
	if (r2 < rc2) {
	  qj  = charge[j];
	  R_1 = invsqrt(r2);
	  pp  = R_1-C;
	  if (r2 > r12) {
	    R   = 1.0/R_1;
	    pp -= A_3*p3(R-r1)+B_4*p4(R-r1);
	  }
	  phi[i] += eps*qj*pp;
	  phi[j] += eps*qi*pp;
	  vsr    += eps*qj*qi*pp;
	  ni++;
	}
      }
    }
  }
  fprintf(stderr,"There were %d short range interactions, vsr=%g\n",ni,vsr);
  
  return vsr;
}

real shiftfunction(real r1,real rc,real R)
{
  real dr;

  if (R <= r1)
    return 0.0;
  else if (R >= rc)
    return -1.0/(R*R);
  
  dr=R-r1;
  
  return A*dr*dr+B*dr*dr*dr;
}

real spreadfunction(real r1,real rc,real R)
{
  real dr;

  if (R <= r1)
    return 0.0;
  else if (R >= rc)
    return 0.0;
  
  dr=R-r1;
  
  return -One_4pi*(dr/R)*(2*A*(2*R-r1)+B*dr*(5*R-2*r1));
}

real potential(real r1,real rc,real R)
{
  if (R < r1)
    return (1.0/R-C);
  else if (R <= rc)
    return (1.0/R - A_3*p3(R-r1) - B_4*p4(R-r1) - C);
  else 
    return 0.0;
}

real calc_selfenergy(FILE *fp,int natoms,real charge[],t_block *excl)
{
  static bool bFirst=TRUE;
  static real Vself;
  int  i,i1,i2,j,k;
  real qq,qi,Vex,Vc;
  
  if (bFirst) {
    qq   = 0;
    for(i=0; (i<natoms); i++) 
      qq  += charge[i]*charge[i];
    Vc = 0.5*C*ONE_4PI_EPS0*qq;
    fprintf(fp,"calc_self: qq = %g, Vc=%g\n",qq,Vc);
  
    qq = 0;
    for(i=0; (i<excl->nr); i++) {
      i1 = excl->index[i];
      i2 = excl->index[i+1];
      qi = charge[i];
      for(j=i1; (j<i2); j++) {
	k = excl->a[j];
	if (k != i)
	  qq+=qi*charge[k];
      }
    }
    Vex   = qq*0.5*C*ONE_4PI_EPS0;
    Vself = Vc + Vex;
    
    fprintf(fp,"calc_self: qq = %g, Vex=%g, Vself=%g\n",qq,Vex,Vself);
    
    bFirst = FALSE;
  }
  
  return Vself;
}


void calc_ener(FILE *fp,char *title,bool bHeader,int nmol,
	       int natoms,real phi[],real charge[],t_block *excl)
{
  int  i,i1,i2,j,k;
  real qq,qi,vv,V,Vex,Vc,Vt;
  
  qq   = 0;
  V    = 0;
  for(i=0; (i<natoms); i++) {
    vv   = charge[i]*phi[i];
    V   += vv;
    qq  += charge[i]*charge[i];
  }
  V  = 0.5*V;
  Vc = 0.5*C*ONE_4PI_EPS0*qq;
  
  qq = 0;
  for(i=0; (i<excl->nr); i++) {
    i1 = excl->index[i];
    i2 = excl->index[i+1];
    qi = charge[i];
    for(j=i1; (j<i2); j++) {
      k = excl->a[j];
      if (k != i)
	qq+=qi*charge[k];
    }
  }
  Vex = qq*0.5*C*ONE_4PI_EPS0;
  
  Vt = V - Vc - Vex;
  
  if (bHeader) {
    fprintf(fp,"%12s  %12s  %12s  %12s  %12s  %12s\n",
	    "","Vphi","Vself","Vexcl","Vtot","1Mol");
    
  }
  fprintf(fp,"%12s  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n",
	  title,V,Vc,Vex,Vt,Vt/nmol);
}

real phi_aver(int natoms,real phi[])
{
  real phitot;
  int  i;

  phitot=0;  
  for(i=0; (i<natoms); i++)
    phitot=phitot+phi[i];
    
  return (phitot/natoms);
}

real symmetrize_phi(FILE *log,int natoms,real phi[],bool bVerbose)
{
  real phitot;
  int  i;

  phitot=phi_aver(natoms,phi); 
  if (bVerbose)
    fprintf(log,"phi_aver = %10.3e\n",phitot);
  
  for(i=0; (i<natoms); i++)
    phi[i]-=phitot;
    
  return phitot;
}

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
  eps=ps_open(fn,0,0,
	      (real)(fac*box[XX]+2*offset),(real)(fac*box[YY]+2*offset));
  ps_translate(eps,offset,offset);
  ps_color(eps,0,0,0);
  ps_box(eps,1,1,(real)(fac*box[XX]-1),(real)(fac*box[YY]-1));
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
    ps_fillbox(eps,(real)(x0-dx),(real)(y0-dx),(real)(x0+dx),(real)(y0+dx));
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
	    10*(dx+x[i][XX]),10*x[i][YY],10*(x[i][ZZ]),0.0,phi[i]);
  }
  fclose(fp);
}

