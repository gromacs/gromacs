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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "shift_util.h"
#include "smalloc.h"
#include "physics.h"
#include "txtdump.h"
#include "futil.h"
#include "names.h"
#include "fftgrid.h"
#include "writeps.h"
#include "macros.h"
#include "xvgr.h"

#define p2(x) ((x)*(x))
#define p3(x) ((x)*(x)*(x)) 
#define p4(x) ((x)*(x)*(x)*(x)) 

static real A,A_3,B,B_4,C,c1,c2,c3,c4,c5,c6,One_4pi,FourPi_V,Vol,N0;

void set_shift_consts(FILE *log,real r1,real rc,rvec box,t_forcerec *fr)
{
  /* A, B and C are recalculated in tables.c */
  if (r1 < rc) {
    A   = (2*r1-5*rc)/(p3(rc)*p2(rc-r1));
    B   = (4*rc-2*r1)/(p3(rc)*p3(rc-r1));
    /*C   = (10*rc*rc-5*rc*r1+r1*r1)/(6*rc*rc); Hermans Eq. not correct */
  }
  else
    fatal_error(0,"r1 (%f) >= rc (%f) in %s, line %d",
		r1,rc,__FILE__,__LINE__);

  A_3 = A/3.0;
  B_4 = B/4.0;
  C   = 1/rc-A_3*p3(rc-r1)-B_4*p4(rc-r1);
  N0  = 2.0*M_PI*p3(rc)*p3(rc-r1);

  Vol     =(box[XX]*box[YY]*box[ZZ]);
  FourPi_V=4.0*M_PI/Vol;

  if(debug)
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

  if(debug)
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

real gknew(real k,real rc,real r1)
{
  real rck,rck2;
  
  rck = rc*k;
  rck2= rck*rck;
  
  return -15.0*((rck2-3.0)*sin(rck) + 3*rck*cos(rck))/(Vol*rck2*rck2*rck);
}

real calc_dx2dx(rvec xi,rvec xj,rvec box,rvec dx)
{
  int  m;
  real ddx,dx2;
  
  dx2=0;
  for(m=0; (m<DIM); m++) {
    ddx=xj[m]-xi[m];
    if (ddx < -0.5*box[m])
      ddx+=box[m];
    else if (ddx >= 0.5*box[m])
      ddx-=box[m];
    dx[m]=ddx;
    dx2 += ddx*ddx;
  }
  return dx2;
}

real calc_dx2(rvec xi,rvec xj,rvec box)
{
  rvec dx;
  
  return calc_dx2dx(xi,xj,box,dx);
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

real new_f(real r,real rc)
{
  real rrc,rrc2,rrc3;
  
  rrc  = r/rc;
  rrc2 = rrc*rrc;
  rrc3 = rrc2*rrc;
  return 1.5*rrc2*rrc3 - 2.5*rrc3 + 1.0;
}

real new_phi(real r,real rc)
{
  real rrr;
  
  rrr = sqr(r/rc);
  
  return 1/r-(0.125/rc)*(15 + 3*rrr*rrr - 10*rrr);
}

real old_f(real r,real rc,real r1)
{
  real dr,r2;

  if (r <= r1)
    return 1.0;
  else if (r >= rc)
    return 0;
  
  dr = r-r1;
  r2 = r*r;
  return 1+A*r2*dr*dr+B*r2*dr*dr*dr;
}

real old_phi(real r,real rc,real r1)
{
  real dr;
  
  if (r <= r1)
    return 1/r-C;
  else if (r >= rc)
    return 0.0;
    
  dr = r-r1;
  
  return 1/r-A_3*dr*dr*dr-B_4*dr*dr*dr*dr - C;
}

real phi_sr(FILE *log,int nj,rvec x[],real charge[],real rc,real r1,rvec box,
	    real phi[],t_block *excl,rvec f_sr[],bool bOld)
{
  int  i,j,k,m,ni,i1,i2;
  real pp,r2,R,R_1,R_2,rc2;
  real qi,qj,vsr,eps,fscal;
  rvec dx;
  
  vsr = 0.0;
  eps = ONE_4PI_EPS0;
  rc2 = rc*rc;
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
	r2=calc_dx2dx(x[i],x[j],box,dx);
	if (r2 < rc2) {
	  qj    = charge[j];
	  R_1   = invsqrt(r2);
	  R_2   = R_1*R_1;
	  R     = invsqrt(R_2);
	  if (bOld) {
	    fscal = old_f(R,rc,r1)*R_2;
	    pp    = old_phi(R,rc,r1);
	  }
	  else {
	    fscal = new_f(R,rc)*R_2;
	    pp    = new_phi(R,rc);
	  }
	  phi[i] += eps*qj*pp;
	  phi[j] += eps*qi*pp;
	  vsr    += eps*qj*qi*pp;
	  for(m=0; (m<DIM); m++) {
	    f_sr[i][m] += dx[m]*fscal;
	    f_sr[j][m] -= dx[m]*fscal;
	  }
	  ni++;
	}
      }
    }
  }
  fprintf(log,"There were %d short range interactions, vsr=%g\n",ni,vsr);
  
  return vsr;
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



real shift_LRcorrection(FILE *fp,t_nsborder *nsb,t_commrec *cr,t_forcerec *fr,
			real charge[],t_block *excl,rvec x[],
			bool bOld,matrix box,matrix lr_vir)
{
  static bool bFirst=TRUE;
  static real Vself;
  int    i,i1,i2,j,k,m,iv,jv;
  int *AA;
  double qq; /* Necessary for precision */
  double isp=0.564189583547756;
  real   qi,dr,dr2,dr_1,dr_3,fscal,Vexcl,qtot=0;
  rvec   df,dx;
  real   r1=fr->rcoulomb_switch;
  real   rc=fr->rcoulomb;
  ivec   shift;     
  int    start=START(nsb);
  int    natoms=HOMENR(nsb);
  
  if (bFirst) {
    qq =0;  
    for(i=start; (i<start+natoms); i++) 
      qq  += charge[i]*charge[i];
    
    /* Obsolete, until we implement dipole and charge corrections.
       qtot=0;
       for(i=0;i<nsb->natoms;i++)
       qtot+=charge[i];
    */
   
    Vself = 0.5*C*ONE_4PI_EPS0*qq;
    if(debug) {
	fprintf(fp,"Long Range corrections for shifted interactions:\n");
	fprintf(fp,"r1 = %g, rc=%g\n",r1,rc);
	fprintf(fp,"start=%d,natoms=%d\n",start,natoms);
	fprintf(fp,"qq = %g, Vself=%g\n",qq,Vself);
    }
    
  }
  AA = excl->a;
  Vexcl = 0;
 
  for(i=start; (i<start+natoms); i++) {
    /* Initiate local variables (for this i-particle) to 0 */
    i1  = excl->index[i];
    i2  = excl->index[i+1];
    qi  = charge[i]*ONE_4PI_EPS0;

    /* Loop over excluded neighbours */
    for(j=i1; (j<i2); j++) {
      k = AA[j];
      /* 
       * First we must test whether k <> i, and then, because the
       * exclusions are all listed twice i->k and k->i we must select
       * just one of the two.
       * As a minor optimization we only compute forces when the charges
       * are non-zero.
       */
      if (k > i) {
	qq = qi*charge[k];
	if (qq != 0.0) {
	  dr2 = 0;
	  rvec_sub(x[i],x[k],dx);
	  for(m=DIM-1; m>=0; m--) {
	    if (dx[m] > 0.5*box[m][m])
	      rvec_dec(dx,box[m]);
	    else if (dx[m] < -0.5*box[m][m])
	      rvec_inc(dx,box[m]);
	    
	    dr2  += dx[m]*dx[m];
	  }
	  dr_1    = invsqrt(dr2);
	  dr      = 1.0/dr_1;
	  dr_3    = dr_1*dr_1*dr_1;
	  /* Compute exclusion energy and scalar force */

	  Vexcl  += qq*(dr_1-potential(r1,rc,dr));
	  fscal   = qq*(-shiftfunction(r1,rc,dr))*dr_3;
	  
	  if ((fscal != 0) && debug)
	    fprintf(debug,"i: %d, k: %d, dr: %.3f fscal: %.3f\n",i,k,dr,fscal);
	    
	  /* The force vector is obtained by multiplication with the 
	   * distance vector 
	   */
	  svmul(fscal,dx,df);
	  rvec_inc(fr->f_el_recip[k],df);
	  rvec_dec(fr->f_el_recip[i],df);
	  for(iv=0;iv<DIM;iv++)
	      for(jv=0;jv<DIM;jv++)
		  lr_vir[iv][jv]+=0.5*dx[iv]*df[jv];
	}
      }
    }
  }
  if (bFirst && debug)
    fprintf(fp,"Long Range correction: Vexcl=%g\n",Vexcl);
  
  bFirst = FALSE;
  /* Return the correction to the energy */
  return (-(Vself+Vexcl));
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
  t_psdata eps;
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
	    "ATOM",(i+1),*atoms->atomname[i],*atoms->resname[rnr],' ',
	    (rnr+1) % 10000,
	    10*(dx+x[i][XX]),10*x[i][YY],10*(x[i][ZZ]),0.0,phi[i]);
  }
  fclose(fp);
}

void write_grid_pqr(char *fn,int nx,int ny,int nz,real ***phi)
{
  FILE *fp;
  int  i,j,k,rnr=0;
  real fac=4.0;
  
  fp=ffopen(fn,"w");
  for(i=0; (i<nx); i++)
    for(j=0; (j<ny); j++)
      for(k=0; (k<nz); k++,rnr++)
	fprintf(fp,"%-6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
		"ATOM",(i+1),"C","C",' ',
		1+((rnr+1) % 10000),fac*i,fac*j,fac*k,0.0,phi[i][j][k]);
  fclose(fp);
}


real analyse_diff(FILE *log,char *label,
		  int natom,rvec ffour[],rvec fpppm[],
		  real phi_f[],real phi_p[],real phi_sr[],
		  char *fcorr,char *pcorr,char *ftotcorr,char *ptotcorr)
{
  int  i,m;
  FILE *fp=NULL,*gp=NULL;
  real f2sum=0,p2sum=0;
  real df,fmax,dp,pmax,rmsf;
  
  fmax = fabs(ffour[0][0]-fpppm[0][0]);
  pmax = fabs(phi_f[0] - phi_p[0]);
  
  for(i=0; (i<natom); i++) {
    dp     = fabs(phi_f[i] - phi_p[i]);
    pmax   = max(dp,pmax);
    p2sum += dp*dp;
    for(m=0; (m<DIM); m++) {
      df     = fabs(ffour[i][m] - fpppm[i][m]);
      fmax   = max(df,fmax);
      f2sum += df*df;
    }
  }
  
  rmsf = sqrt(f2sum/(3.0*natom));
  fprintf(log,"\n********************************\nERROR ANALYSIS for %s\n",
	  label);
  fprintf(log,"%-10s%12s%12s\n","Error:","Max Abs","RMS");
  fprintf(log,"%-10s  %10.3f  %10.3f\n","Force",fmax,rmsf);
  fprintf(log,"%-10s  %10.3f  %10.3f\n","Potential",pmax,sqrt(p2sum/(natom)));

  if (fcorr) {  
    fp = xvgropen(fcorr,"LR Force Correlation","Four-Force","PPPM-Force");
    for(i=0; (i<natom); i++) {
      for(m=0; (m<DIM); m++) {
	fprintf(fp,"%10.3f  %10.3f\n",ffour[i][m],fpppm[i][m]);
      }
    }
    ffclose(fp);
    do_view(fcorr,NULL);
  }
  if (pcorr)  
    fp = xvgropen(pcorr,"LR Potential Correlation","Four-Pot","PPPM-Pot");
  if (ptotcorr)
    gp = xvgropen(ptotcorr,"Total Potential Correlation","Four-Pot","PPPM-Pot");
  for(i=0; (i<natom); i++) {
    if (pcorr)
      fprintf(fp,"%10.3f  %10.3f\n",phi_f[i],phi_p[i]);
    if (ptotcorr)
      fprintf(gp,"%10.3f  %10.3f\n",phi_f[i]+phi_sr[i],phi_p[i]+phi_sr[i]);
  }
  if (pcorr) {
    ffclose(fp);
    do_view(pcorr,NULL);
  }
  if (ptotcorr) {
    ffclose(gp);
    do_view(ptotcorr,NULL);
  }

  return rmsf;
}

