/*
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
#include "coulomb.h"
#include "smalloc.h"
#include "physics.h"
#include "txtdump.h"
#include "futil.h"
#include "names.h"
#include "writeps.h"
#include "macros.h"
#include "xvgr.h"
#include "pppm.h"
#include "gmxfio.h"

#include "thread_mpi.h"

#define p2(x) ((x)*(x))
#define p3(x) ((x)*(x)*(x)) 
#define p4(x) ((x)*(x)*(x)*(x)) 

static real A,A_3,B,B_4,C,c1,c2,c3,c4,c5,c6,One_4pi,FourPi_V,Vol,N0;
#ifdef GMX_THREADS
static tMPI_Thread_mutex_t shift_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif


void set_shift_consts(FILE *log,real r1,real rc,rvec box,t_forcerec *fr)
{
#ifdef GMX_THREADS
  /* at the very least we shouldn't allow multiple threads to set these 
     simulataneously */
  tMPI_Thread_mutex_lock(&shift_mutex);
#endif
  /* A, B and C are recalculated in tables.c */
  if (r1 < rc) {
    A   = (2*r1-5*rc)/(p3(rc)*p2(rc-r1));
    B   = (4*rc-2*r1)/(p3(rc)*p3(rc-r1));
    /*C   = (10*rc*rc-5*rc*r1+r1*r1)/(6*rc*rc); Hermans Eq. not correct */
  }
  else
    gmx_fatal(FARGS,"r1 (%f) >= rc (%f) in %s, line %d",
	      r1,rc,__FILE__,__LINE__);

  A_3 = A/3.0;
  B_4 = B/4.0;
  C   = 1/rc-A_3*p3(rc-r1)-B_4*p4(rc-r1);
  N0  = 2.0*M_PI*p3(rc)*p3(rc-r1);

  Vol     =(box[XX]*box[YY]*box[ZZ]);
  FourPi_V=4.0*M_PI/Vol;

  if (debug) {
      fprintf(debug,"Constants for short-range and fourier stuff:\n"
	      "r1 = %10.3f,  rc = %10.3f\n"
	      "A  = %10.3e,  B  = %10.3e,  C  = %10.3e, FourPi_V = %10.3e\n",
	      r1,rc,A,B,C,FourPi_V);
  }

  /* Constants derived by Mathematica */
  c1 = -40*rc*rc    + 50*rc*r1    - 16*r1*r1;
  c2 =  60*rc       - 30*r1;
  c3 = -10*rc*rc*rc + 20*rc*rc*r1 - 13*rc*r1*r1 + 3*r1*r1*r1;
  c4 = -20*rc*rc    + 40*rc*r1    - 14*r1*r1;
  c5 = -c2;
  c6 = -5*rc*rc*r1  +  7*rc*r1*r1 - 2*r1*r1*r1;

  if (debug) {
      fprintf(debug,"c1 = %10.3e,  c2 = %10.3e,  c3 = %10.3e\n"
	      "c4 = %10.3e,  c5 = %10.3e,  c6 = %10.3e,  N0 = %10.3e\n",
	      c1,c2,c3,c4,c5,c6,N0);
  }
    
  One_4pi = 1.0/(4.0*M_PI);
#ifdef GMX_THREADS
  tMPI_Thread_mutex_unlock(&shift_mutex);
#endif
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



real shift_LRcorrection(FILE *fp,int start,int natoms,
			t_commrec *cr,t_forcerec *fr,
			real charge[],t_blocka *excl,rvec x[],
			gmx_bool bOld,matrix box,matrix lr_vir)
{
  static gmx_bool bFirst=TRUE;
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
	  dr_1    = gmx_invsqrt(dr2);
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
	  rvec_inc(fr->f_novirsum[k],df);
	  rvec_dec(fr->f_novirsum[i],df);
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
  
real phi_aver(int natoms,real phi[])
{
  real phitot;
  int  i;

  phitot=0;  
  for(i=0; (i<natoms); i++)
    phitot=phitot+phi[i];
    
  return (phitot/natoms);
}

real symmetrize_phi(FILE *log,int natoms,real phi[],gmx_bool bVerbose)
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



real analyse_diff(FILE *log,char *label,const output_env_t oenv,
		  int natom,rvec ffour[],rvec fpppm[],
		  real phi_f[],real phi_p[],real phi_sr[],
		  char *fcorr,char *pcorr,char *ftotcorr,char *ptotcorr)
/* Analyse difference between forces from fourier (_f) and other (_p)
 * LR solvers (and potential also).
 * If the filenames are given, xvgr files are written.
 * returns the root mean square error in the force.
 */
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
    fp = xvgropen(fcorr,"LR Force Correlation","Four-Force","PPPM-Force",oenv);
    for(i=0; (i<natom); i++) {
      for(m=0; (m<DIM); m++) {
	fprintf(fp,"%10.3f  %10.3f\n",ffour[i][m],fpppm[i][m]);
      }
    }
    gmx_fio_fclose(fp);
    do_view(oenv,fcorr,NULL);
  }
  if (pcorr)  
    fp = xvgropen(pcorr,"LR Potential Correlation","Four-Pot","PPPM-Pot",oenv);
  if (ptotcorr)
    gp = xvgropen(ptotcorr,"Total Potential Correlation","Four-Pot","PPPM-Pot",
                  oenv);
  for(i=0; (i<natom); i++) {
    if (pcorr)
      fprintf(fp,"%10.3f  %10.3f\n",phi_f[i],phi_p[i]);
    if (ptotcorr)
      fprintf(gp,"%10.3f  %10.3f\n",phi_f[i]+phi_sr[i],phi_p[i]+phi_sr[i]);
  }
  if (pcorr) {
    gmx_fio_fclose(fp);
    do_view(oenv,pcorr,NULL);
  }
  if (ptotcorr) {
    gmx_fio_fclose(gp);
    do_view(oenv,ptotcorr,NULL);
  }

  return rmsf;
}

