/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.2
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#include "math.h"
#include "string.h"
#include "copyrite.h"
#include "stdio.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "statutil.h"
#include "ewald.h"

enum { mGuillot, mAB1, mLjc, mMaaren, mGuillot_Maple, mHard_Wall, mNR };

static double erf2(double x)
{
  return -(4*x/(sqrt(M_PI)))*exp(-x*x);
}

static double erf1(double x)
{
  return (2/sqrt(M_PI))*exp(-x*x);
}

void do_hard(FILE *fp,double resolution,double efac,double delta)
{
  int    i,k,imax;
  double x,vr,vr2,vc,vc2;
  
  if (delta < 0)
    gmx_fatal(FARGS,"Delta should be >= 0 rather than %f\n",delta);
    
  imax     = 3.0/resolution;
  for(i=0; (i<=imax); i++) {
    x   =  i*resolution;
    
    if (x < delta) {
      /* Avoid very high numbers */
      vc = vc2 = 1/delta;
    }
    else {
      vc  = 1/(x);
      vc2 = 2/pow(x,3);
    }
    vr  = erfc(efac*(x-delta))/2;
    vr2 = (1-erf2(efac*(x-delta)))/2;
    fprintf(fp,"%12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n",
	    x,vr,vr2,0.0,0.0,vc,vc2);
  }

}

void do_AB1(FILE *fp,int eel,double resolution,int ndisp,int nrep)
{
  int    i,k,imax;
  double myfac[3] = { 1, -1, 1 };
  double myexp[3] = { 1, 6, 0 };
  double x,v,v2;
  
  myexp[1] = ndisp;
  myexp[2] = nrep;
  imax     = 3.0/resolution;
  for(i=0; (i<=imax); i++) {
    x   =  i*resolution;
    
    fprintf(fp,"%12.5e",x);
    
    for(k=0; (k<3); k++) {
      if (x < 0.04) {
	/* Avoid very high numbers */
	v = v2 = 0;
      }
      else {
	v  =  myfac[k]*pow(x,-myexp[k]);
	v2 = (myexp[k]+1)*(myexp[k])*v/(x*x); 
      }
      fprintf(fp,"  %12.5e  %12.5e",v,v2);
    }
    fprintf(fp,"\n");
  }
}

void lo_do_ljc(double r,
	       double *vc,double *vc2,
	       double *vd,double *vd2,
	       double *vr,double *vr2)
{
  double r2,r6,r12;
  
  r2    = r*r;
  r6    = 1.0/(r2*r2*r2);
  r12   = r6*r6;

  *vc   = 1.0/r;
  *vc2  = 2.0/(r*r2);

  *vd   = -r6;
  *vd2  = 42.0*(*vd)/r2;

  *vr  = r12;
  *vr2 = 156.0*(*vr)/r2;
}

/* use with coulombtype = user */
void lo_do_ljc_pme(double r,
		   double rcoulomb, double ewald_rtol,
		   double *vc,double *vc2,
		   double *vd,double *vd2,
		   double *vr,double *vr2)
{
  double r2,r6,r12;
  double isp= 0.564189583547756;
  double ewc;

  ewc = calc_ewaldcoeff(rcoulomb,ewald_rtol);

  r2    = r*r;
  r6    = 1.0/(r2*r2*r2);
  r12   = r6*r6;
  
  *vc   = erfc(ewc*r)/r;
  *vc2  = 2*erfc(ewc*r)/(r*r2)+4*exp(-(ewc*ewc*r2))*ewc*isp/r2+
          4*ewc*ewc*ewc*exp(-(ewc*ewc*r2))*isp;

  *vd   = -r6;
  *vd2  = 42.0*(*vd)/r2;

  *vr  = r12;
  *vr2 = 156.0*(*vr)/r2;
}

void lo_do_guillot(double r,double xi, double xir,
			    double *vc,double *vc2,
			    double *vd,double *vd2,
			    double *vr,double *vr2)
{
  double qO     = -0.888;
  double qOd    = 0.226;
  double f0     = qOd/qO;
  double sqpi   = sqrt(M_PI);
  double r1,r2,z;
  
  r1    = r/(2*xi);
  r2    = r/(sqrt(2)*xi);
  *vc   = (1+sqr(f0)*erf(r1) + 2*f0*erf(r2))/r;
  *vc2  = ((2/sqr(r))*(*vc -
		       sqr(f0)*erf1(r1)/(2*xi) -
		       4*f0*erf1(r2)/sqrt(2)*xi) + 
	   (1/r)*(sqr(f0/(2.0*xi))*erf2(r1) + (2*f0/sqr(xi)))*erf2(r2));
  *vd   = -1.0/(r*r*r*r*r*r);
  *vd2  = 42.0*(*vd)/(r*r);
  z     = r/(2.0*xir);
  *vr   = erfc(z)/z;
  *vr2  = (sqpi*(*vr)/(2.0*z*z)+(1.0/(z*z)+1)*exp(-z*z))/(sqpi*sqr(xir));
}

void lo_do_guillot_maple(double r,double xi,double xir,
			 double *vc,double *vc2,
			 double *vd,double *vd2,
			 double *vr,double *vr2)
{
  double qO     = -0.888;
  double qOd    = 0.226;
  double f0     = qOd/qO;
  double sqpi   = sqrt(M_PI);

  *vc = pow(-f0/(1.0+f0)+1.0,2.0)/r+pow(-f0/(1.0+f0)+1.0,2.0)*f0*f0*erf(r/xi/2.0)/r+2.0*pow(-f0/(1.0+f0)+1.0,2.0)*f0*erf(r*sqrt(2.0)/xi/2.0)/r;
  *vc2 = 2.0*pow(-f0/(1.0+f0)+1.0,2.0)/(r*r*r)-pow(-f0/(1.0+f0)+1.0,2.0)*f0*f0/sqrt(M_PI)/(xi*xi*xi)*exp(-r*r/(xi*xi)/4.0)/2.0-2.0*pow(-f0/(1.0+f0)+1.0,2.0)*f0*f0/sqrt(M_PI)*exp(-r*r/(xi*xi)/4.0)/xi/(r*r)+2.0*pow(-f0/(1.0+f0)+1.0,2.0)*f0*f0*erf(r/xi/2.0)/(r*r*r)-2.0*pow(-f0/(1.0+f0)+1.0,2.0)*f0/sqrt(M_PI)/(xi*xi*xi)*exp(-r*r/(xi*xi)/2.0)*sqrt(2.0)-4.0*pow(-f0/(1.0+f0)+1.0,2.0)*f0/sqrt(M_PI)*exp(-r*r/(xi*xi)/2.0)*sqrt(2.0)/xi/(r*r)+4.0*pow(-f0/(1.0+f0)+1.0,2.0)*f0*erf(r*sqrt(2.0)/xi/2.0)/(r*r*r);
  
  *vd   = -1.0/(r*r*r*r*r*r);
  *vd2  = -42.0/(r*r*r*r*r*r*r*r);
  *vr   = 2.0*erfc(r/xir/2.0)/r*xir;
  *vr2  = 1.0/sqrt(M_PI)/(xir*xir)*exp(-r*r/(xir*xir)/4.0)+4.0/sqrt(M_PI)*exp(-r*r/(xir*xir)/4.0)/(r*r)+4.0*erfc(r/xir/2.0)/(r*r*r)*xir;
}

static void do_guillot(FILE *fp,int eel,double resolution,double rc,double rtol,double xi,double xir)
{
  int    i,i0,imax;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3/resolution;
  for(i=0; (i<=imax); i++) {
    r     = i*resolution;
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else 
      lo_do_guillot(r,xi,xir,&vc,&vc2,&vd,&vd2,&vr,&vr2);
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
}

static void do_ljc(FILE *fp,int eel,double resolution,real rc,real rtol)
{
  int    i,i0,imax;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3/resolution;
  for(i=0; (i<=imax); i++) {
    r     = i*resolution;
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    } else {
      if (eel == eelPME) {
	lo_do_ljc_pme(r,rc,rtol,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      } else if (eel == eelCUT) { 
	lo_do_ljc(r,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      }
    }
    fprintf(fp,"%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
}

static void do_guillot_maple(FILE *fp,int eel,double resolution,double rc,double rtol,double xi,double xir)
{
  int    i,i0,imax;
  //  double xi     = 0.15;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3/resolution;
  for(i=0; (i<=imax); i++) {
    r     = i*resolution;
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else
      if (eel == eelPME) {
	fprintf(fp, "Not implemented\n");
      } else if (eel == eelCUT) { 
	lo_do_guillot_maple(r,xi,xir,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      }
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
} 

static void do_maaren(FILE *fp,int eel,double resolution,int npow)
{
  int    i,i0,imax;
  double xi     = 0.05;
  double xir     = 0.0615;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3/resolution;
  for(i=0; (i<=imax); i++) {
    r     = i*resolution;
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else {
      lo_do_guillot_maple(r,xi,xir,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      vr  =  pow(r,-1.0*npow);
      vr2 = (npow+1.0)*(npow)*vr/sqr(r); 
    }
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gen_table generates tables for mdrun for use with the USER defined",
    "potentials."
  };
  static char *opt[]     = { NULL, "cut", "rf", "pme", NULL };
  static char *model[]   = { NULL, "guillot", "AB1", "ljc", "maaren", "guillot_maple", "hard_wall", NULL };
  static real resolution = 0.001,delta=0,efac=500,rc=0.9,rtol=1e-05,xi=0.15,xir=0.0615;
  static int  nrep       = 12;
  static int  ndisp      = 6;
  t_pargs pa[] = {
    { "-el",     FALSE, etENUM, {opt},
      "Electrostatics type: cut, rf or pme" },
    { "-rc",     FALSE, etREAL, {&rc},
      "Cut-off required for rf or pme" },
    { "-rtol",   FALSE, etREAL, {&rtol},
      "Ewald tolerance required for pme" },
    { "-xi",   FALSE, etREAL, {&xi},
      "Width of the Gaussian diffuse charge of the G&G model" },
    { "-xir",   FALSE, etREAL, {&xir},
      "Width of erfc(z)/z repulsion of the G&G model (z=0.5 rOO/xir)" },
    { "-m",      FALSE, etENUM, {model},
      "Model for the tables" },
    { "-resol",  FALSE, etREAL, {&resolution},
      "Resolution of the table (nm)" },
    { "-delta",  FALSE, etREAL, {&delta},
      "Displacement in the Coulomb functions (nm), used as 1/(r+delta). Only for hard wall potential." },
    { "-efac",   FALSE, etREAL, {&efac},
      "Number indicating the steepness of the hardwall potential." },
    { "-nrep",   FALSE, etINT,  {&nrep},
      "Power for the repulsion potential (with model AB1 or maaren)" },
    { "-ndisp",   FALSE, etINT,  {&ndisp},
      "Power for the dispersion potential (with model AB1 or maaren)" }
  };
#define NPA asize(pa)
  t_filenm fnm[] = {
    { efXVG, "-o", "table", ffWRITE }
  };
#define NFILE asize(fnm)
  FILE *fp;
  int  eel=0,m=0;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);
  
  if (strcmp(opt[0],"cut") == 0) 
    eel = eelCUT;
  else if (strcmp(opt[0],"rf") == 0) 
    eel = eelRF;
  else if (strcmp(opt[0],"pme") == 0) 
    eel = eelPME;
  else 
    gmx_fatal(FARGS,"Invalid argument %s for option -e",opt[0]);
  if (strcmp(model[0],"maaren") == 0) 
    m = mMaaren;
  else if (strcmp(model[0],"AB1") == 0) 
    m = mAB1;
  else if (strcmp(model[0],"ljc") == 0) 
    m = mLjc;
  else if (strcmp(model[0],"guillot") == 0) 
    m = mGuillot;
  else if (strcmp(model[0],"guillot_maple") == 0) 
    m = mGuillot_Maple;
  else if (strcmp(model[0],"hard_wall") == 0) 
    m = mHard_Wall;
  else 
    gmx_fatal(FARGS,"Invalid argument %s for option -m",opt[0]);
    
  fp = ffopen(opt2fn("-o",NFILE,fnm),"w");
  switch (m) {
  case mGuillot:
    do_guillot(fp,eel,resolution,rc,rtol,xi,xir);
    break;
  case mGuillot_Maple:
    fprintf(fp, "#\n# Table Guillot_Maple: rc=%g, rtol=%g, xi=%g, xir=%g\n#\n",rc,rtol,xi,xir);
    do_guillot_maple(fp,eel,resolution,rc,rtol,xi,xir);
    break;
  case mMaaren:
    do_maaren(fp,eel,resolution,nrep);
    break;
  case mAB1:
    fprintf(fp, "#\n# Table AB1: ndisp=%d nrep=%d\n#\n",ndisp,nrep);
    do_AB1(fp,eel,resolution,ndisp,nrep);
    break;
  case mLjc:
    fprintf(fp, "#\n# Table LJC(12-6-1): rc=%g, rtol=%g\n#\n",rc,rtol);
    do_ljc(fp,eel,resolution,rc,rtol);
    break;
  case mHard_Wall:
    do_hard(fp,resolution,efac,delta);
    break;
  default:
    gmx_fatal(FARGS,"Model %s not supported yet",model[0]);
  }  
  fclose(fp);
  
  thanx(stdout);
  
  return 0;
}
