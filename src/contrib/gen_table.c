/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "copyrite.h"
#include "typedefs.h"
#include "macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/commandline/pargs.h"
#include "coulomb.h"

enum { mGuillot2001a, mAB1, mLjc, mMaaren, mGuillot_Maple, mHard_Wall, mGG, mGG_qd_q, mGG_qd_qd, mGG_q_q, mNR };

static double erf2(double x)
{
  return -(2*x*M_2_SQRTPI)*exp(-x*x);
}

static double erf1(double x)
{
  return M_2_SQRTPI*exp(-x*x);
}

static void do_hard(FILE *fp,int pts_nm,double efac,double delta)
{
  int    i,k,imax;
  double x,vr,vr2,vc,vc2;
  
  if (delta < 0)
    gmx_fatal(FARGS,"Delta should be >= 0 rather than %f\n",delta);
    
  imax     = 3.0*pts_nm;
  for(i=0; (i<=imax); i++) {
    x   =  i*(1.0/pts_nm);
    
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

static void do_AB1(FILE *fp,int eel,int pts_nm,int ndisp,int nrep)
{
  int    i,k,imax;
  double myfac[3] = { 1, -1, 1 };
  double myexp[3] = { 1, 6, 0 };
  double x,v,v2;
  
  myexp[1] = ndisp;
  myexp[2] = nrep;
  imax     = 3.0*pts_nm;
  for(i=0; (i<=imax); i++) {
    x   =  i*(1.0/pts_nm);
    
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

static void lo_do_ljc(double r,
		      double *vc,double *fc,
		      double *vd,double *fd,
		      double *vr,double *fr)
{
  double r2,r_6,r_12;
  
  r2    = r*r;
  r_6   = 1.0/(r2*r2*r2);
  r_12  = r_6*r_6;

  *vc   = 1.0/r;            /*  f(x)     Coulomb    */
  *fc   = 1.0/(r2);         /* -f'(x)               */
  
  *vd   = -r_6;             /*  g(c)     Dispersion */
  *fd   =  6.0*(*vd)/r;     /* -g'(x)               */

  *vr   = r_12;             /*  h(x)     Repulsion  */
  *fr   = 12.0*(*vr)/r;     /* -h'(x)               */
}

/* use with coulombtype = user */
static void lo_do_ljc_pme(double r,
			  double rcoulomb, double ewald_rtol,
			  double *vc,double *fc,
			  double *vd,double *fd,
			  double *vr,double *fr)
{
  double r2,r_6,r_12;
  double ewc;

  ewc  = calc_ewaldcoeff(rcoulomb,ewald_rtol);

  r2   = r*r;
  r_6  = 1.0/(r2*r2*r2);
  r_12 = r_6*r_6;
  
  *vc   = erfc(ewc*r)/r;
  /* *vc2  = 2*erfc(ewc*r)/(r*r2)+2*exp(-(ewc*ewc*r2))*ewc*M_2_SQRTPI/r2+
     2*ewc*ewc*ewc*exp(-(ewc*ewc*r2))*M_2_SQRTPI;*/
  *fc  = ewc*exp(-ewc*ewc*r2)*M_2_SQRTPI/r + erfc(ewc*r)/r2;

  *vd  = -r_6;
  *fd  = -6.0*(*vd)/r;

  *vr  = r_12;
  *fr  = 12.0*(*vr)/r;
}

static void lo_do_guillot(double r,double xi, double xir,
			  double *vc,double *fc,
			  double *vd,double *fd,
			  double *vr,double *fr)
{
  double qO     = -0.888;
  double qOd    =  0.226;
  double f0     = qOd/qO;
  double sqpi   = sqrt(M_PI);
  double rxi1,rxi2,z;
  double r2,r_6;

  r2   = r*r;
  r_6  = 1.0/(r2*r2*r2);
  
  rxi1    = r/(2*xi);
  rxi2    = r/(sqrt(2)*xi);
  *vc   = (1 + f0*f0*erf(r/(2*xi)) + 2*f0*erf(r/(sqrt(2)*xi)) )/r;

  *fc   =  f0*f0*erf(r/(2*xi)) + 2*f0*erf(r/(sqrt(2)*xi));
    ;
 /* MuPad: Uc := erf(r/(2*xi))/r +  

     Mathematica:
     r1 := r/(2*xi);
     r2 := r/(Sqrt[2] * xi);
     Uc[r_] := (1 + f0 * f0 * Erf[r/(2*xi)] + 2 * f0 * Erf[r/(Sqrt[2]*xi)]) / r;
     -D[Uc[r],r]
     CForm= 
     -(((2*f0*Sqrt(2/Pi))/(Power(E,Power(r,2)/(2.*Power(xi,2)))*xi) + 
     Power(f0,2)/(Power(E,Power(r,2)/(4.*Power(xi,2)))*Sqrt(Pi)*xi))/r) + 
     (1 + Power(f0,2)*Erf(r/(2.*xi)) + 2*f0*Erf(r/(Sqrt(2)*xi)))/Power(r,2)

     
Uc1[r_] := 1/r;
-D[Uc1[r],r]
          -2
Out[20]= r

Uc2[r_] := f0^2 * Erf[r1] / r;
-D[Uc2[r],r]


Uc3[r_] := 2 * f0 * Erf[r2]/ r;
-D[Uc3[r],r]

Uc[r_] := Erf[r/(2*xi)] / r

D[Uc[r],r]


D[Erf[r],r]

*/
    *vc   = (1 + sqr(f0)*erf(rxi1) + 2*f0*erf(rxi2))/r;
    *fc   = 
      (1/r 
	+ (- f0 * (2 * sqrt(2) + exp(r2/4*xi*xi)*f0)/(exp(r2/(2*xi*xi))*sqrt(M_PI)*xi) + f0*f0*erf(r/(2*xi)) + 2 *f0 * erf(r/(sqrt(2 * xi)))  )/r2)
      ;


  /*  *vc2  = ((2/sqr(r))*(*vc -
		       sqr(f0)*erf1(r1)/(2*xi) -
		       4*f0*erf1(r2)/sqrt(2)*xi) + 
		       (1/r)*(sqr(f0/(2.0*xi))*erf2(r1) + (2*f0/sqr(xi)))*erf2(r2)); */

  *vd  = -r_6;
  *fd  = -6.0*(*vd)/r;

  z     = r/(2.0*xir);
  *vr   = erfc(z)/z;
  *fr   = 0.0;
  //  *vr2  = (sqpi*(*vr)/(2.0*z*z)+(1.0/(z*z)+1)*exp(-z*z))/(sqpi*sqr(xir));
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

static void lo_do_GG(double r,double xi,double xir,
		     double *vc,double *fc,
		     double *vd,double *fd,
		     double *vr,double *fr)
{
  double qO     = -0.888;
  double qOd    =  0.226;
  double f0     = qOd/qO;
  double sqpi   = sqrt(M_PI);
  double r2,xi2;

  r2 = r*r;
  xi2 = xi*xi;

  *vc = 1.0/r + f0*f0*erf(r/(2*xi))/r + 2*f0*erf(r/(sqrt(2)*xi))/r;

  // -D[1/r,r] -D[f0*f0*Erf[r/(2*xi)]/r,r] -D[2*f0*Erf[r/(Sqrt[2]*xi)]/r,r]
  *fc  = (
    1.0/r2 +
    f0*f0*(-exp(-r2/(4*xi2))/(sqrt(M_PI) * r * xi) + erf(r/(2*xi))/r2) +
    2*f0*(-sqrt(2.0/M_PI)*exp(-r2/(2*xi2))/ (r*xi) + erf(r/(sqrt(2)*xi))/r2)
    );

  // -D[1/r^6,r]
  *vd  = -1.0/(r*r*r*r*r*r);
  *fd  = 6.0*(*vd)/r;
  
  //  -D[2*xir*Erfc[r/(2*xir)]/r,r]
  *vr  = 2.*xir*erfc(r/(2.*xir))/r;
  *fr  = -(-2.*exp(-r2/(4*xir*xir)) / (sqrt(M_PI)*r)  - 2*xir*erfc(r/(2*xir))/r2  );

}

/* Guillot2001 diffuse charge - diffuse charge interaction
   Mathematica

In[19]:= Uc[r_] := Erf[r/(2*xi)]/r

In[20]:= -D[Uc[r],r]

                                             r
                                        Erf[----]
                       1                    2 xi
Out[20]= -(-------------------------) + ---------
             2      2                       2
            r /(4 xi )                     r
           E           Sqrt[Pi] r xi
*/
void lo_do_GG_qd_qd(double r,double xi,double xir,
		    double *vc,double *fc,
		    double *vd,double *fd,
		    double *vr,double *fr)
{
  double sqpi   = sqrt(M_PI);

  *vc = erf(r/(2*xi))/r; 
    //erf((r*(1.0/2.0))/xi)/r;
  *fc = -(1.0/(exp(r*r/(4*xi*xi))*sqpi*r*xi)) + (erf(r/2*xi)/(r*r));

    //2.0*pow(r, -3.0)*erf((r*(1.0/2.0))/xi) - (1.0/2.0)*pow(M_PI, -1.0/2.0)*pow(xi, -3.0)*exp((-1.0/4.0)*(r*r)*pow(xi, -2.0)) - (2.0*pow(M_PI, -1.0/2.0)*pow(r, -2.0)*exp((-1.0/4.0)*(r*r)*pow(xi, -2.0)))/xi ;
  *vd  = 0.0;
  *fd  = 0.0;
  *vr  = 0.0;
  *fr  = 0.0;
}

/* Guillot2001 charge - diffuse charge interaction eqn 4 & 5
   Mathematica
In[17]:= Uc[r_] := Erf[r/(Sqrt[2]*xi)]/r

In[18]:= -D[Uc[r],r]

                    2                  r
               Sqrt[--]        Erf[----------]
                    Pi             Sqrt[2] xi
Out[18]= -(----------------) + ---------------
             2      2                 2
            r /(2 xi )               r
           E           r xi
*/
void lo_do_GG_q_qd(double r,double xi,double xir,
		   double *vc,double *fc,
		   double *vd,double *fd,
		   double *vr,double *fr)
{
  double sqpi   = sqrt(M_PI);

  *vc = erf(r/(sqrt(2)*xi)) / r;
    //erf(((1.0/2.0)*pow(2.0, 1.0/2.0)*r)/xi)/r ;
  *fc = -(sqrt(2/M_PI)/(exp(r*r/(2*xi*xi))*r*xi)) + (erf(r/(sqrt(2)*xi))/(r*r));
    //2.0*pow(r, -3.0)*erf(((1.0/2.0)*pow(2.0, 1.0/2.0)*r)/xi) - pow(2.0, 1.0/2.0)*pow(M_PI, -1.0/2.0)*pow(xi, -3.0)*exp((-1.0/2.0)*(r*r)*pow(xi, -2.0)) - (2.0*pow(2.0, 1.0/2.0)*pow(M_PI, -1.0/2.0)*pow(r, -2.0)*exp((-1.0/2.0)*(r*r)*pow(xi, -2.0)))/xi ;

  *vd  = 0.0;
  *fd  = 0.0;
  *vr  = 0.0;
  *fr  = 0.0;
}

/* Guillot2001 charge - charge interaction (normal coulomb), repulsion and dispersion
   Mathematica

In[6]:= Uc[r_] := 1.0/r

In[7]:= -D[Uc[r],r]

        1.
Out[7]= --
         2
        r

In[8]:= Ud[r_] := -1.0/r^6

In[9]:= -D[Ud[r],r]

        -6.
Out[9]= ---
         7
        r

In[13]:= Ur[r_] := (2*xir)*Erfc[r/(2*xir)]/r

In[14]:= -D[Ur[r],r]
                                                r
                                   2 xir Erfc[-----]
                    2                         2 xir
Out[16]= ----------------------- + -----------------
           2       2                       2
          r /(4 xir )                     r
         E            Sqrt[Pi] r


*/
void lo_do_GG_q_q(double r,double xi,double xir,
		  double *vc,double *fc,
		  double *vd,double *fd,
		  double *vr,double *fr)
{
  double sqpi   = sqrt(M_PI);

  *vc  = 1.0/r;
  *fc  = 1.0/(r*r);

  *vd  = -1.0/(r*r*r*r*r*r);
  *fd  = -6.0/(r*r*r*r*r*r*r);

  *vr  = (2.0*xir*erfc(r/(2.0*xir)))/r;
  *fr  = 2.0/(exp((r*r)/(4*xir*xir)) * sqpi *r) + (2*xir*erfc((r*xir)/2.0))/(r*r);
    //4.0*pow(M_PI, -1.0/2.0)*pow(r, -2.0)*exp((-1.0/4.0)*(r*r)*pow(xir, -2.0)) + pow(M_PI, -1.0/2.0)*pow(xir, -2.0)*exp((-1.0/4.0)*(r*r)*pow(xir, -2.0)) + 4.0*pow(r, -3.0)*xir*erfc((r*(1.0/2.0))/xir);
}

static void do_guillot(FILE *fp,int eel,int pts_nm,double rc,double rtol,double xi,double xir)
{
  int    i,i0,imax;
  double r,vc,fc,vd,fd,vr,fr;

  imax = 3*pts_nm;
  for(i=0; (i<=imax); i++) {
    r     = i*(1.0/pts_nm);
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = fc = vd = fd = vr = fr = 0;
    }
    else 
      if (eel == eelPME) {
	fprintf(fp, "Not implemented\n");
      } else if (eel == eelCUT) { 
	lo_do_guillot(r,xi,xir,&vc,&fc,&vd,&fd,&vr,&fr);
      }
    fprintf(fp,"%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n",
	    r,vc,fc,vd,fd,vr,fr);

  }
}

/* TODO: 
   PvM: Everything is hardcoded, we should fix that. How?
*/
static void do_guillot2001a(const char *file,int eel,int pts_nm,double rc,double rtol,double xi,double xir)
{
  FILE *fp=NULL;
  static char buf[256];
  static char *atype[]   = { "HW", "OW", "HWd", "OWd", NULL };
  int    i,j,k,i0,imax,atypemax=4;
  double r,vc,fc,vd,fd,vr,fr;

  /* For Guillot2001a we have four types: HW, OW, HWd and OWd. */

  for (j=0;(j<atypemax);j++) {           /* loops over types */
    for (k=0; (k<=j); k++) {                    
      sprintf(buf,"table_%s_%s.xvg",atype[k],atype[j]);
      
      printf("%d %d %s\n", j, k, buf);
      /* Guillot2001a eqn 2, 6 and 7 */
      if (((strcmp(atype[j],"HW") == 0) && (strcmp(atype[k],"HW") == 0)) ||
	  ((strcmp(atype[j],"OW") == 0) && (strcmp(atype[k],"HW") == 0)) ||
	  ((strcmp(atype[j],"OW") == 0) && (strcmp(atype[k],"OW") == 0))) {

	fp = gmx_ffopen(buf,"w");
  
	imax = 3*pts_nm;
	for(i=0; (i<=imax); i++) {
	  r     = i*(1.0/pts_nm);
	  /* Avoid very large numbers */
	  if (r < 0.04) {
	    vc = fc = vd = fd = vr = fr = 0;
	  }
	  else 
	    if (eel == eelPME || eel == eelRF) {
	      fprintf(stderr, "Not implemented\n");
	      exit(1);
	    } else if (eel == eelCUT) { 
	      lo_do_GG_q_q(r,xi,xir,&vc,&fc,&vd,&fd,&vr,&fr);
	    }
	  fprintf(fp,"%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n",
		  r,vc,fc,vd,fd,vr,fr);
	  
	}
	gmx_ffclose(fp);
     
	/* Guillot eqn 4 and 5 */
      } else if (((strcmp(atype[j],"HWd") == 0) && (strcmp(atype[k],"HW") == 0)) ||
		 ((strcmp(atype[j],"HWd") == 0) && (strcmp(atype[k],"OW") == 0)) ||
		 ((strcmp(atype[j],"OWd") == 0) && (strcmp(atype[k],"HW") == 0)) ||
		 ((strcmp(atype[j],"OWd") == 0) && (strcmp(atype[k],"OW") == 0))) {
	
	fp = gmx_ffopen(buf,"w");
  
	imax = 3*pts_nm;
	for(i=0; (i<=imax); i++) {
	  r     = i*(1.0/pts_nm);
	  /* Avoid very large numbers */
	  if (r < 0.04) {
	    vc = fc = vd = fd = vr = fr = 0;
	  }
	  else 
	    if (eel == eelPME || eel == eelRF) {
	      fprintf(stderr, "Not implemented\n");
	      exit(1);
	    } else if (eel == eelCUT) { 
	      lo_do_GG_q_qd(r,xi,xir,&vc,&fc,&vd,&fd,&vr,&fr);
	    }
	  fprintf(fp,"%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n",
		  r,vc,fc,vd,fd,vr,fr);
	  
	}
	gmx_ffclose(fp);

	/* Guillot2001a eqn 3 */
      } else if (((strcmp(atype[j],"HWd") == 0) && (strcmp(atype[k],"HWd") == 0)) ||
		 ((strcmp(atype[j],"OWd") == 0) && (strcmp(atype[k],"HWd") == 0)) ||
		 ((strcmp(atype[j],"OWd") == 0) && (strcmp(atype[k],"OWd") == 0))) {

	fp = gmx_ffopen(buf,"w");
  
	imax = 3*pts_nm;
	for(i=0; (i<=imax); i++) {
	  r     = i*(1.0/pts_nm);
	  /* Avoid very large numbers */
	  if (r < 0.04) {
	    vc = fc = vd = fd = vr = fr = 0;
	  }
	  else 
	    if (eel == eelPME || eel == eelRF) {
	      fprintf(stderr, "Not implemented\n");
	      exit(1);
	    } else if (eel == eelCUT) { 
	      lo_do_GG_qd_qd(r,xi,xir,&vc,&fc,&vd,&fd,&vr,&fr);
	    }
	  fprintf(fp,"%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n",
		  r,vc,fc,vd,fd,vr,fr);
	  
	}
	gmx_ffclose(fp);

      } else 
	gmx_fatal(FARGS,"Invalid atom type: %s %s", atype[j], atype[k]);
      
      
    }
  }
}

static void do_ljc(FILE *fp,int eel,int pts_nm,real rc,real rtol)
{ 
  int    i,i0,imax;
  double r,vc,fc,vd,fd,vr,fr;

  imax = 3*pts_nm;
  for(i=0; (i<=imax); i++) {
    r     = i*(1.0/pts_nm);
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = fc = vd = fd = vr = fr = 0;
    } else {
      if (eel == eelPME) {
	lo_do_ljc_pme(r,rc,rtol,&vc,&fc,&vd,&fd,&vr,&fr);
      } else if (eel == eelCUT) { 
	lo_do_ljc(r,&vc,&fc,&vd,&fd,&vr,&fr);
      }
    }
    fprintf(fp,"%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n",
	    r,vc,fc,vd,fd,vr,fr);
  }
}

static void do_guillot_maple(FILE *fp,int eel,int pts_nm,double rc,double rtol,double xi,double xir)
{
  int    i,i0,imax;
  /*  double xi     = 0.15;*/
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3*pts_nm;
  for(i=0; (i<=imax); i++) {
    r     = i*(1.0/pts_nm);
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

static void do_GG(FILE *fp,int eel,int pts_nm,double rc,double rtol,double xi,double xir)
{
  int    i,i0,imax;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3*pts_nm;
  for(i=0; (i<=imax); i++) {
    r     = i*(1.0/pts_nm);
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else
      if (eel == eelPME) {
	fprintf(fp, "Not implemented\n");
      } else if (eel == eelCUT) { 
	lo_do_GG(r,xi,xir,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      }
    fprintf(fp,"%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
} 

static void do_GG_q_q(FILE *fp,int eel,int pts_nm,double rc,double rtol,double xi,double xir)
{
  int    i,i0,imax;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3*pts_nm;
  for(i=0; (i<=imax); i++) {
    r     = i*(1.0/pts_nm);
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else
      if (eel == eelPME) {
	fprintf(fp, "Not implemented\n");
      } else if (eel == eelCUT) { 
	lo_do_GG_q_q(r,xi,xir,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      }
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
} 

static void do_GG_q_qd(FILE *fp,int eel,int pts_nm,double rc,double rtol,double xi,double xir)
{
  int    i,i0,imax;
  /*  double xi     = 0.15;*/
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3*pts_nm;
  for(i=0; (i<=imax); i++) {
    r     = i*(1.0/pts_nm);
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else
      if (eel == eelPME) {
	fprintf(fp, "Not implemented\n");
      } else if (eel == eelCUT) { 
	lo_do_GG_q_qd(r,xi,xir,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      }
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
} 

static void do_GG_qd_qd(FILE *fp,int eel,int pts_nm,double rc,double rtol,double xi,double xir)
{
  int    i,i0,imax;
  /*  double xi     = 0.15;*/
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3*pts_nm;
  for(i=0; (i<=imax); i++) {
    r     = i*(1.0/pts_nm);
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else
      if (eel == eelPME) {
	fprintf(fp, "Not implemented\n");
      } else if (eel == eelCUT) { 
	lo_do_GG_qd_qd(r,xi,xir,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      }
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
} 

static void do_maaren(FILE *fp,int eel,int pts_nm,int npow)
{
  int    i,i0,imax;
  double xi     = 0.05;
  double xir     = 0.0615;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3*pts_nm;
  for(i=0; (i<=imax); i++) {
    r     = i*(1.0/pts_nm);
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
    "[TT]gen_table[tt] generates tables for [TT]mdrun[tt] for use with the USER defined",
    "potentials. Note that the format has been update for higher",
    "accuracy in the forces starting with version 4.0. Using older",
    "tables with 4.0 will silently crash your simulations, as will",
    "using new tables with an older GROMACS version. This is because in the",
    "old version the second derevative of the potential was specified",
    "whereas in the new version the first derivative of the potential",
    "is used instead.[PAR]"
  };
  static char *opt[]     = { NULL, "cut", "rf", "pme", NULL };
  /*  static char *model[]   = { NULL, "guillot", "AB1", "ljc", "maaren", "guillot_maple", "hard_wall", "gg_q_q", "gg_qd_q", "gg_qd_qd", NULL }; */
  static char *model[]   = { NULL, "ljc", "gg", "guillot2001a",  
			     NULL };
  static real delta=0,efac=500,rc=0.9,rtol=1e-05,xi=0.15,xir=0.0615;
  static real w1=20,w2=20;
  static int  nrow1=1,nrow2=1;
  static int  nrep       = 12;
  static int  ndisp      = 6;
  static int  pts_nm     = 500;
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
    { "-resol",  FALSE, etINT,  {&pts_nm},
      "Resolution of the table (points per nm)" },
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
  char *fn;
  int  eel=0,m=0;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,
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
  else if (strcmp(model[0],"guillot2001a") == 0) 
    m = mGuillot2001a;
  else if (strcmp(model[0],"guillot_maple") == 0) 
    m = mGuillot_Maple;
  else if (strcmp(model[0],"hard_wall") == 0) 
    m = mHard_Wall;
  else if (strcmp(model[0],"gg") == 0) 
    m = mGG;
  else if (strcmp(model[0],"gg_qd_q") == 0) 
    m = mGG_qd_q;
  else if (strcmp(model[0],"gg_qd_qd") == 0) 
    m = mGG_qd_qd;
  else if (strcmp(model[0],"gg_q_q") == 0) 
    m = mGG_q_q;
  else 
    gmx_fatal(FARGS,"Invalid argument %s for option -m",opt[0]);
    
  fn = opt2fn("-o",NFILE,fnm);
  if ((m != mGuillot2001a)) 
    fp = gmx_ffopen(fn,"w");
  switch (m) {
  case mGuillot2001a:
    do_guillot2001a(fn,eel,pts_nm,rc,rtol,xi,xir);
    break;
  case mGuillot_Maple:
    fprintf(fp, "#\n# Table Guillot_Maple: rc=%g, rtol=%g, xi=%g, xir=%g\n#\n",rc,rtol,xi,xir);
    do_guillot_maple(fp,eel,pts_nm,rc,rtol,xi,xir);
    break;
  case mGG_q_q:
    fprintf(fp, "#\n# Table GG_q_q: rc=%g, rtol=%g, xi=%g, xir=%g\n#\n",rc,rtol,xi,xir);
    do_GG_q_q(fp,eel,pts_nm,rc,rtol,xi,xir);
    break;
  case mGG:
    fprintf(fp, "#\n# Table GG: rc=%g, rtol=%g, xi=%g, xir=%g\n#\n",rc,rtol,xi,xir);
    do_GG(fp,eel,pts_nm,rc,rtol,xi,xir);
    break;
  case mGG_qd_q:
    fprintf(stdout, "case mGG_qd_q");
    fprintf(fp, "#\n# Table GG_qd_q: rc=%g, rtol=%g, xi=%g, xir=%g\n#\n",rc,rtol,xi,xir);
    do_GG_q_qd(fp,eel,pts_nm,rc,rtol,xi,xir);
    break;
  case mGG_qd_qd:
    fprintf(stdout, "case mGG_qd_qd");
    fprintf(fp, "#\n# Table GG_qd_qd: rc=%g, rtol=%g, xi=%g, xir=%g\n#\n",rc,rtol,xi,xir);
    do_GG_qd_qd(fp,eel,pts_nm,rc,rtol,xi,xir);
    break;
  case mMaaren:
    do_maaren(fp,eel,pts_nm,nrep);
    break;
  case mAB1:
    fprintf(fp, "#\n# Table AB1: ndisp=%d nrep=%d\n#\n",ndisp,nrep);
    do_AB1(fp,eel,pts_nm,ndisp,nrep);
    break;
  case mLjc:
    fprintf(fp, "#\n# Table LJC(12-6-1): rc=%g, rtol=%g\n#\n",rc,rtol);
    do_ljc(fp,eel,pts_nm,rc,rtol);
    break;
  case mHard_Wall:
    do_hard(fp,pts_nm,efac,delta);
    break;
  default:
    gmx_fatal(FARGS,"Model %s not supported yet",model[0]);
  }  
  if ((m != mGuillot2001a)) 
    gmx_ffclose(fp);
  
  gmx_thanx(stdout);
  
  return 0;
}
