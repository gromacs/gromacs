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
 * S  C  A  M  O  R  G
 */
static char *SRCID_clincs_c = "$Id$";

#include <math.h>
#include "main.h"
#include "update.h"
#include "physics.h"
#include "vec.h"

void cconerr(real *max,real *rms,int *imax,rvec *xprime,
	     int ncons,int *bla1,int *bla2,real *bllen)
     
{
  real      len,d,ma,ms,tmp0,tmp1,tmp2;
  int       b,i,j,im;
  
  ma=0;
  ms=0;
  im=0;
  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    tmp0=xprime[i][0]-xprime[j][0];
    tmp1=xprime[i][1]-xprime[j][1];
    tmp2=xprime[i][2]-xprime[j][2];
    len=sqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2);
    d=fabs(len/bllen[b]-1);
    if (d > ma) {
      ma=d;
      im=b;
    }
    ms=ms+d*d;
  }
  *max=ma;
  *rms=sqrt(ms/ncons);
  *imax=im;
}

void lincs_warning(rvec *x,rvec *xprime,
		   int ncons,int *bla1,int *bla2,real *bllen,real wangle)
{
  int b,i,j;
  rvec v0,v1;
  real wfac,d0,d1,cosine;
  char buf[STRLEN];
  
  wfac=cos(DEG2RAD*wangle);
  
  sprintf(buf,"bonds that rotated more than %g degrees:\n"
	  " atom 1 atom 2  angle  previous, current, constraint length\n",
	  wangle);
  fprintf(stderr,buf);
  fprintf(stdlog,buf); 

  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    rvec_sub(x[i],x[j],v0);
    rvec_sub(xprime[i],xprime[j],v1);
    d0=norm(v0);
    d1=norm(v1);
    cosine=iprod(v0,v1)/(d0*d1);
    if (cosine<wfac) {
      sprintf(buf," %6d %6d  %5.1f  %8.4f %8.4f    %8.4f\n",
	      i+1,j+1,RAD2DEG*acos(cosine),d0,d1,bllen[b]);
      fprintf(stderr,buf);
      fprintf(stdlog,buf);
    }
  }
}
