/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_csettle_c = "$Id$";

#include <math.h>
#include <stdio.h>
#include "vec.h"
#include "constr.h"
#include "fatal.h"

#ifdef DEBUG
static void check_cons(FILE *fp,char *title,real x[],int OW1,int HW2,int HW3)
{
  rvec dOH1,dOH2,dHH;
  int  m;
  
  for(m=0; (m<DIM); m++) {
    dOH1[m]=x[OW1+m]-x[HW2+m];
    dOH2[m]=x[OW1+m]-x[HW3+m];
    dHH[m]=x[HW2+m]-x[HW3+m];
  }
  fprintf(fp,"%10s, OW1=%3d, HW2=%3d, HW3=%3d,  dOH1: %8.3f, dOH2: %8.3f, dHH: %8.3f\n",
	  title,OW1/DIM,HW2/DIM,HW3/DIM,norm(dOH1),norm(dOH2),norm(dHH));
}
#endif

void csettle(FILE *fp,int nshake, int owptr[],real b4[], real after[],
	     real dOH,real dHH,real mO,real mH,int *error)
{
  /* ***************************************************************** */
  /*                                                               ** */
  /*    Subroutine : setlep - reset positions of TIP3P waters      ** */
  /*    Author : Shuichi Miyamoto                                  ** */
  /*    Date of last update : Oct. 1, 1992                         ** */
  /*                                                               ** */
  /*    Reference for the SETTLE algorithm                         ** */
  /*           S. Miyamoto et al., J. Comp. Chem., 13, 952 (1992). ** */
  /*                                                               ** */
  /* ***************************************************************** */
  
  /* Initialized data */
  static bool bFirst=TRUE;
  static real wo,wh,wohh;
  static real ra,rb,rc,rc2,rone;
#ifdef DEBUG_PRES
  static int step = 0;
#endif
  
    
  /* Local variables */
  real gama, beta, alpa, xcom, ycom, zcom, al2be2, tmp, tmp2;
  real axlng, aylng, azlng, trns11, trns21, trns31, trns12, trns22, 
    trns32, trns13, trns23, trns33, cosphi, costhe, sinphi, sinthe, 
    cospsi, xaksxd, yaksxd, xakszd, yakszd, zakszd, zaksxd, xaksyd, 
    xb0, yb0, zb0, xc0, yc0, zc0, xa1;
  real ya1, za1, xb1, yb1;
  real zb1, xc1, yc1, zc1, yaksyd, zaksyd, sinpsi, xa3, ya3, za3, 
    xb3, yb3, zb3, xc3, yc3, zc3, xb0d, yb0d, xc0d, yc0d, 
    za1d, xb1d, yb1d, zb1d, xc1d, yc1d, zc1d, ya2d, xb2d, yb2d, yc2d, 
    xa3d, ya3d, za3d, xb3d, yb3d, zb3d, xc3d, yc3d, zc3d;
  real t1,t2;
  
  int i, ow1, hw2, hw3;

  *error=-1;
  if (bFirst) {
    fprintf(fp,"Going to use C-settle (%d waters)\n",nshake);
    wo     = mO;
    wh     = mH;
    wohh   = mO+2.0*mH;
    rc     = dHH/2.0;
    ra     = 2.0*wh*sqrt(dOH*dOH-rc*rc)/wohh;
    rb     = sqrt(dOH*dOH-rc*rc)-ra;
    rc2    = dHH;
    rone   = 1.0;

    wo    /= wohh;
    wh    /= wohh;
        
    fprintf(fp,"wo = %g, wh =%g, wohh = %g, rc = %g, ra = %g\n",
	    wo,wh,wohh,rc,ra);
    fprintf(fp,"rb = %g, rc2 = %g, rone = %g, dHH = %g, dOH = %g\n",
	    rb,rc2,rone,dHH,dOH);
	
    bFirst = FALSE;
  }
#ifdef PRAGMAS
#pragma ivdep
#endif
  for (i = 0; i < nshake; ++i) {
    /*    --- Step1  A1' ---      */
    ow1 = owptr[i] * 3;
    hw2 = ow1 + 3;
    hw3 = ow1 + 6;
    xb0 = b4[hw2    ] - b4[ow1];
    yb0 = b4[hw2 + 1] - b4[ow1 + 1];
    zb0 = b4[hw2 + 2] - b4[ow1 + 2];
    xc0 = b4[hw3    ] - b4[ow1];
    yc0 = b4[hw3 + 1] - b4[ow1 + 1];
    zc0 = b4[hw3 + 2] - b4[ow1 + 2];
    /* 6 flops */
    
    xcom = (after[ow1    ] * wo + (after[hw2    ] + after[hw3    ]) * wh);
    ycom = (after[ow1 + 1] * wo + (after[hw2 + 1] + after[hw3 + 1]) * wh);
    zcom = (after[ow1 + 2] * wo + (after[hw2 + 2] + after[hw3 + 2]) * wh);
    /* 12 flops */
    
    xa1 = after[ow1    ] - xcom;
    ya1 = after[ow1 + 1] - ycom;
    za1 = after[ow1 + 2] - zcom;
    xb1 = after[hw2    ] - xcom;
    yb1 = after[hw2 + 1] - ycom;
    zb1 = after[hw2 + 2] - zcom;
    xc1 = after[hw3    ] - xcom;
    yc1 = after[hw3 + 1] - ycom;
    zc1 = after[hw3 + 2] - zcom;
    /* 9 flops */
    
    xakszd = yb0 * zc0 - zb0 * yc0;
    yakszd = zb0 * xc0 - xb0 * zc0;
    zakszd = xb0 * yc0 - yb0 * xc0;
    xaksxd = ya1 * zakszd - za1 * yakszd;
    yaksxd = za1 * xakszd - xa1 * zakszd;
    zaksxd = xa1 * yakszd - ya1 * xakszd;
    xaksyd = yakszd * zaksxd - zakszd * yaksxd;
    yaksyd = zakszd * xaksxd - xakszd * zaksxd;
    zaksyd = xakszd * yaksxd - yakszd * xaksxd;
    /* 27 flops */
    
    axlng = invsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
    aylng = invsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
    azlng = invsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);
      
    trns11 = xaksxd * axlng;
    trns21 = yaksxd * axlng;
    trns31 = zaksxd * axlng;
    trns12 = xaksyd * aylng;
    trns22 = yaksyd * aylng;
    trns32 = zaksyd * aylng;
    trns13 = xakszd * azlng;
    trns23 = yakszd * azlng;
    trns33 = zakszd * azlng;
    /* 24 flops */
    
    xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
    yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
    xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
    yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;
    /*
    xa1d = trns11 * xa1 + trns21 * ya1 + trns31 * za1;
    ya1d = trns12 * xa1 + trns22 * ya1 + trns32 * za1;
    */
    za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1;
    xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
    yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
    zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
    xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
    yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
    zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;
    /* 65 flops */
        
    sinphi = za1d / ra;
    tmp    = rone - sinphi * sinphi;
    if (tmp <= 0) {
      *error = i;
      cosphi = 0;
    }
    else
      cosphi = sqrt(tmp);
    sinpsi = (zb1d - zc1d) / (rc2 * cosphi);
    tmp2   = rone - sinpsi * sinpsi;
    if (tmp2 <= 0) {
      *error = i;
      cospsi = 0;
    }
    else
      cospsi = sqrt(tmp2);
    /* 46 flops */
    
    ya2d =  ra * cosphi;
    xb2d = -rc * cospsi;
    t1   = -rb * cosphi;
    t2   =  rc * sinpsi * sinphi;
    yb2d =  t1 - t2;
    yc2d =  t1 + t2;
    /* 7 flops */
        
    /*     --- Step3  al,be,ga 		      --- */
    alpa   = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
    beta   = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
    gama   = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
    al2be2 = alpa * alpa + beta * beta;
    sinthe = (alpa * gama - beta * sqrt(al2be2 - gama * gama)) / al2be2;
    /* 47 flops */
    
    /*  --- Step4  A3' --- */
    costhe = sqrt(rone - sinthe * sinthe);
    xa3d = -ya2d * sinthe;
    ya3d = ya2d * costhe;
    za3d = za1d;
    xb3d = xb2d * costhe - yb2d * sinthe;
    yb3d = xb2d * sinthe + yb2d * costhe;
    zb3d = zb1d;
    xc3d = -xb2d * costhe - yc2d * sinthe;
    yc3d = -xb2d * sinthe + yc2d * costhe;
    zc3d = zc1d;
    /* 26 flops */
    
    /*    --- Step5  A3 --- */
    xa3 = trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
    ya3 = trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
    za3 = trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
    xb3 = trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
    yb3 = trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
    zb3 = trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
    xc3 = trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
    yc3 = trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
    zc3 = trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;
    /* 45 flops */
    
    after[ow1] = xcom + xa3;
    after[ow1 + 1] = ycom + ya3;
    after[ow1 + 2] = zcom + za3;
    after[hw2] = xcom + xb3;
    after[hw2 + 1] = ycom + yb3;
    after[hw2 + 2] = zcom + zb3;
    after[hw3] = xcom + xc3;
    after[hw3 + 1] = ycom + yc3;
    after[hw3 + 2] = zcom + zc3;
    /* 9 flops */
#ifdef DEBUG
    check_cons(fp,"settle",after,ow1,hw2,hw3);
#endif
  }
}

