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
static char *SRCID_calch_c = "$Id$";
#include "macros.h"
#include "calch.h"
#include "maths.h"
#include "vec.h"
#include "physics.h"
	
void gen_waterhydrogen(rvec x[3])
{
#define AA 0.081649
#define BB 0.0
#define CC 0.0577350
  const  rvec   XH1[6] = {
    { AA,     BB,     CC },
    { AA,     BB,     CC },
    { AA,     BB,     CC },
    { -AA,    BB,     CC },
    { -AA,    BB,     CC },
    { BB,     AA,    -CC }
  };
  const  rvec   XH2[6] = {
    { -AA,   BB,   CC },
    { BB,    AA,  -CC },
    { BB,   -AA,  -CC },
    { BB,    AA,  -CC },
    { BB,   -AA,  -CC },
    { BB,   -AA,  -CC }
  };
#undef AA
#undef BB
#undef CC
  static int l=0;
  int        m;
  
  /* This was copied from Gromos */
  for(m=0; (m<DIM); m++) {
    x[1][m]=x[0][m]+XH1[l][m];
    x[2][m]=x[0][m]+XH2[l][m];
  }
      
  l=(l+1) % 6;
}

void calc_h_pos(int nht,int nh[],int a[],rvec x[])
{
#define alfaH   (DEG2RAD*109.5)
#define distH   0.1

#define alfaCOM (DEG2RAD*117)
#define alfaCO  (DEG2RAD*121)
#define alfaCOA (DEG2RAD*115)

#define distO   0.123
#define distOA  0.125
#define distOM  0.136

  rvec sa,sb,sij;
  real s6,rij,ra,rb,xh;
  int  d;
  
  if ((nht < 1) || (nht > 9))
    fatal_error(0,"Invalid argument (%d) for nht in routine genh\n",nht);
  
  /* from macros.h:  
     #define AI a[0]
     #define AJ a[1]
     #define AK a[2]
     #define AL a[3] */
     
#define H1 nh[0]
#define H2 nh[1]
#define H3 nh[2]

  s6=0.5*sqrt(3.e0);

  /* construct one planar hydrogen (peptide,rings) */
  if (nht == 1) {
    rij = 0.e0;
    rb  = 0.e0;
    for(d=0; (d<DIM); d++) {
      sij[d] = x[AI][d]-x[AJ][d];
      sb[d]  = x[AI][d]-x[AK][d];
      rij   += sqr(sij[d]);
      rb    += sqr(sb[d]);
    }
    rij = sqrt(rij);
    rb  = sqrt(rb);
    ra  = 0.e0;
    for(d=0; (d<DIM); d++) {
      sa[d] = sij[d]/rij+sb[d]/rb;
      ra   += sqr(sa[d]);
    }
    ra = sqrt(ra);
    for(d=0; (d<DIM); d++)
      x[H1][d] = x[AI][d]+distH*sa[d]/ra;

    return;
  }
  
  /* construct one, two or three dihedral hydrogens */
  if ((nht != 5) && (nht !=6) && (nht != 7)) {
    rij = 0.e0;
    for(d=0; (d<DIM); d++) {
      xh     = x[AJ][d];
      sij[d] = x[AI][d]-xh;
      sb[d]  = xh-x[AK][d];
      rij   += sqr(sij[d]);
    }
    rij = sqrt(rij);
    sa[XX] = sij[YY]*sb[ZZ]-sij[ZZ]*sb[YY];
    sa[YY] = sij[ZZ]*sb[XX]-sij[XX]*sb[ZZ];
    sa[ZZ] = sij[XX]*sb[YY]-sij[YY]*sb[XX];
    ra = 0.e0;
    for(d=0; (d<DIM); d++) {
      sij[d] = sij[d]/rij;
      ra    += sqr(sa[d]);
    }
    ra = sqrt(ra);
    for(d=0; (d<DIM); d++) 
      sa[d] = sa[d]/ra;
    
    sb[XX] = sa[YY]*sij[ZZ]-sa[ZZ]*sij[YY];
    sb[YY] = sa[ZZ]*sij[XX]-sa[XX]*sij[ZZ];
    sb[ZZ] = sa[XX]*sij[YY]-sa[YY]*sij[XX];
  }

  switch (nht) {
  case 2: /* one single hydrogen, e.g. hydroxyl */
    for(d=0; (d<DIM); d++) {
      x[H1][d] = x[AI][d]+distH*sin(alfaH)*sb[d]-distH*cos(alfaH)*sij[d];
    }
    break;
  case 3: /* two planar hydrogens, e.g. -NH2 */
    for(d=0; (d<DIM); d++) {
      x[H1][d] = x[AI][d]-distH*sin(alfaH)*sb[d]-distH*cos(alfaH)*sij[d];
      x[H2][d] = x[AI][d]+distH*sin(alfaH)*sb[d]-distH*cos(alfaH)*sij[d];
    }
    break;
  case 4: /* two or three tetrahedral hydrogens, e.g. -CH3 */
    for(d=0; (d<DIM); d++) {
      x[H1][d] = x[AI][d]+distH*sin(alfaH)*sb[d]-distH*cos(alfaH)*sij[d];
      x[H2][d] = ( x[AI][d] 
		   - distH*sin(alfaH)*0.5*sb[d]
		   + distH*sin(alfaH)*s6*sa[d]
		   - distH*cos(alfaH)*sij[d] );
      if (H3 != -1) 
	x[H3][d] = ( x[AI][d]
		     - distH*sin(alfaH)*0.5*sb[d]
		     - distH*sin(alfaH)*s6*sa[d]
		     - distH*cos(alfaH)*sij[d] );
    }
    break;
  case 5: { /* one tetrahedral hydrogen, e.g. C3CH */
    real center;
    rvec dxc;
    
    for(d=0; (d<DIM); d++) {
      center=(x[AJ][d]+x[AK][d]+x[AL][d])/3.0;
      dxc[d]=x[AI][d]-center;
    }
    center=norm(dxc);
    for(d=0; (d<DIM); d++)
      x[H1][d]=x[AI][d]+dxc[d]*distH/center;
    break;
  }
  case 6: { /* two tetrahedral hydrogens, e.g. C-CH2-C */
    rvec BB,CC1,CC2,NN;
    real bb,nn;
    
    for(d=0; (d<DIM); d++) 
      BB[d]=x[AI][d]-0.5*(x[AJ][d]+x[AK][d]);
    bb=norm(BB);

    rvec_sub(x[AI],x[AJ],CC1);
    rvec_sub(x[AI],x[AK],CC2);
    oprod(CC1,CC2,NN);
    nn=norm(NN);
    
    for(d=0; (d<DIM); d++) {
      x[H1][d]=x[AI][d]+distH*(cos(alfaH/2.0)*BB[d]/bb+
			       sin(alfaH/2.0)*NN[d]/nn);
      x[H2][d]=x[AI][d]+distH*(cos(alfaH/2.0)*BB[d]/bb-
			       sin(alfaH/2.0)*NN[d]/nn);
    }
    break;
  }
  case 7: /* two water hydrogens */
    gen_waterhydrogen(&(x[AI]));
    break;
  case 8: { /* two carboxyl oxygens, -COO- */
    for(d=0; (d<DIM); d++) {
      x[H1][d] = x[AI][d]-distOM*sin(alfaCOM)*sb[d]-distOM*cos(alfaCOM)*sij[d];
      x[H2][d] = x[AI][d]+distOM*sin(alfaCOM)*sb[d]-distOM*cos(alfaCOM)*sij[d];
    }
    break;
  }
  case 9: { /* carboxyl oxygens and hydrogen, -COOH */
    int na2[4]; /* i,j,k,l   */
    int nh2[3]; /* new atoms */

    /* first add two oxygens */
    for(d=0; (d<DIM); d++) {
      x[H1][d] = x[AI][d]-distO *sin(alfaCO )*sb[d]-distO *cos(alfaCO )*sij[d];
      x[H2][d] = x[AI][d]+distOA*sin(alfaCOA)*sb[d]-distOA*cos(alfaCOA)*sij[d];
    }
    
    /* now use rule 2 to add hydrogen to 2nd oxygen */
    na2[0]=H2; /* new i = n' */
    na2[1]=AI; /* new j = i  */
    na2[2]=AJ; /* new k = j  */
    na2[3]=AK; /* new l = k, not used */
    nh2[0]=H3; /* n' is third new atom */
    calc_h_pos(2,nh2,na2,x);
    
    break;
  }
  } /* end switch */
}
