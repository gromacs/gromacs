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
 * Gyas ROwers Mature At Cryogenic Speed
 */
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

void calc_h_pos(int nht,int nh[],int na[],real d,real alfa,rvec x[])
{
  rvec sa,sb,sij;
  int  ai,aj,ak,al;
  int  h1,h2,h3;
  real s6,sad,cad,rij,ra,rb,xh;
  int  m;
  
  if ((nht < 1) || (nht > 7))
    fatal_error(0,"Invalid argument (%d) for nht in routine genh\n",nht);
  
  ai  = na[0];
  aj  = na[1];
  ak  = na[2];
  h1  = nh[0];
  h2  = nh[1];
  h3  = nh[2];

  s6=0.5*sqrt(3.e0);
  sad=d* sin(alfa);
  cad=d* cos(alfa);

  /* construct one planar hydrogen (peptide,rings) */
  if (nht == 1) {
    rij = 0.e0;
    rb  = 0.e0;
    for(m=0; (m<DIM); m++) {
      sij[m] = x[ai][m]-x[aj][m];
      sb[m]  = x[ai][m]-x[ak][m];
      rij   += sqr(sij[m]);
      rb    += sqr(sb[m]);
    }
    rij = sqrt(rij);
    rb  = sqrt(rb);
    ra  = 0.e0;
    for(m=0; (m<DIM); m++) {
      sa[m] = sij[m]/rij+sb[m]/rb;
      ra   += sqr(sa[m]);
    }
    ra = sqrt(ra);
    for(m=0; (m<DIM); m++)
      x[h1][m] = x[ai][m]+d*sa[m]/ra;

    return;
  }
  
  /* construct one, two or three dihedral hydrogens */
  rij = 0.e0;
  for(m=0; (m<DIM); m++) {
    xh     = x[aj][m];
    sij[m] = x[ai][m]-xh;
    sb[m]  = xh-x[ak][m];
    rij   += sqr(sij[m]);
  }
  
  rij = sqrt(rij);
  sa[XX] = sij[YY]*sb[ZZ]-sij[ZZ]*sb[YY];
  sa[YY] = sij[ZZ]*sb[XX]-sij[XX]*sb[ZZ];
  sa[ZZ] = sij[XX]*sb[YY]-sij[YY]*sb[XX];
  ra = 0.e0;
  for(m=0; (m<DIM); m++) {
    sij[m] = sij[m]/rij;
    ra    += sqr(sa[m]);
  }
  ra = sqrt(ra);
  for(m=0; (m<DIM); m++) 
    sa[m] = sa[m]/ra;

  sb[XX] = sa[YY]*sij[ZZ]-sa[ZZ]*sij[YY];
  sb[YY] = sa[ZZ]*sij[XX]-sa[XX]*sij[ZZ];
  sb[ZZ] = sa[XX]*sij[YY]-sa[YY]*sij[XX];

  switch (nht) {
  case 2:
    for(m=0; (m<DIM); m++) {
      x[h1][m] = x[ai][m]+sad*sb[m]-cad*sij[m];
    }
    break;
  case 3:
    for(m=0; (m<DIM); m++) {
      x[h1][m] = x[ai][m]-sad*sb[m]-cad*sij[m];
      x[h2][m] = x[ai][m]+sad*sb[m]-cad*sij[m];
    }
    break;
  case 4:
    for(m=0; (m<DIM); m++) {
      x[h1][m] = x[ai][m]+sad*sb[m]-cad*sij[m];
      x[h2][m] = x[ai][m]-sad*0.5*sb[m]+sad*s6*sa[m]-cad*sij[m];
      if (h3 != -1) 
	x[h3][m] = x[ai][m]-sad*0.5*sb[m]-sad*s6*sa[m]-cad*sij[m];
    }
    break;
  case 5: {
    real center;
    rvec dxc;
    
    al  = na[3];
    for(m=0; (m<DIM); m++) {
      center=(x[aj][m]+x[ak][m]+x[al][m])/3.0;
      dxc[m]=x[ai][m]-center;
    }
    center=norm(dxc);
    for(m=0; (m<DIM); m++)
      x[h1][m]=x[ai][m]+dxc[m]*d/center;
    break;
  }
  case 6: {
    rvec BB,CC1,CC2,NN;
    real bb,nn,hoek;
    
    for(m=0; (m<DIM); m++) 
      BB[m]=x[ai][m]-0.5*(x[aj][m]+x[ak][m]);
    bb=norm(BB);

    rvec_sub(x[ai],x[aj],CC1);
    rvec_sub(x[ai],x[ak],CC2);
    oprod(CC1,CC2,NN);
    nn=norm(NN);
    
    hoek=109.5*DEG2RAD/2.0;
    
    for(m=0; (m<DIM); m++) {
      x[h1][m]=x[ai][m]+d*(cos(hoek)*BB[m]/bb+sin(hoek)*NN[m]/nn);
      x[h2][m]=x[ai][m]+d*(cos(hoek)*BB[m]/bb-sin(hoek)*NN[m]/nn);
    }
    break;
  }
  case 7:
    gen_waterhydrogen(&(x[ai]));
    break;
  }
}

