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

void calc_h_pos(int nht,int nh[],int na[],rvec x[])
{
#define alfaH   (DEG2RAD*109.5)
#define alfaCOM (DEG2RAD*117)
#define alfaCO  (DEG2RAD*121)
#define alfaCOA (DEG2RAD*115)

#define distH   0.1
#define distO   0.123
#define distOA  0.125
#define distOM  0.136

  rvec sa,sb,sij;
  int  ai,aj,ak,al;
  int  h1,h2,h3;
  real s6,rij,ra,rb,xh;
  int  m;
  
  if ((nht < 1) || (nht > 9))
    fatal_error(0,"Invalid argument (%d) for nht in routine genh\n",nht);
  
  ai  = na[0];
  aj  = na[1];
  ak  = na[2];
  al  = na[3];
  h1  = nh[0];
  h2  = nh[1];
  h3  = nh[2];

  s6=0.5*sqrt(3.e0);

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
      x[h1][m] = x[ai][m]+distH*sa[m]/ra;

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
      x[h1][m] = x[ai][m]+distH*sin(alfaH)*sb[m]-distH*cos(alfaH)*sij[m];
    }
    break;
  case 3:
    for(m=0; (m<DIM); m++) {
      x[h1][m] = x[ai][m]-distH*sin(alfaH)*sb[m]-distH*cos(alfaH)*sij[m];
      x[h2][m] = x[ai][m]+distH*sin(alfaH)*sb[m]-distH*cos(alfaH)*sij[m];
    }
    break;
  case 4:
    for(m=0; (m<DIM); m++) {
      x[h1][m] = x[ai][m]+distH*sin(alfaH)*sb[m]-distH*cos(alfaH)*sij[m];
      x[h2][m] = x[ai][m]-distH*sin(alfaH)*0.5*sb[m]+distH*sin(alfaH)*s6*sa[m]-distH*cos(alfaH)*sij[m];
      if (h3 != -1) 
	x[h3][m] = x[ai][m]-distH*sin(alfaH)*0.5*sb[m]-distH*sin(alfaH)*s6*sa[m]-distH*cos(alfaH)*sij[m];
    }
    break;
  case 5: {
    real center;
    rvec dxc;
    
    for(m=0; (m<DIM); m++) {
      center=(x[aj][m]+x[ak][m]+x[al][m])/3.0;
      dxc[m]=x[ai][m]-center;
    }
    center=norm(dxc);
    for(m=0; (m<DIM); m++)
      x[h1][m]=x[ai][m]+dxc[m]*distH/center;
    break;
  }
  case 6: {
    rvec BB,CC1,CC2,NN;
    real bb,nn;
    
    for(m=0; (m<DIM); m++) 
      BB[m]=x[ai][m]-0.5*(x[aj][m]+x[ak][m]);
    bb=norm(BB);

    rvec_sub(x[ai],x[aj],CC1);
    rvec_sub(x[ai],x[ak],CC2);
    oprod(CC1,CC2,NN);
    nn=norm(NN);
    
    for(m=0; (m<DIM); m++) {
      x[h1][m]=x[ai][m]+distH*(cos(alfaH/2.0)*BB[m]/bb+
			       sin(alfaH/2.0)*NN[m]/nn);
      x[h2][m]=x[ai][m]+distH*(cos(alfaH/2.0)*BB[m]/bb-
			       sin(alfaH/2.0)*NN[m]/nn);
    }
    break;
  }
  case 7:
    gen_waterhydrogen(&(x[ai]));
    break;
  case 8: {
    for(m=0; (m<DIM); m++) {
      x[h1][m] = x[ai][m]-distOM*sin(alfaCOM)*sb[m]-distOM*cos(alfaCOM)*sij[m];
      x[h2][m] = x[ai][m]+distOM*sin(alfaCOM)*sb[m]-distOM*cos(alfaCOM)*sij[m];
    }
    break;
  }
  case 9: {
    int na2[4]; /* i,j,k,l   */
    int nh2[3]; /* new atoms */

    /* first add two oxygens */
    for(m=0; (m<DIM); m++) {
      x[h1][m] = x[ai][m]-distO *sin(alfaCO )*sb[m]-distO *cos(alfaCO )*sij[m];
      x[h2][m] = x[ai][m]+distOA*sin(alfaCOA)*sb[m]-distOA*cos(alfaCOA)*sij[m];
    }
    
    /* now use rule 2 to add hydrogen to 2nd oxygen */
    na2[0]=nh[1]; /* new i = n' */
    na2[1]=na[0]; /* new j = i  */
    na2[2]=na[1]; /* new k = j  */
    na2[3]=na[2]; /* new l = k, not used */
    nh2[0]=nh[2]; /* n' is third new atom */
    calc_h_pos(2,nh2,na2,x);
    
    break;
  }
  }
}

