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
static char *SRCID_inloopc_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "main.h"
#include "futil.h"
#include "typedefs.h"
#include "vec.h"
#include "fatal.h"
#include "inloop.h"

void c_coul(real ix,real iy,real iz,real qi,
	    real pos[],int nj,t_nl_j jjnr[],
	    real charge[],
	    real faction[],real fip[3],
	    real *Vc)
{
  int     k,jnr,j3;
  real    fX,fY,fZ;
  real    rijX,rijY,rijZ;
  real    fijscal,rsq,vijcoul,vctot;
  real    rinv1,rinv2;
  real    fjx,fjy,fjz;
  real    tx,ty,tz;

  fX=0;
  fY=0;
  fZ=0; 

  vctot = 0;
  
  /* Cray or PGC compiler directive: ignore vector dependency */
#pragma ivdep
  for(k=0; (k<nj); k++) {
    jnr            = jjnr[k];
    j3             = 3*jnr;
    fjx            = faction[j3+XX];
    fjy            = faction[j3+YY];
    fjz            = faction[j3+ZZ];
    rijX           = ix - pos[j3+XX];
    rijY           = iy - pos[j3+YY];
    rijZ           = iz - pos[j3+ZZ];
    rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ);
    
    rinv1          = invsqrt(rsq);
    vijcoul        = qi*charge[jnr]*rinv1;
    rinv2          = rinv1*rinv1;
    vctot          = vctot+vijcoul;
    fijscal        = vijcoul*rinv2;
    tx             = rijX*fijscal;
    ty             = rijY*fijscal;
    tz             = rijZ*fijscal;
    fX             = fX+tx;
    fY             = fY+ty;
    fZ             = fZ+tz;
    faction[j3+XX] = fjx-tx;
    faction[j3+YY] = fjy-ty;
    faction[j3+ZZ] = fjz-tz;
  }
  
  fip[XX]=fX;
  fip[YY]=fY;
  fip[ZZ]=fZ;
  *Vc += vctot;
}

void c_ljc(real ix,real iy,real iz,real qi,
	   real pos[],int nj,int type[],t_nl_j jjnr[],
	   real charge[],real nbfp[],
	   real faction[],real fip[],
	   real *Vc,real *Vnb)
{      
  const  real twelve=12.0;
  const  real six=6.0;
  int    k,jnr,j3,tj;
  real   fX,fY,fZ;
  real   rijX,rijY,rijZ;
  real   fijscal,vijcoul,vctot,vnbtot;
  real   rinv1,rinv2,rinv6;
  real   fjx,fjy,fjz;
  real   tx,ty,tz,vnb6,vnb12;
  
  fX     = 0;
  fY     = 0;
  fZ     = 0;
  vctot  = 0;
  vnbtot = 0;
  
  /* Cray or PGC compiler directive: ignore vector dependency */
#pragma ivdep
  for(k=0; (k<nj); k++) {
    jnr            = jjnr[k];
    j3             = 3*jnr;
    rijX           = ix - pos[j3+XX];
    rijY           = iy - pos[j3+YY];
    rijZ           = iz - pos[j3+ZZ];
    
    rinv1          = invsqrt((rijX*rijX)+(rijY*rijY)+(rijZ*rijZ));
    rinv2          = rinv1*rinv1;
    rinv6          = rinv2*rinv2*rinv2;
    
    tj             = 2*type[jnr];
    vnb6           = nbfp[tj]   * rinv6;
    vnb12          = nbfp[tj+1] * rinv6*rinv6;
    vijcoul        = qi*charge[jnr]*rinv1;
    vctot          = vctot+vijcoul;
    vnbtot         = vnbtot+vnb12-vnb6;
    fijscal        = (twelve*vnb12-six*vnb6+vijcoul)*rinv2;
    
    fjx            = faction[j3+XX];
    tx             = rijX*fijscal;
    fX             = fX + tx;
    faction[j3+XX] = fjx - tx;
    fjy            = faction[j3+YY];
    ty             = rijY*fijscal;
    fY             = fY + ty;
    faction[j3+YY] = fjy - ty;
    fjz            = faction[j3+ZZ];
    tz             = rijZ*fijscal;
    fZ             = fZ + tz;
    faction[j3+ZZ] = fjz - tz;
  }

  fip[XX]  = fX;
  fip[YY]  = fY;
  fip[ZZ]  = fZ;
  *Vc     += vctot;
  *Vnb    += vnbtot;
}

void c_bham(real ix,real iy,real iz,real qi,
	    real pos[],int nj,int type[],t_nl_j jjnr[],
	    real charge[],real nbfp[],
	    real faction[],real fip[],
	    real *Vc,real *Vnb)
{      
  const  real six=6.0;
  int    k,jnr,j3,tj;
  real   fX,fY,fZ;
  real   rijX,rijY,rijZ;
  real   fijscal,vijcoul,vctot,vnbtot;
  real   r1,r_sq,rinv1,rinv2,rinv6;
  real   fjx,fjy,fjz;
  real   tx,ty,tz,vnb6,vnbexp,br;
  
  fX=0;
  fY=0;
  fZ=0;
  vctot  = 0;
  vnbtot = 0;
  
#ifdef _860_  
  fatal_error(0,"Sorry no Buckingham");
#else

  /* Cray or PGC compiler directive: ignore vector dependency */
#pragma ivdep
  for(k=0; (k<nj); k++) {
    jnr            = jjnr[k];
    j3             = 3*jnr;
    rijX           = ix - pos[j3+XX];
    rijY           = iy - pos[j3+YY];
    rijZ           = iz - pos[j3+ZZ];
    r_sq           = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ);
    rinv1          = invsqrt(r_sq);
    rinv2          = rinv1*rinv1;
    rinv6          = rinv2*rinv2*rinv2;
    r1             = r_sq*rinv1;
    
    tj             = 3*type[jnr];
    br             = nbfp[tj+1]*r1;
    vnbexp         = nbfp[tj]*exp(-br);
    vnb6           = nbfp[tj+2]*rinv6;
    vijcoul        = qi*charge[jnr]*rinv1;
    
    vctot          = vctot+vijcoul;
    vnbtot         = vnbtot+vnbexp-vnb6;
    fijscal        = (br*vnbexp-six*vnb6+vijcoul)*rinv2;
    
    fjx            = faction[j3+XX];
    tx             = rijX*fijscal;
    fX             = fX + tx;
    faction[j3+XX] = fjx - tx;
    fjy            = faction[j3+YY];
    ty             = rijY*fijscal;
    fY             = fY + ty;
    faction[j3+YY] = fjy - ty;
    fjz            = faction[j3+ZZ];
    tz             = rijZ*fijscal;
    fZ             = fZ + tz;
    faction[j3+ZZ] = fjz - tz;
  }
#endif

  fip[XX]=fX;
  fip[YY]=fY;
  fip[ZZ]=fZ;
  *Vc     += vctot;
  *Vnb    += vnbtot;
}

void c_water(int i0,real xw[],real fudge,
	     real pos[],int nj,int type[],t_nl_j jjnr[],
	     real charge[],real nbfp[],
	     real faction[],real fw[],
	     real *Vc,real *Vnb)
{      
  const  real twelve=12.0;
  const  real six=6.0;
  int    k,jnr,j3,tj;
  real   rinv6;
  real   fxJ,fyJ,fzJ;
  real   vnb6,vnb12,vctot,vnbtot;
  real   fxO,fyO,fzO,fxH1,fyH1,fzH1,fxH2,fyH2,fzH2;
  real   ixO,iyO,izO,ixH1,iyH1,izH1,ixH2,iyH2,izH2;
  real   dxO,dyO,dzO,dxH1,dyH1,dzH1,dxH2,dyH2,dzH2;
  real   txO,tyO,tzO,txH1,tyH1,tzH1,txH2,tyH2,tzH2;
  real   rinv1O,rinv1H1,rinv1H2,rinv2O,rinv2H1,rinv2H2;
  real   vcO,vcH1,vcH2,fsO,fsH1,fsH2;
  real   qO,qH,qj;
  real   jx,jy,jz;
  
  /* initiate local force variables */
  fxO=fyO=fzO=fxH1=fyH1=fzH1=fxH2=fyH2=fzH2=0.0;
  
  /* Copy i-particles into temps */
  ixO    = xw[0];
  iyO    = xw[1];
  izO    = xw[2];
  ixH1   = xw[3];
  iyH1   = xw[4];
  izH1   = xw[5];
  ixH2   = xw[6];
  iyH2   = xw[7];
  izH2   = xw[8];
  
  vctot  = 0;
  vnbtot = 0;
  
  qO     = charge[i0  ]*fudge;
  qH     = charge[i0+1]*fudge;
    
  /* Currently this loop costs 93 floating point operations when
   * invsqrt = 5 (as when using GROMACS version).
   * On the i860 we use (40 MHz) the loop takes 214 cycles, which gives
   * a performance of 17.3 MFlop/s
   */
#ifdef _860_
  fatal_error(0,"Sorry, no waterloop here");
#else

  /* Cray or PGC compiler directive: ignore vector dependency */
#pragma ivdep
  for(k=0; (k<nj); k++) {
    jnr = jjnr[k];
    j3  = 3*jnr;
    qj  = charge[jnr];
    
    jx  = pos[j3+XX];
    jy  = pos[j3+YY];
    jz  = pos[j3+ZZ];

    /* First one is for oxygen, with LJ */    
    dxO           = ixO - jx;
    dyO           = iyO - jy;
    dzO           = izO - jz;
    
    rinv1O         = invsqrt((dxO*dxO)+(dyO*dyO)+(dzO*dzO));
    /* Now for first hydrogen */
    dxH1           = ixH1 - jx;
    dyH1           = iyH1 - jy;
    dzH1           = izH1 - jz;
    rinv1H1        = invsqrt((dxH1*dxH1)+(dyH1*dyH1)+(dzH1*dzH1));
    
    /* Now for second hydrogen */
    dxH2           = ixH2 - jx;
    dyH2           = iyH2 - jy;
    dzH2           = izH2 - jz;
    rinv1H2        = invsqrt((dxH2*dxH2)+(dyH2*dyH2)+(dzH2*dzH2));
    
    /* O block */
    rinv2O         = rinv1O*rinv1O;
    rinv6          = rinv2O*rinv2O*rinv2O;
    tj             = 2*type[jnr];
    vnb6           = nbfp[tj]   * rinv6;
    vnb12          = nbfp[tj+1] * rinv6*rinv6;
    vcO            = qO*qj*rinv1O;
    vctot          = vctot+vcO;
    vnbtot         = vnbtot+vnb12-vnb6;
    fsO            = (twelve*vnb12-six*vnb6+vcO)*rinv2O;
    
    txO            = dxO*fsO;
    fxO            = fxO + txO;
    fxJ            = faction[j3+XX] - txO;
    tyO            = dyO*fsO;
    fyO            = fyO + tyO;
    fyJ            = faction[j3+YY] - tyO;
    tzO            = dzO*fsO;
    fzO            = fzO + tzO;
    fzJ            = faction[j3+ZZ] - tzO;

    /* H1 block */    
    rinv2H1        = rinv1H1*rinv1H1;
    vcH1           = qH*qj*rinv1H1;
    vctot          = vctot+vcH1;
    fsH1           = vcH1*rinv2H1;
    
    txH1           = dxH1*fsH1;
    fxH1           = fxH1 + txH1;
    fxJ            = fxJ - txH1;
    tyH1           = dyH1*fsH1;
    fyH1           = fyH1 + tyH1;
    fyJ            = fyJ - tyH1;
    tzH1           = dzH1*fsH1;
    fzH1           = fzH1 + tzH1;
    fzJ            = fzJ - tzH1;
    
    /* H2 block */
    rinv2H2        = rinv1H2*rinv1H2;
    vcH2           = qH*qj*rinv1H2;
    vctot          = vctot+vcH2;
    fsH2           = vcH2*rinv2H2;
    
    txH2           = dxH2*fsH2;
    fxH2           = fxH2 + txH2;
    faction[j3+XX] = fxJ - txH2;
    tyH2           = dyH2*fsH2;
    fyH2           = fyH2 + tyH2;
    faction[j3+YY] = fyJ - tyH2;
    tzH2           = dzH2*fsH2;
    fzH2           = fzH2 + tzH2;
    faction[j3+ZZ] = fzJ - tzH2;
  }
#endif
  fw[0]=fxO;
  fw[1]=fyO;
  fw[2]=fzO;
  fw[3]=fxH1;
  fw[4]=fyH1;
  fw[5]=fzH1;
  fw[6]=fxH2;
  fw[7]=fyH2;
  fw[8]=fzH2;

  *Vc     += vctot;
  *Vnb    += vnbtot;
}

void c_wcoul(int i0,real xw[],real fudge,
	     real pos[],int nj,t_nl_j jjnr[],
	     real charge[],real faction[],real fw[],
	     real *Vc)
{      
  int    k,jnr,j3;
  real   fxJ,fyJ,fzJ;
  real   fxO,fyO,fzO,fxH1,fyH1,fzH1,fxH2,fyH2,fzH2;
  real   ixO,iyO,izO,ixH1,iyH1,izH1,ixH2,iyH2,izH2;
  real   dxO,dyO,dzO,dxH1,dyH1,dzH1,dxH2,dyH2,dzH2;
  real   txO,tyO,tzO,txH1,tyH1,tzH1,txH2,tyH2,tzH2;
  real   rinv1O,rinv1H1,rinv1H2,rinv2O,rinv2H1,rinv2H2;
  real   vcO,vcH1,vcH2,fsO,fsH1,fsH2,vctot;
  real   qO,qH,qj;
  real   jx,jy,jz;
  
  /* initiate local force variables */
  fxO=fyO=fzO=fxH1=fyH1=fzH1=fxH2=fyH2=fzH2=0.0;
  
  /* Copy i-particles into temps */
  ixO  = xw[0];
  iyO  = xw[1];
  izO  = xw[2];
  ixH1 = xw[3];
  iyH1 = xw[4];
  izH1 = xw[5];
  ixH2 = xw[6];
  iyH2 = xw[7];
  izH2 = xw[8];
  
  vctot= 0;
  
  qO   = charge[i0  ]*fudge;
  qH   = charge[i0+1]*fudge;
    
#ifdef _860_
  fatal_error(0,"Sorry, no waterloop here");
#else

  /* Cray or PGC compiler directive: ignore vector dependency */
#pragma ivdep
  for(k=0; (k<nj); k++) {
    jnr = jjnr[k];
    j3  = 3*jnr;
    qj  = charge[jnr];
    
    jx  = pos[j3+XX];
    jy  = pos[j3+YY];
    jz  = pos[j3+ZZ];

    /* First one is for oxygen, with LJ */    
    dxO           = ixO - jx;
    dyO           = iyO - jy;
    dzO           = izO - jz;
    
    rinv1O         = invsqrt((dxO*dxO)+(dyO*dyO)+(dzO*dzO));
    /* Now for first hydrogen */
    dxH1           = ixH1 - jx;
    dyH1           = iyH1 - jy;
    dzH1           = izH1 - jz;
    rinv1H1        = invsqrt((dxH1*dxH1)+(dyH1*dyH1)+(dzH1*dzH1));
    
    /* Now for second hydrogen */
    dxH2           = ixH2 - jx;
    dyH2           = iyH2 - jy;
    dzH2           = izH2 - jz;
    rinv1H2        = invsqrt((dxH2*dxH2)+(dyH2*dyH2)+(dzH2*dzH2));
    
    /* O block */
    rinv2O         = rinv1O*rinv1O;
    vcO            = qO*qj*rinv1O;
    vctot          = vctot+vcO;
    fsO            = vcO*rinv2O;
    
    txO            = dxO*fsO;
    fxO            = fxO + txO;
    fxJ            = faction[j3+XX] - txO;
    tyO            = dyO*fsO;
    fyO            = fyO + tyO;
    fyJ            = faction[j3+YY] - tyO;
    tzO            = dzO*fsO;
    fzO            = fzO + tzO;
    fzJ            = faction[j3+ZZ] - tzO;

    /* H1 block */    
    rinv2H1        = rinv1H1*rinv1H1;
    vcH1           = qH*qj*rinv1H1;
    vctot          = vctot+vcH1;
    fsH1           = vcH1*rinv2H1;
    
    txH1           = dxH1*fsH1;
    fxH1           = fxH1 + txH1;
    fxJ            = fxJ - txH1;
    tyH1           = dyH1*fsH1;
    fyH1           = fyH1 + tyH1;
    fyJ            = fyJ - tyH1;
    tzH1           = dzH1*fsH1;
    fzH1           = fzH1 + tzH1;
    fzJ            = fzJ - tzH1;
    
    /* H2 block */
    rinv2H2        = rinv1H2*rinv1H2;
    vcH2           = qH*qj*rinv1H2;
    vctot          = vctot+vcH2;
    fsH2           = vcH2*rinv2H2;
    
    txH2           = dxH2*fsH2;
    fxH2           = fxH2 + txH2;
    faction[j3+XX] = fxJ - txH2;
    tyH2           = dyH2*fsH2;
    fyH2           = fyH2 + tyH2;
    faction[j3+YY] = fyJ - tyH2;
    tzH2           = dzH2*fsH2;
    fzH2           = fzH2 + tzH2;
    faction[j3+ZZ] = fzJ - tzH2;
  }
#endif
  fw[0]=fxO;
  fw[1]=fyO;
  fw[2]=fzO;
  fw[3]=fxH1;
  fw[4]=fyH1;
  fw[5]=fzH1;
  fw[6]=fxH2;
  fw[7]=fyH2;
  fw[8]=fzH2;
  
  *Vc += vctot;
}

void c_ljcfree(real ix,real iy,real iz,int inr,
	       real pos[],int nj,t_nl_j jjnr[],
	       int  typeA[],  int typeB[],
	       real eps,
	       real chargeA[],real chargeB[],
	       real nbfpA[],  real nbfpB[],
	       real faction[],real fip[],
	       real *Vc,      real *Vnb,
	       real lambda,   real *dvdlambda,
	       real k_rf,     real c_rf,
	       real tfac,     real trunctab[])
{      
  static bool bFirst=TRUE;
  const  real twelve=12.0;
  const  real six=6.0;
  int    k,jnr,j3,tjA,tjB,ir2;
  real   fX,fY,fZ;
  real   rijX,rijY,rijZ;
  real   fijscal,rsq,vijcoul,vctot,vnbtot;
  real   rinv1,rinv2,rinv6,rinv12;
  real   fjx,fjy,fjz;
  real   tx,ty,tz,vnb6,vnb12;
  real   L1,dvdl;
  real   lam2,lam3,lam4,lam1_2,lam1_3,lam1_4;
  real   qiA,qiB,qqA,qqB,c6,c12,c6a,c6b,c12a,c12b;
  real   rffac2,qqq;
  
  fX     = 0.0;
  fY     = 0.0;
  fZ     = 0.0;
  vctot  = 0.0;
  vnbtot = 0.0;
  L1     = 1.0-lambda;
  dvdl   = 0.0;
  rffac2 = k_rf*2.0;
  
  lam2   = lambda*lambda;
  lam1_2 = L1*L1;
  lam3   = lambda*lam2;
  lam1_3 = L1*L1*L1;
  lam4   = lam2*lam2;
  lam1_4 = lam1_2*lam1_2;
  
  qiA    = eps*chargeA[inr];
  qiB    = eps*chargeB[inr];
  
  /* Cray or PGC compiler directive: ignore vector dependency */
#pragma ivdep
  for(k=0; (k<nj); k++) {
    jnr            = jjnr[k];
    j3             = 3*jnr;
    rijX           = ix - pos[j3+XX];
    rijY           = iy - pos[j3+YY];
    rijZ           = iz - pos[j3+ZZ];
    rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ);
    ir2            = tfac*rsq;
    
    rinv1          = invsqrt(rsq)*trunctab[ir2];
    rinv2          = rinv1*rinv1;
    rinv6          = rinv2*rinv2*rinv2;
  
    qqA            = qiA*chargeA[jnr];
    qqB            = qiB*chargeB[jnr];
    dvdl          += 2.0*(lambda*qqB - L1*qqA)*rinv1;
    qqq            = (lam1_2*qqA + lam2*qqB);
    vijcoul        = qqq*(rinv1 + rsq*k_rf - c_rf);
    vctot          = vctot+vijcoul;
    
    tjA            = 2*typeA[jnr];
    tjB            = 2*typeB[jnr];
    
    c6a            = nbfpA[tjA];
    c6b            = nbfpB[tjB];
    c12a           = nbfpA[tjA+1];
    c12b           = nbfpB[tjB+1];
    c6             = lam1_3*c6a  + lam3*c6b;
    c12            = lam1_4*c12a + lam4*c12b;
    vnb6           = c6  * rinv6;
    rinv12         = rinv6*rinv6;
    vnb12          = c12 * rinv12;
    dvdl          += (4.0*(lam3*c12b - lam1_3*c12a)*rinv12 - 
		      3.0*(lam2*c6b  - lam1_2*c6a)*rinv6);
    
    vnbtot         = vnbtot+vnb12-vnb6;
    fijscal        = ((twelve*vnb12-six*vnb6+qqq*rinv1)*rinv2
		      - qqq*rffac2);
    
    fjx            = faction[j3+XX];
    tx             = rijX*fijscal;
    fX             = fX + tx;
    faction[j3+XX] = fjx - tx;
    fjy            = faction[j3+YY];
    ty             = rijY*fijscal;
    fY             = fY + ty;
    faction[j3+YY] = fjy - ty;
    fjz            = faction[j3+ZZ];
    tz             = rijZ*fijscal;
    fZ             = fZ + tz;
    faction[j3+ZZ] = fjz - tz;
  }
  
  fip[XX]     = fX;
  fip[YY]     = fY;
  fip[ZZ]     = fZ;
  *Vc        += vctot;
  *Vnb       += vnbtot;
  *dvdlambda += dvdl;

  if (bFirst) {
    fprintf(stdlog,"lambda: %10g, lam2:   %10g, lam3:   %10g, lam4:   %10g\n",
	    lambda,lam2,lam3,lam4);
    fprintf(stdlog,"L1:     %10g, lam1_2: %10g, lam1_3: %10g, lam1_4: %10g\n",
	    L1,lam1_2,lam1_3,lam1_4);
    fprintf(stdlog,"qiA:    %10g, qiB:    %10g\n",
	    qiA,qiB);
    fprintf(stdlog,"qqA:    %10g, qqB:    %10g\n",
	    qqA,qqB);
    fprintf(stdlog,"eps:    %10g, nj:     %10d, dvdl:   %10g\n",eps,nj,dvdl);
    fprintf(stdlog,"fip[X]: %10g, fip[Y]: %10g, fip[Z]: %10g, Vc: %10g, Vnb: %10g\n",
	    fip[XX],fip[YY],fip[ZZ],vctot,vnbtot);
    bFirst=FALSE;
  }
}

void c_tab(real ix,real iy,real iz,real qi,
	   real pos[],int nj,int type[],int jjnr[],real charge[],
	   real nbfp[],real faction[],real fip[],
	   real *Vc,real *Vnb,int ntab,real tabscale,real VFtab[])
{
#ifdef DEBUG
  static FILE   *fp=NULL,*gp;
#endif
  int       k,jnr,j3,tj;
  real      fX,fY,fZ;
  real      rijX,rijY,rijZ;
  real      vijcoul,fijD,fijR,fijC,fijscal;
  real      fjx,fjy,fjz;
  real      tx,ty,tz,vnb6,vnb12;
  real      vctot,vnbtot;
  real      qq,c6,c12,rsq;
  real      r1,r1t,h_1;
  real      eps,eps2,Y,F,Fp,Geps,Heps2,two=2.0,VV,FF;
  int       n0,n1,nnn;
  
  fX     = 0;
  fY     = 0;
  fZ     = 0;
  vctot  = 0;
  vnbtot = 0;
  h_1    = tabscale;

#ifdef DEBUG
  if (!fp) {
    fp=ffopen("fdump.dat","w");
    gp=ffopen("ftab.dat","w");
  }
#endif
  /*
   * Inner loop using table lookup. The tables are interpolated using a cubic
   * spline alogrithm. There are separate tables for force and potential
   * but for the sake of caching performance these have been combined
   * into a single array. 
   *
   * The cubic spline interpolation looks like this:
   *
   *                         2      3           3
   * Y = (a Yi + b Yi+1) + (h /6)((a -a)Y"i + (b -b)Y"i+1)
   *
   * where a = 1-b, and Yi and Y"i are the tabulated values of a function
   * Y(x) and its second derivative. See Numerical recipes.
   *
   * The algorithm goes a little something like this:
   * 1. Calc distance vector (rijX,rijY,rijZ)
   * 2. Calc distance squared rsq
   * 3. Multiply rsq by a constant to get a table index
   * 4. Calculate fractional component (a,b) and a^3 and b^3
   * 5. Do the interpolation to calc the potential and the
   *    the scalar force
   * 6. Calc the vector force by multiplying with the distance vector
   *
   * The tables are stored as 
   * (h^2)/6  Vdisp  Vdisp" Fdisp Fdisp" Vrep Vrep" Frep Frep"
   * Vcoul Vcoul" Fcoul Fcoul"
   * In all there are 13 values in each table entry.
   * The table entries for the forces are actually divided by r so that
   * we don't have to do a sqrt calculation.
   *         
   */
  for(k=0; (k<nj); k++) {  
    jnr            = jjnr[k];
    j3             = 3*jnr;
    rijX           = ix - pos[j3];
    rijY           = iy - pos[j3+1];
    rijZ           = iz - pos[j3+2];
         
    rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ);
    r1             = sqrt(rsq);
    r1t            = r1*tabscale;
    n0             = r1t;
    n1             = 12*n0;
    eps            = (r1t-n0);
    eps2           = eps*eps;
        
#define EXTRACT(nn) { nnn = nn; Y=VFtab[nnn]; F=VFtab[nnn+1]; Geps=VFtab[nnn+2]*eps; Heps2=VFtab[nnn+3]*eps2; Fp=F+Geps+Heps2; VV=Y+eps*Fp; FF=Fp+Geps+two*Heps2; }

    /* Coulomb */
    EXTRACT(n1)
    qq             = qi*charge[jnr];
    vijcoul        = qq*VV;
    fijC           = qq*FF;
    vctot          = vctot  + vijcoul;
    
    /* Dispersion */
    tj             = 2*type[jnr];
    c6             = nbfp[tj];
    EXTRACT(n1+4)
    vnb6           = c6*VV;
    fijD           = c6*FF;
			 
    /* Repulsion */
    c12            = nbfp[tj+1];
    EXTRACT(n1+8)
    vnb12          = c12*VV;
    fijR           = c12*FF;
    vnbtot         = vnbtot + vnb12 + vnb6;
#ifdef DEBUG
    fprintf(fp,"%10g  %10g  %10g\n",r1,FF*h_1,-12.0*pow(r1,-13.0));
    
    /*    fprintf(fp,"%10g  %10g  %10g  %10g  %10g  %10d  %10d  %10g  %10g\n",r1,FF*h_1,-12.0*pow(r1,-13.0),r1,r1t,n0,n1,eps,eps2);
	  fprintf(gp,"%10g  %10g  %10g  %10g  %10g  %10g  %10g  %10g\n",r1,Y,F,Geps,Heps2,Fp,VV,FF);
	  */
#endif
    /* Total force */
    fijscal        = -(fijD + fijR + fijC)*h_1/r1;
    
    fjx            = faction[j3+XX];
    tx             = rijX*fijscal;
    fX             = fX + tx;
    faction[j3+XX] = fjx - tx;
    fjy            = faction[j3+YY];
    ty             = rijY*fijscal;
    fY             = fY + ty;
    faction[j3+YY] = fjy - ty;
    fjz            = faction[j3+ZZ];
    tz             = rijZ*fijscal;
    fZ             = fZ + tz;
    faction[j3+ZZ] = fjz - tz;
  }

  fip[XX] = fX;
  fip[YY] = fY;
  fip[ZZ] = fZ;
  *Vc   += vctot;
  *Vnb  += vnbtot;

#ifdef DEBUG  
  fflush(fp);
#endif
}

void c_bhamtab(real ix,real iy,real iz,real qi,
	       real pos[],int nj,int type[],int jjnr[],real charge[],
	       real nbfp[],real faction[],real fip[],
	       real *Vc,real *Vnb,int ntab,
	       real tabscale,real tabscale_exp,real VFtab[])
{
#ifdef DEBUG
  static FILE   *fp=NULL,*gp;
#endif
  int       k,jnr,j3,tj;
  real      fX,fY,fZ;
  real      rijX,rijY,rijZ;
  real      vijcoul,fijD,fijR,fijC,fijscal;
  real      fjx,fjy,fjz;
  real      tx,ty,tz,vnb6,vnbexp;
  real      vctot,vnbtot;
  real      qq,c6,a,b,rsq;
  real      r1,r1t,invh,invh_exp;
  real      eps,eps2,Y,F,Fp,Geps,Heps2,two=2.0,VV,FF;
  int       n0,n1,nnn;
  
  fX     = 0;
  fY     = 0;
  fZ     = 0;
  vctot  = 0;
  vnbtot = 0;
  invh     = tabscale;
  invh_exp = tabscale_exp;

#ifdef DEBUG
  if (!fp) {
    fp=ffopen("fdump.dat","w");
    gp=ffopen("ftab.dat","w");
  }
#endif
  for(k=0; (k<nj); k++) {  
    jnr            = jjnr[k];
    j3             = 3*jnr;
    rijX           = ix - pos[j3];
    rijY           = iy - pos[j3+1];
    rijZ           = iz - pos[j3+2];
         
    rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ);
    r1             = sqrt(rsq);
    r1t            = r1*tabscale;
    n0             = r1t;
    n1             = 12*n0;
    eps            = (r1t-n0);
    eps2           = eps*eps;
        
#define EXTRACT(nn) { nnn = nn; Y=VFtab[nnn]; F=VFtab[nnn+1]; Geps=VFtab[nnn+2]*eps; Heps2=VFtab[nnn+3]*eps2; Fp=F+Geps+Heps2; VV=Y+eps*Fp; FF=Fp+Geps+two*Heps2; }

    /* Coulomb */
    EXTRACT(n1)
    qq             = qi*charge[jnr];
    vijcoul        = qq*VV;
    fijC           = qq*FF;
    vctot          = vctot  + vijcoul;
    
    /* Dispersion */
    tj             = 3*type[jnr];
    c6             = nbfp[tj+2];
    EXTRACT(n1+4)
    vnb6           = c6*VV;
    fijD           = c6*FF;
			 
    /* Repulsion */
    b              = nbfp[tj+1];
    r1t            = b*r1*tabscale_exp;
    n0             = r1t;
    n1             = 12*n0;
    eps            = (r1t-n0);
    eps2           = eps*eps;
    a              = nbfp[tj];
    EXTRACT(n1+8)
    vnbexp         = a*VV;
    fijR           = a*b*FF;
    vnbtot         = vnbtot + vnbexp + vnb6;

    /* Total force */
    fijscal        = -((fijD + fijC)*invh + fijR*invh_exp)/r1;

    fjx            = faction[j3+XX];
    tx             = rijX*fijscal;
    fX             = fX + tx;
    faction[j3+XX] = fjx - tx;
    fjy            = faction[j3+YY];
    ty             = rijY*fijscal;
    fY             = fY + ty;
    faction[j3+YY] = fjy - ty;
    fjz            = faction[j3+ZZ];
    tz             = rijZ*fijscal;
    fZ             = fZ + tz;
    faction[j3+ZZ] = fjz - tz;
  }

  fip[XX] = fX;
  fip[YY] = fY;
  fip[ZZ] = fZ;
  *Vc   += vctot;
  *Vnb  += vnbtot;

#ifdef DEBUG  
  fflush(fp);
#endif
}

void c_coultab(real ix,real iy,real iz,real qi,
	       real pos[],int nj,int type[],int jjnr[],real charge[],
	       real nbfp[],real faction[],real fip[],
	       real *Vc,real *Vnb,int ntab,real tabscale,real VFtab[])
{
  int       k,jnr,j3;
  real      fX,fY,fZ;
  real      rijX,rijY,rijZ;
  real      vijcoul,fijC,fijscal;
  real      fjx,fjy,fjz;
  real      tx,ty,tz;
  real      vctot,vnbtot;
  real      qq,rsq;
  real      r1,r1t,h_1;
  real      eps,eps2,Y,F,Fp,Geps,Heps2,two=2.0,VV,FF;
  int       n0,n1,nnn;
  
  fX     = 0;
  fY     = 0;
  fZ     = 0;
  vctot  = 0;
  vnbtot = 0;
  h_1    = tabscale;
  
  /*
   * Inner loop using table lookup. The tables are interpolated using a cubic
   * spline alogrithm. There are separate tables for force and potential
   * but for the sake of caching performance these have been combined
   * into a single array. 
   *
   * The cubic spline interpolation looks like this:
   *
   *                         2      3           3
   * Y = (a Yi + b Yi+1) + (h /6)((a -a)Y"i + (b -b)Y"i+1)
   *
   * where a = 1-b, and Yi and Y"i are the tabulated values of a function
   * Y(x) and its second derivative. See Numerical recipes.
   *
   * The algorithm goes a little something like this:
   * 1. Calc distance vector (rijX,rijY,rijZ)
   * 2. Calc distance squared rsq
   * 3. Multiply rsq by a constant to get a table index
   * 4. Calculate fractional component (a,b) and a^3 and b^3
   * 5. Do the interpolation to calc the potential and the
   *    the scalar force
   * 6. Calc the vector force by multiplying with the distance vector
   *
   * The tables are stored as 
   * (h^2)/6  Vdisp  Vdisp" Fdisp Fdisp" Vrep Vrep" Frep Frep"
   * Vcoul Vcoul" Fcoul Fcoul"
   * In all there are 13 values in each table entry.
   * The table entries for the forces are actually divided by r so that
   * we don't have to do a sqrt calculation.
   *         
   */
  for(k=0; (k<nj); k++) {  
    jnr            = jjnr[k];
    j3             = 3*jnr;
    rijX           = ix - pos[j3];
    rijY           = iy - pos[j3+1];
    rijZ           = iz - pos[j3+2];
         
    rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ);
    r1             = sqrt(rsq);
    r1t            = r1*tabscale;
    n0             = r1t;
    n1             = 12*n0;
    eps            = (r1t-n0);
    eps2           = eps*eps;
        
#define EXTRACT(nn) { nnn = nn; Y=VFtab[nnn]; F=VFtab[nnn+1]; Geps=VFtab[nnn+2]*eps; Heps2=VFtab[nnn+3]*eps2; Fp=F+Geps+Heps2; VV=Y+eps*Fp; FF=Fp+Geps+two*Heps2; }

    /* Coulomb */
    EXTRACT(n1)
    qq             = qi*charge[jnr];
    vijcoul        = qq*VV;
    fijC           = qq*FF;
    vctot          = vctot  + vijcoul;
    
    /* Total force */
    fijscal        = -fijC*h_1/r1;
    
    fjx            = faction[j3+XX];
    tx             = rijX*fijscal;
    fX             = fX + tx;
    faction[j3+XX] = fjx - tx;
    fjy            = faction[j3+YY];
    ty             = rijY*fijscal;
    fY             = fY + ty;
    faction[j3+YY] = fjy - ty;
    fjz            = faction[j3+ZZ];
    tz             = rijZ*fijscal;
    fZ             = fZ + tz;
    faction[j3+ZZ] = fjz - tz;
  }

  fip[XX] = fX;
  fip[YY] = fY;
  fip[ZZ] = fZ;
  *Vc   += vctot;
  *Vnb  += vnbtot;

}

void c_free(real ix,real iy,real iz,int inr,
	    real pos[],int nj,t_nl_j jjnr[],
	    int  typeA[],  int typeB[],
	    real epsilon,
	    real chargeA[],real chargeB[],
	    real nbfpA[],  real nbfpB[],
	    real faction[],real fip[],
	    real *Vc,      real *Vnb,
	    real lambda,   real *dvdlambda,
	    int  ntab,     real tabscale,
	    real VFtab[])
{      
  static bool bFirst=TRUE;
  int    k,jnr,j3,tjA,tjB;
  real   fX,fY,fZ;
  real   rijX,rijY,rijZ;
  real   fijscal,rsq,vijcoul,vctot,vnbtot;
  real   fjx,fjy,fjz;
  real   tx,ty,tz,vnb6,vnb12;
  real   L1,dvdl;
  real   lam2,lam3,lam4,lam1_2,lam1_3,lam1_4;
  real   qiA,qiB,qqA,qqB,c6,c12,c6a,c6b,c12a,c12b;
  real   qqq,fijC,fijD,fijR;
  real   eps,eps2,Y,F,Fp,Geps,Heps2,two=2.0,VV,FF;
  real   r1,r1t,h_1;
  int    n0,n1,nnn;

  fX     = 0.0;
  fY     = 0.0;
  fZ     = 0.0;
  vctot  = 0.0;
  vnbtot = 0.0;
  L1     = 1.0-lambda;
  dvdl   = 0.0;
  
  lam2   = lambda*lambda;
  lam1_2 = L1*L1;
  lam3   = lambda*lam2;
  lam1_3 = L1*L1*L1;
  lam4   = lam2*lam2;
  lam1_4 = lam1_2*lam1_2;
  
  qiA    = epsilon*chargeA[inr];
  qiB    = epsilon*chargeB[inr];
  h_1    = tabscale;
  
  /* Cray or PGC compiler directive: ignore vector dependency */
#pragma ivdep
  for(k=0; (k<nj); k++) {
    jnr            = jjnr[k];
    j3             = 3*jnr;
    rijX           = ix - pos[j3+XX];
    rijY           = iy - pos[j3+YY];
    rijZ           = iz - pos[j3+ZZ];
    rsq            = (rijX*rijX)+(rijY*rijY)+(rijZ*rijZ);
    
    r1             = sqrt(rsq);
    r1t            = r1*tabscale;
    n0             = r1t;
    n1             = 12*n0;
    eps            = (r1t-n0);
    eps2           = eps*eps;
    
    tjA            = 2*typeA[jnr];
    tjB            = 2*typeB[jnr];
    
    qqA            = qiA*chargeA[jnr];
    qqB            = qiB*chargeB[jnr];
    qqq            = (lam1_2*qqA + lam2*qqB);
    c6a            = nbfpA[tjA];
    c6b            = nbfpB[tjB];
    c12a           = nbfpA[tjA+1];
    c12b           = nbfpB[tjB+1];
    c6             = lam1_3*c6a  + lam3*c6b;
    c12            = lam1_4*c12a + lam4*c12b;
    /* Exact shape of these functions  is not important
     * for dV/dL, therefore we continue calling it rinv6 and rinv12
     */
    EXTRACT(n1)
    vijcoul        = qqq*VV;
    fijC           = qqq*FF;
    dvdl          += 2.0*(lambda*qqB - L1*qqA)*VV;
    vctot          = vctot  + vijcoul;
    
    /* Dispersion */
    EXTRACT(n1+4)
    vnb6           = c6*VV;
    fijD           = c6*FF;
    dvdl          += 3.0*(lam2*c6b  - lam1_2*c6a)*VV;
			 
    /* Repulsion */
    EXTRACT(n1+8)
    vnb12          = c12*VV;
    fijR           = c12*FF;
    dvdl          += (4.0*(lam3*c12b - lam1_3*c12a)*VV);
    vnbtot         = vnbtot + vnb12 + vnb6;

    
    /* Total force */
    fijscal        = -(fijD + fijR + fijC)*h_1/r1;

    fjx            = faction[j3+XX];
    tx             = rijX*fijscal;
    fX             = fX + tx;
    faction[j3+XX] = fjx - tx;
    fjy            = faction[j3+YY];
    ty             = rijY*fijscal;
    fY             = fY + ty;
    faction[j3+YY] = fjy - ty;
    fjz            = faction[j3+ZZ];
    tz             = rijZ*fijscal;
    fZ             = fZ + tz;
    faction[j3+ZZ] = fjz - tz;
  }
  
  fip[XX]     = fX;
  fip[YY]     = fY;
  fip[ZZ]     = fZ;
  *Vc        += vctot;
  *Vnb       += vnbtot;
  *dvdlambda += dvdl;

  if (bFirst) {
    fprintf(stdlog,"lambda: %10g, lam2:   %10g, lam3:   %10g, lam4:   %10g\n",
	    lambda,lam2,lam3,lam4);
    fprintf(stdlog,"L1:     %10g, lam1_2: %10g, lam1_3: %10g, lam1_4: %10g\n",
	    L1,lam1_2,lam1_3,lam1_4);
    fprintf(stdlog,"qiA:    %10g, qiB:    %10g\n",
	    qiA,qiB);
    fprintf(stdlog,"qqA:    %10g, qqB:    %10g\n",
	    qqA,qqB);
    fprintf(stdlog,"eps:    %10g, nj:     %10d, dvdl:   %10g\n",eps,nj,dvdl);
    fprintf(stdlog,"fip[X]: %10g, fip[Y]: %10g, fip[Z]: %10g, Vc: %10g, Vnb: %10g\n",
	    fip[XX],fip[YY],fip[ZZ],vctot,vnbtot);
    bFirst=FALSE;
  }
}

