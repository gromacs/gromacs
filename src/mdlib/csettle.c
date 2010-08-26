/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include "vec.h"
#include "constr.h"
#include "gmx_fatal.h"
#include "smalloc.h"

typedef struct
{
    real   mO;
    real   mH;
    double wo;
    double wh;
    double wohh;
    real   dOH;
    real   dHH;
    real   ra;
    real   rb;
    real   rc;
    real   rc2;
    /* For projection */
    real   imO;
    real   imH;
    real   invdOH;
    real   invdHH;
    matrix invmat;
} settleparam_t;

typedef struct gmx_settledata
{
    settleparam_t massw;
    settleparam_t mass1;
} t_gmx_settledata;


static void init_proj_matrix(settleparam_t *p,
                             real invmO,real invmH,real dOH,real dHH)
{
    real   imOn,imHn;
    matrix mat;

    p->imO = invmO;
    p->imH = invmH;
    /* We normalize the inverse masses with imO for the matrix inversion.
     * so we can keep using masses of almost zero for frozen particles,
     * without running out of the float range in m_inv.
     */
    imOn = 1;
    imHn = p->imH/p->imO;

    /* Construct the constraint coupling matrix */
    mat[0][0] = imOn + imHn;
    mat[0][1] = imOn*(1 - 0.5*dHH*dHH/(dOH*dOH));
    mat[0][2] = imHn*0.5*dHH/dOH;
    mat[1][1] = mat[0][0];
    mat[1][2] = mat[0][2];
    mat[2][2] = imHn + imHn;
    mat[1][0] = mat[0][1];
    mat[2][0] = mat[0][2];
    mat[2][1] = mat[1][2];

    m_inv(mat,p->invmat);

    msmul(p->invmat,1/p->imO,p->invmat);

    p->invdOH = 1/dOH;
    p->invdHH = 1/dHH;
}

static void settleparam_init(settleparam_t *p,
                             real mO,real mH,real invmO,real invmH,
                             real dOH,real dHH)
{
    p->mO     = mO;
    p->mH     = mH;
    p->wo     = mO;
    p->wh     = mH;
    p->wohh   = mO + 2.0*mH;
    p->dOH    = dOH;
    p->dHH    = dHH;
    p->rc     = dHH/2.0;
    p->ra     = 2.0*p->wh*sqrt(dOH*dOH - p->rc*p->rc)/p->wohh;
    p->rb     = sqrt(dOH*dOH - p->rc*p->rc) - p->ra;
    p->rc2    = dHH;

    p->wo    /= p->wohh;
    p->wh    /= p->wohh;

    /* For projection: connection matrix inversion */
    init_proj_matrix(p,invmO,invmH,dOH,dHH);

    if (debug)
    {
        fprintf(debug,"wo = %g, wh =%g, wohh = %g, rc = %g, ra = %g\n",
                p->wo,p->wh,p->wohh,p->rc,p->ra);
        fprintf(debug,"rb = %g, rc2 = %g, dHH = %g, dOH = %g\n",
                p->rb,p->rc2,p->dHH,p->dOH);
    }
}

gmx_settledata_t settle_init(real mO,real mH,real invmO,real invmH,
                             real dOH,real dHH)
{
    gmx_settledata_t settled;

    snew(settled,1);

    settleparam_init(&settled->massw,mO,mH,invmO,invmH,dOH,dHH);

    settleparam_init(&settled->mass1,1.0,1.0,1.0,1.0,dOH,dHH);

    return settled;
}

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


void settle_proj(FILE *fp,
                 gmx_settledata_t settled,int econq,
                 int nsettle, t_iatom iatoms[],rvec x[],
                 rvec *der,rvec *derp,
                 gmx_bool bCalcVir,tensor rmdder,t_vetavars *vetavar)
{
    /* Settle for projection out constraint components
     * of derivatives of the coordinates.
     * Berk Hess 2008-1-10
     */
    
    settleparam_t *p;
    real   imO,imH,dOH,dHH,invdOH,invdHH;
    matrix invmat;
    int    i,m,m2,ow1,hw2,hw3;
    rvec   roh2,roh3,rhh,dc,fc,fcv;
    rvec   derm[3],derpm[3];
    real   invvscale,vscale_nhc,veta;
    real   kfacOH,kfacHH;

    if (econq == econqForce)
    {
        p = &settled->mass1;
    }
    else
    {
        p = &settled->massw;
    }
    imO    = p->imO;
    imH    = p->imH;
    copy_mat(p->invmat,invmat);
    dOH    = p->dOH;
    dHH    = p->dHH;
    invdOH = p->invdOH;
    invdHH = p->invdHH;
    
    veta = vetavar->veta;     
    vscale_nhc = vetavar->vscale_nhc[0]; /* assume the first temperature control group. */

#ifdef PRAGMAS
#pragma ivdep
#endif
    
    for (i=0; i<nsettle; i++)
    {
        ow1 = iatoms[i*2+1];
        hw2 = ow1 + 1;
        hw3 = ow1 + 2;


        for(m=0; m<DIM; m++)
        {
            /* in the velocity case, these are the velocities, so we 
               need to modify with the pressure control velocities! */
            
            derm[0][m]  = vscale_nhc*der[ow1][m] + veta*x[ow1][m];
            derm[1][m]  = vscale_nhc*der[hw2][m] + veta*x[hw2][m];
            derm[2][m]  = vscale_nhc*der[hw3][m] + veta*x[hw3][m];
            
        }
        /* 27 flops */
        
        for(m=0; m<DIM; m++)
        {
            roh2[m] = (x[ow1][m] - x[hw2][m])*invdOH;
            roh3[m] = (x[ow1][m] - x[hw3][m])*invdOH;
            rhh [m] = (x[hw2][m] - x[hw3][m])*invdHH;
        }
        /* 18 flops */
        
        /* Determine the projections of der(modified) on the bonds */
        clear_rvec(dc);
        for(m=0; m<DIM; m++)
        {
            dc[0] += (derm[0][m] - derm[1][m])*roh2[m];
            dc[1] += (derm[0][m] - derm[2][m])*roh3[m];
            dc[2] += (derm[1][m] - derm[2][m])*rhh [m];
        }
        /* 27 flops */
        
        /* Determine the correction for the three bonds */
        mvmul(invmat,dc,fc);
        /* 15 flops */
        
        /* divide velocity by vscale_nhc for determining constrained velocities, since they 
           have not yet been multiplied */
        svmul(1.0/vscale_nhc,fc,fcv);
        /* 7? flops */
        
        /* Subtract the corrections from derp */
        for(m=0; m<DIM; m++)
        {
            derp[ow1][m] -= imO*( fcv[0]*roh2[m] + fcv[1]*roh3[m]);
            derp[hw2][m] -= imH*(-fcv[0]*roh2[m] + fcv[2]*rhh [m]);
            derp[hw3][m] -= imH*(-fcv[1]*roh3[m] - fcv[2]*rhh [m]);
        }
        
        /* 45 flops */

        if (bCalcVir)
        {
            /* Determining r \dot m der is easy,
             * since fc contains the mass weighted corrections for der.
             */
            
            for(m=0; m<DIM; m++)
            {
                for(m2=0; m2<DIM; m2++)
                {
                    rmdder[m][m2] +=
                        dOH*roh2[m]*roh2[m2]*fcv[0] +
                        dOH*roh3[m]*roh3[m2]*fcv[1] +
                        dHH*rhh [m]*rhh [m2]*fcv[2]; 
                }
            }
        }
    }
    /* conrect rmdder, which will be used to calcualate the virial; we need to use 
       the unscaled multipliers in the virial */
    msmul(rmdder,1.0/vetavar->vscale,rmdder);
}


/* Our local shake routine to be used when settle breaks down due to a zero determinant */
static int xshake(real b4[], real after[], real dOH, real dHH, real mO, real mH) 
{  
  real bondsq[3];
  real bond[9];
  real invmass[3];
  real M2[3];
  int iconv;
  int iatom[3]={0,0,1};
  int jatom[3]={1,2,2};
  real rijx,rijy,rijz,tx,ty,tz,im,jm,acor,rp,diff;
  int i,ll,ii,jj,l3,ix,iy,iz,jx,jy,jz,conv;

  invmass[0]=1.0/mO;
  invmass[1]=1.0/mH;
  invmass[2]=1.0/mH;

  bondsq[0]=dOH*dOH;
  bondsq[1]=bondsq[0];
  bondsq[2]=dHH*dHH;
  
  M2[0]=1.0/(2.0*(invmass[0]+invmass[1]));
  M2[1]=M2[0];
  M2[2]=1.0/(2.0*(invmass[1]+invmass[2]));

  for(ll=0;ll<3;ll++) {
    l3=3*ll;
    ix=3*iatom[ll];
    jx=3*jatom[ll];
    for(i=0;i<3;i++) 
      bond[l3+i]= b4[ix+i] - b4[jx+i];
  }

  for(i=0,iconv=0;i<1000 && iconv<3; i++) {
    for(ll=0;ll<3;ll++) {
      ii = iatom[ll];
      jj = jatom[ll];
      l3 = 3*ll;
      ix = 3*ii;
      jx = 3*jj;
      iy = ix+1;
      jy = jx+1;
      iz = ix+2;
      jz = jx+2;

      rijx = bond[l3];
      rijy = bond[l3+1];
      rijz = bond[l3+2];  

      
      tx   = after[ix]-after[jx];
      ty   = after[iy]-after[jy];
      tz   = after[iz]-after[jz];
      
      rp   = tx*tx+ty*ty+tz*tz;
      diff = bondsq[ll] - rp;

      if(fabs(diff)<1e-8) {
	iconv++;
      } else {
	rp = rijx*tx+rijy*ty+rijz*tz;
	if(rp<1e-8) {
	  return -1;
	}
	acor = diff*M2[ll]/rp;
	im           = invmass[ii];
	jm           = invmass[jj];
	tx           = rijx*acor;
	ty           = rijy*acor;
	tz           = rijz*acor;
	after[ix] += tx*im;
	after[iy] += ty*im;
	after[iz] += tz*im;
	after[jx] -= tx*jm;
	after[jy] -= ty*jm;
	after[jz] -= tz*jm;
      }
    }
  }
  return 0;
}


void csettle(gmx_settledata_t settled,
             int nsettle, t_iatom iatoms[],real b4[], real after[],
             real invdt,real *v,gmx_bool bCalcVir,tensor rmdr,int *error,t_vetavars *vetavar)
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
    settleparam_t *p;
    real   mO,mH,mOs,mHs,invdts;
    /* These three weights need have double precision. Using single precision
     * can result in huge velocity and pressure deviations. */
    double wo,wh,wohh;
    real   ra,rb,rc,rc2,dOH,dHH;
    
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
    real dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
    real mdax, mday, mdaz, mdbx, mdby, mdbz, mdcx, mdcy, mdcz;
    
    int doshake;
    
    int i, shakeret, ow1, hw2, hw3;
    
    *error = -1;

    p = &settled->massw;
    mO   = p->mO;
    mH   = p->mH;
    wo   = p->wo;
    wh   = p->wh;
    wohh = p->wohh;
    rc   = p->rc;
    ra   = p->ra;
    rb   = p->rb;
    rc2  = p->rc2;
    dOH  = p->dOH;
    dHH  = p->dHH;
    
    mOs  = mO / vetavar->rvscale;
    mHs  = mH / vetavar->rvscale;
    invdts = invdt/(vetavar->rscale);
    
#ifdef PRAGMAS
#pragma ivdep
#endif
  for (i = 0; i < nsettle; ++i) {
    doshake = 0;
    /*    --- Step1  A1' ---      */
    ow1 = iatoms[i*2+1] * 3;
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
    
    axlng = gmx_invsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
    aylng = gmx_invsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
    azlng = gmx_invsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);
      
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
    tmp    = 1.0 - sinphi * sinphi;
    if (tmp <= 0) {
      *error = i;
      doshake = 1;
      cosphi = 0;
    }
    else
      cosphi = tmp*gmx_invsqrt(tmp);
    sinpsi = (zb1d - zc1d) / (rc2 * cosphi);
    tmp2   = 1.0 - sinpsi * sinpsi;
    if (tmp2 <= 0) {
      *error = i;
      doshake = 1;
      cospsi = 0;
    }
    else
      cospsi = tmp2*gmx_invsqrt(tmp2);
    /* 46 flops */
    
    if(!doshake) {
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
      tmp2   = (al2be2 - gama * gama);
      sinthe = (alpa * gama - beta * tmp2*gmx_invsqrt(tmp2)) / al2be2;
      /* 47 flops */
      
      /*  --- Step4  A3' --- */
      tmp2  = 1.0 - sinthe *sinthe;
      costhe = tmp2*gmx_invsqrt(tmp2);
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

      dax = xa3 - xa1;
      day = ya3 - ya1;
      daz = za3 - za1;
      dbx = xb3 - xb1;
      dby = yb3 - yb1;
      dbz = zb3 - zb1;
      dcx = xc3 - xc1;
      dcy = yc3 - yc1;
      dcz = zc3 - zc1;
      /* 9 flops, counted with the virial */

      if (v) {
          v[ow1]     += dax*invdts;
          v[ow1 + 1] += day*invdts;
          v[ow1 + 2] += daz*invdts;
          v[hw2]     += dbx*invdts;
          v[hw2 + 1] += dby*invdts;
          v[hw2 + 2] += dbz*invdts;
          v[hw3]     += dcx*invdts;
          v[hw3 + 1] += dcy*invdts;
          v[hw3 + 2] += dcz*invdts;
          /* 3*6 flops */
      }

      if (bCalcVir) {
          mdax = mOs*dax;
          mday = mOs*day;
          mdaz = mOs*daz;
          mdbx = mHs*dbx;
          mdby = mHs*dby;
          mdbz = mHs*dbz;
          mdcx = mHs*dcx;
          mdcy = mHs*dcy;
          mdcz = mHs*dcz;
          rmdr[XX][XX] -= b4[ow1]*mdax + b4[hw2]*mdbx + b4[hw3]*mdcx;
          rmdr[XX][YY] -= b4[ow1]*mday + b4[hw2]*mdby + b4[hw3]*mdcy;
          rmdr[XX][ZZ] -= b4[ow1]*mdaz + b4[hw2]*mdbz + b4[hw3]*mdcz;
          rmdr[YY][XX] -= b4[ow1+1]*mdax + b4[hw2+1]*mdbx + b4[hw3+1]*mdcx;
          rmdr[YY][YY] -= b4[ow1+1]*mday + b4[hw2+1]*mdby + b4[hw3+1]*mdcy;
          rmdr[YY][ZZ] -= b4[ow1+1]*mdaz + b4[hw2+1]*mdbz + b4[hw3+1]*mdcz;
          rmdr[ZZ][XX] -= b4[ow1+2]*mdax + b4[hw2+2]*mdbx + b4[hw3+2]*mdcx;
          rmdr[ZZ][YY] -= b4[ow1+2]*mday + b4[hw2+2]*mdby + b4[hw3+2]*mdcy;
          rmdr[ZZ][ZZ] -= b4[ow1+2]*mdaz + b4[hw2+2]*mdbz + b4[hw3+2]*mdcz;
	/* 3*24 - 9 flops */
      }
    } else {
      /* If we couldn't settle this water, try a simplified iterative shake instead */
        /* no pressure control in here yet */
     if(xshake(b4+ow1,after+ow1,dOH,dHH,mO,mH)!=0)
	*error=i;
    }
#ifdef DEBUG
    if (debug)
    {
        check_cons(debug,"settle",after,ow1,hw2,hw3);
    }
#endif
  }
}
