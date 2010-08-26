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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "vec.h"
#include "typedefs.h"

void
gmx_nb_free_energy_kernel(int                  icoul,
                          int                  ivdw,
                          int                  nri,
                          int *                iinr,
                          int *                jindex,
                          int *                jjnr,
                          int *                shift,
                          real *               shiftvec,
                          real *               fshift,
                          int *                gid,
                          real *               x,
                          real *               f,
                          real *               chargeA,
                          real *               chargeB,
                          real                 facel,
                          real                 krf,
                          real                 crf,
                          real                 ewc,
                          real *               Vc,
                          int *                typeA,
                          int *                typeB,
                          int                  ntype,
                          real *               nbfp,
                          real *               Vvdw,
                          real                 tabscale,
                          real *               VFtab,
                          real                 lambda,
                          real *               dvdlambda,
                          real                 alpha,
                          int                  lam_power,
                          real                 sigma6_def,
                          real                 sigma6_min,
                          gmx_bool                 bDoForces,
                          int *                outeriter,
                          int *                inneriter)
{
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    real          shX,shY,shZ;
    real          Fscal,FscalA,FscalB,tx,ty,tz;
    real          VcoulA,VcoulB,VvdwA,VvdwB;
    real          rinv6,r,rt;
    real          iqA,iqB;
    real          qqA,qqB,vcoul,vctot,krsq;
    int           ntiA,ntiB;
    int           tjA,tjB;
    real          rinvsix;
    real          Vvdw6,Vvdwtot;
    real          Vvdw12;
    real          ix,iy,iz,fix,fiy,fiz;
    real          dx,dy,dz,rsq,r4,r6,rinv;
    real          c6A,c12A,c6B,c12B;
    real          dvdl,L1,alfA,alfB,dalfA,dalfB;
    real          sigma6a,sigma6b;
    real          rA,rinvA,rinv4A,rB,rinvB,rinv4B;
    int           do_coultab,do_vdwtab,do_tab,tab_elemsize;
    int           n0,n1,nnn;
    real          Y,F,G,H,Fp,Geps,Heps2,eps,eps2,VV,FF;
    double        isp=0.564189583547756;


    /* fix compiler warnings */
    nj1 = 0;
    n1  = 0;
    eps = 0;
    eps2 = 0;
   
    dvdl = 0;
    L1   = 1.0 - lambda;

    alfA  = alpha*(lam_power==2 ? lambda*lambda : lambda);
    alfB  = alpha*(lam_power==2 ? L1*L1 : L1);
    dalfA = alpha*lam_power/6.0*(lam_power==2 ? lambda : 1); 
    dalfB = alpha*lam_power/6.0*(lam_power==2 ? L1 : 1); 

    /* Ewald table is special (icoul==5) */
    
    do_coultab = (icoul==3);
    do_vdwtab  = (ivdw==3);
    
    do_tab = do_coultab || do_vdwtab;
    
    /* we always use the combined table here */
    tab_elemsize = 12;
    
    for(n=0; (n<nri); n++)
    {
        is3              = 3*shift[n];     
        shX              = shiftvec[is3];  
        shY              = shiftvec[is3+1];
        shZ              = shiftvec[is3+2];
        nj0              = jindex[n];      
        nj1              = jindex[n+1];    
        ii               = iinr[n];        
        ii3              = 3*ii;           
        ix               = shX + x[ii3+0];
        iy               = shY + x[ii3+1];
        iz               = shZ + x[ii3+2];
        iqA              = facel*chargeA[ii];
        iqB              = facel*chargeB[ii];
        ntiA             = 2*ntype*typeA[ii];
        ntiB             = 2*ntype*typeB[ii];
        vctot            = 0;              
        Vvdwtot          = 0;              
        fix              = 0;              
        fiy              = 0;              
        fiz              = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            dx               = ix - x[j3];      
            dy               = iy - x[j3+1];      
            dz               = iz - x[j3+2];      
            rsq              = dx*dx+dy*dy+dz*dz;
            rinv             = gmx_invsqrt(rsq);
            r                = rsq*rinv;
            tjA              = ntiA+2*typeA[jnr];
            tjB              = ntiB+2*typeB[jnr];
            c6A              = nbfp[tjA];
            c6B              = nbfp[tjB];
            c12A             = nbfp[tjA+1];
            c12B             = nbfp[tjB+1];
            qqA              = iqA*chargeA[jnr]; 
            qqB              = iqB*chargeB[jnr]; 
            
            if((c6A > 0) && (c12A > 0)) 
            {
                sigma6a      = c12A/c6A;

                if (sigma6a < sigma6_min)
                {
                    sigma6a  = sigma6_min;
                }
            }
            else 
            {
                sigma6a      = sigma6_def;
            }
            if((c6B > 0) && (c12B > 0))
            {
                sigma6b      = c12B/c6B;

                if (sigma6b < sigma6_min)
                {
                    sigma6b  = sigma6_min;
                }
            }
            else
            {
                sigma6b      = sigma6_def;
            }
                        
            r4               = rsq*rsq;
            r6               = r4*rsq;
            
            FscalA           = 0;
            VcoulA           = 0;
            VvdwA            = 0;
            rinv4A           = 0;
            
            /* Only spend time on A state if it is non-zero */
            if( (qqA != 0) || (c6A != 0) || (c12A != 0) ) 
            {
                rA             = pow(alfA*sigma6a+r6,1.0/6.0);
                rinvA          = 1.0/rA;
                rinv4A         = rinvA*rinvA;
                rinv4A         = rinv4A*rinv4A;

                
                if(do_tab)
                {
                    rt         = rA*tabscale;
                    n0         = rt;
                    eps        = rt-n0;
                    eps2       = eps*eps;
                    n1         = tab_elemsize*n0;
                }
                
                if(icoul==1 || icoul==5)
                {
                    /* simple cutoff */
                    VcoulA     = qqA*rinvA;
                    FscalA     = VcoulA*rinvA*rinvA;
                }
                else if(icoul==2)
                {
                    /* reaction-field */
                    krsq       = krf*rA*rA;      
                    VcoulA     = qqA*(rinvA+krsq-crf);
                    FscalA     = qqA*(rinvA-2.0*krsq)*rinvA*rinvA;
                }
                else if(icoul==3)
                {
                    /* non-Ewald tabulated coulomb */
                    nnn        = n1;
                    Y          = VFtab[nnn];
                    F          = VFtab[nnn+1];
                    Geps       = eps*VFtab[nnn+2];
                    Heps2      = eps2*VFtab[nnn+3];
                    Fp         = F+Geps+Heps2;
                    VV         = Y+eps*Fp;
                    FF         = Fp+Geps+2.0*Heps2;
                    VcoulA     = qqA*VV;
                    FscalA     = -qqA*tabscale*FF*rinvA;                    
                }
                
                if(ivdw==1)
                {
                    /* cutoff LJ */
                    rinv6            = rinvA*rinvA*rinv4A;
                    Vvdw6            = c6A*rinv6;     
                    Vvdw12           = c12A*rinv6*rinv6;
                    VvdwA            = Vvdw12-Vvdw6;
                    FscalA          += (12.0*Vvdw12-6.0*Vvdw6)*rinvA*rinvA;                    
                }
                else if(ivdw==3)
                {
                    /* Table LJ */
		    nnn = n1+4;
                    
                    /* dispersion */
                    Y          = VFtab[nnn];
                    F          = VFtab[nnn+1];
                    Geps       = eps*VFtab[nnn+2];
                    Heps2      = eps2*VFtab[nnn+3];
                    Fp         = F+Geps+Heps2;
                    VV         = Y+eps*Fp;
                    FF         = Fp+Geps+2.0*Heps2;
                    VvdwA     += c6A*VV;
                    FscalA    -= c6A*tabscale*FF*rinvA;                    
                    
                    /* repulsion */
                    Y          = VFtab[nnn+4];
                    F          = VFtab[nnn+5];
                    Geps       = eps*VFtab[nnn+6];
                    Heps2      = eps2*VFtab[nnn+7];
                    Fp         = F+Geps+Heps2;
                    VV         = Y+eps*Fp;
                    FF         = Fp+Geps+2.0*Heps2;
                    VvdwA     += c12A*VV;
                    FscalA    -= c12A*tabscale*FF*rinvA;
                }           
                /* Buckingham vdw free energy not supported */
            }
            
            FscalB           = 0;
            VcoulB           = 0;
            VvdwB            = 0;
            rinv4B           = 0;
            
            /* Only spend time on B state if it is non-zero */
            if( (qqB != 0) || (c6B != 0) || (c12B != 0) ) 
            {
                rB             = pow(alfB*sigma6b+r6,1.0/6.0);
                rinvB          = 1.0/rB;
                rinv4B         = rinvB*rinvB;
                rinv4B         = rinv4B*rinv4B;
                
                
                if(do_tab)
                {
                    rt         = rB*tabscale;
                    n0         = rt;
                    eps        = rt-n0;
                    eps2       = eps*eps;
                    n1         = tab_elemsize*n0;
                }
                
                if(icoul==1 || icoul==5)
                {
                    /* simple cutoff */
                    VcoulB     = qqB*rinvB;
                    FscalB     = VcoulB*rinvB*rinvB;
                }
                else if(icoul==2)
                {
                    /* reaction-field */
                    krsq       = krf*rB*rB;      
                    VcoulB     = qqB*(rinvB+krsq-crf);
                    FscalB     = qqB*(rinvB-2.0*krsq)*rinvB*rinvB;                    
                }
                else if(icoul==3)
                {
                    /* non-Ewald tabulated coulomb */
                    nnn        = n1;
                    Y          = VFtab[nnn];
                    F          = VFtab[nnn+1];
                    Geps       = eps*VFtab[nnn+2];
                    Heps2      = eps2*VFtab[nnn+3];
                    Fp         = F+Geps+Heps2;
                    VV         = Y+eps*Fp;
                    FF         = Fp+Geps+2.0*Heps2;
                    VcoulB     = qqB*VV;
                    FscalB     = -qqB*tabscale*FF*rinvB;                    
                }
                
                if(ivdw==1)
                {
                    /* cutoff LJ */
                    rinv6            = rinvB*rinvB*rinv4B;
                    Vvdw6            = c6B*rinv6;     
                    Vvdw12           = c12B*rinv6*rinv6;
                    VvdwB            = Vvdw12-Vvdw6;
                    FscalB          += (12.0*Vvdw12-6.0*Vvdw6)*rinvB*rinvB;                    
                }
                else if(ivdw==3)
                {
                    /* Table LJ */
                    nnn = n1+4;
                    
                    /* dispersion */
                    Y          = VFtab[nnn];
                    F          = VFtab[nnn+1];
                    Geps       = eps*VFtab[nnn+2];
                    Heps2      = eps2*VFtab[nnn+3];
                    Fp         = F+Geps+Heps2;
                    VV         = Y+eps*Fp;
                    FF         = Fp+Geps+2.0*Heps2;
                    VvdwB     += c6B*VV;
                    FscalB    -= c6B*tabscale*FF*rinvB;                    
                    
                    /* repulsion */
                    Y          = VFtab[nnn+4];
                    F          = VFtab[nnn+5];
                    Geps       = eps*VFtab[nnn+6];
                    Heps2      = eps2*VFtab[nnn+7];
                    Fp         = F+Geps+Heps2;
                    VV         = Y+eps*Fp;
                    FF         = Fp+Geps+2.0*Heps2;
                    VvdwB     += c12B*VV;
                    FscalB    -= c12B*tabscale*FF*rinvB;                    
                }           
                /* Buckingham vdw free energy not supported */
            }

            Fscal = 0;
            
            if(icoul==5)
            {
                /* Soft-core Ewald interactions are special:
                 * For the direct space interactions we effectively want the
                 * normal coulomb interaction (added above when icoul==5),
                 * but need to subtract the part added in reciprocal space.
                 */
                if (r != 0) 
                {
                    VV    = gmx_erf(ewc*r)*rinv;
                    FF    = rinv*rinv*(VV - 2.0*ewc*isp*exp(-ewc*ewc*rsq));
                }
                else 
                {
                    VV    = ewc*2.0/sqrt(M_PI);
                    FF    = 0;
                }
                vctot  -= (lambda*qqB + L1*qqA)*VV;
                Fscal  -= (lambda*qqB + L1*qqA)*FF;
                dvdl   -= (qqB - qqA)*VV;
            }
	    
            /* Assemble A and B states */
            vctot         += lambda*VcoulB + L1*VcoulA;
            Vvdwtot       += lambda*VvdwB  + L1*VvdwA;
                
            Fscal         += (L1*FscalA*rinv4A + lambda*FscalB*rinv4B)*r4;
            dvdl          += (VcoulB + VvdwB) - (VcoulA + VvdwA);
            dvdl          += lambda*dalfB*FscalB*sigma6b*rinv4B
	                       - L1*dalfA*FscalA*sigma6a*rinv4A;
                
            if (bDoForces)
            {
                tx         = Fscal*dx;     
                ty         = Fscal*dy;     
                tz         = Fscal*dz;     
                fix        = fix + tx;      
                fiy        = fiy + ty;      
                fiz        = fiz + tz;      
                f[j3]      = f[j3]   - tx;
                f[j3+1]    = f[j3+1] - ty;
                f[j3+2]    = f[j3+2] - tz;
            }
        }
        
        if (bDoForces)
        {
            f[ii3]         = f[ii3]        + fix;
            f[ii3+1]       = f[ii3+1]      + fiy;
            f[ii3+2]       = f[ii3+2]      + fiz;
            fshift[is3]    = fshift[is3]   + fix;
            fshift[is3+1]  = fshift[is3+1] + fiy;
            fshift[is3+2]  = fshift[is3+2] + fiz;
        }
        ggid               = gid[n];
        Vc[ggid]           = Vc[ggid] + vctot;
        Vvdw[ggid]         = Vvdw[ggid] + Vvdwtot;
    }
    
    *dvdlambda      += dvdl;
    *outeriter       = nri;            
    *inneriter       = nj1;            
}
