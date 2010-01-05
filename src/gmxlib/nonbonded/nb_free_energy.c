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
                          real *               Vv,
                          real                 tabscale,
                          real *               VFtab,
                          real                 lambda_coul,
                          real                 lambda_vdw,
                          real *               dvdl,
                          real                 alpha_coul,
                          real                 alpha_vdw,
                          int                  lam_power,
                          real                 def_sigma6,
                          bool                 bDoForces,
                          int *                outeriter,
                          int *                inneriter)
{
#define  STATE_A  0
#define  STATE_B  1    
#define  NSTATES  2
    int           i,j,n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    real          shX,shY,shZ;
    real          Fscal,FscalC[NSTATES],FscalV[NSTATES],tx,ty,tz;
    real          Vcoul[NSTATES],Vvdw[NSTATES];
    real          rinv6,r,rt,rtC,rtV;
    real          iqA,iqB;
    real          qq[NSTATES],vctot,krsq;
    int           ntiA,ntiB,tj[NSTATES];
    real          Vvdw6, Vvdw12,vvtot;
    real          ix,iy,iz,fix,fiy,fiz;
    real          dx,dy,dz,rsq,r4,r6,rinv;
    real          c6[NSTATES],c12[NSTATES];
    real          dvdl_vdw,dvdl_coul,LFC[NSTATES],LFV[NSTATES],DLF[NSTATES];
    real          alf_coul[NSTATES],dalf_coul[NSTATES],alf_vdw[NSTATES],dalf_vdw[NSTATES];
    real          sigma6[NSTATES],alpha_vdw_eff,alpha_coul_eff;
    real          rC,rV,rinvC,rinvV,rinv4C[NSTATES],rinv4V[NSTATES];
    int           do_coultab,do_vdwtab,do_tab,tab_elemsize;
    int           n0,n1C,n1V,nnn;
    real          Y,F,G,H,Fp,Geps,Heps2,epsC,eps2C,epsV,eps2V,VV,FF;
    double        isp=0.564189583547756;

    /* fix compiler warnings */
    nj1   = 0;
    n1C   = n1V   = 0;
    epsC  = epsV  = 0;
    eps2C = eps2V = 0;
    
    dvdl_coul  = 0;
    dvdl_vdw   = 0;

    LFC[STATE_A] = 1.0 - lambda_coul;
    LFV[STATE_A] = 1.0 - lambda_vdw;

    LFC[STATE_B] = lambda_coul;
    LFV[STATE_B] = lambda_vdw;

    DLF[STATE_A] = -1;
    DLF[STATE_B] = 1;

    /* Ewald (not PME) table is special (icoul==enbcoulFEWALD) */
    
    do_coultab = (icoul==enbcoulTAB);
    do_vdwtab  = (ivdw==enbcoulTAB);
    
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
        vvtot            = 0;              
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
            r4               = rsq*rsq;
            r6               = r4*rsq;

            tj[STATE_A]      = ntiA+2*typeA[jnr];
            tj[STATE_B]      = ntiB+2*typeB[jnr];
            qq[STATE_A]      = iqA*chargeA[jnr]; 
            qq[STATE_B]      = iqB*chargeB[jnr]; 

            for (i=0;i<NSTATES;i++) 
            {
                c6[i]              = nbfp[tj[i]];
                c12[i]             = nbfp[tj[i]+1];
                if((c6[i] > 0) && (c12[i] > 0)) 
                {
                    sigma6[i]           = c12[i]/c6[i];
                }
                else 
                {
                    sigma6[i]           = def_sigma6;
                }
            }

            /* only use softcore if one of the states has a zero endstate - it's for avoiding infinities!*/
            if((c12[STATE_A] > 0) && (c12[STATE_B] > 0)) {
                alpha_vdw_eff    = 0;
                alpha_coul_eff   = 0;
            } else {
                alpha_vdw_eff    = alpha_vdw;
                alpha_coul_eff   = alpha_coul;
            }

            j=0;
            for (i=0;i<NSTATES;i++) 
            {
                if (i==STATE_A) {
                    j=STATE_B;
                } else if (i==STATE_B) {
                    j=STATE_A;
                }
                alf_coul[i]  = alpha_coul_eff*(lam_power==2 ? LFC[j]*LFC[j] : LFC[j]);
                dalf_coul[i] = alpha_coul_eff*lam_power/6.0*(lam_power==2 ? DLF[j]*LFC[j] : DLF[j]); 
                alf_vdw[i]   = alpha_vdw_eff *(lam_power==2 ? LFV[j]*LFV[j] : LFV[j]);
                dalf_vdw[i]  = alpha_vdw_eff *lam_power/6.0*(lam_power==2 ? DLF[j]*LFV[j] : DLF[j]); 
                
                FscalC[i]    = 0;
                FscalV[i]    = 0;
                Vcoul[i]     = 0;
                Vvdw[i]      = 0;
                rinv4C[i]    = 0;
                rinv4V[i]    = 0;
            
                /* Only spend time on A or B state if it is non-zero */
                if( (qq[i] != 0) || (c6[i] != 0) || (c12[i] != 0) ) 
                {
                    rC             = pow(alf_coul[i]*sigma6[i]+r6,1.0/6.0);
                    rinvC          = 1.0/rC;
                    rinv4C[i]      = rinvC*rinvC;
                    rinv4C[i]      = rinv4C[i]*rinv4C[i];
                    
                    rV             = pow(alf_vdw[i]*sigma6[i]+r6,1.0/6.0);
                    rinvV          = 1.0/rV;
                    rinv4V[i]      = rinvV*rinvV;
                    rinv4V[i]      = rinv4V[i]*rinv4V[i];
                    
                    if (do_tab)
                    {
                        rtC        = rC*tabscale;
                        n0         = rtC;
                        epsC       = rtC-n0;
                        eps2C      = epsC*epsC;
                        n1C        = tab_elemsize*n0;
                        
                        rtV        = rV*tabscale;
                        n0         = rtV;
                        epsV       = rtV-n0;
                        eps2V      = epsV*epsV;
                        n1V        = tab_elemsize*n0;
                    }
                    
                    if(icoul==enbcoulOOR || icoul==enbcoulFEWALD)
                    {
                        /* simple cutoff */
                        Vcoul[i]   = qq[i]*rinvC;
                        FscalC[i]  = Vcoul[i]*rinvC*rinvC;
                    }
                    else if(icoul==enbcoulRF)
                    {
                        /* reaction-field */
                        krsq       = krf*rC*rC;      
                        Vcoul[i]   = qq[i]*(rinvC+krsq-crf);
                        FscalC[i]  = qq[i]*(rinvC-2.0*krsq)*rinvC*rinvC;
                    }
                    else if (icoul==enbcoulTAB)
                    {
                        /* non-Ewald tabulated coulomb */
                        nnn        = n1C;
                        Y          = VFtab[nnn];
                        F          = VFtab[nnn+1];
                        Geps       = epsC*VFtab[nnn+2];
                        Heps2      = eps2C*VFtab[nnn+3];
                        Fp         = F+Geps+Heps2;
                        VV         = Y+epsC*Fp;
                        FF         = Fp+Geps+2.0*Heps2;
                        Vcoul[i]   = qq[i]*VV;
                        FscalC[i]  = -qq[i]*tabscale*FF*rinvC;                    
                    }
                    
                    if(ivdw==enbvdwLJ)
                    {
                        /* cutoff LJ */
                        rinv6            = rinvV*rinvV*rinv4V[i];
                        Vvdw6            = c6[i]*rinv6;     
                        Vvdw12           = c12[i]*rinv6*rinv6;
                        Vvdw[i]          = Vvdw12-Vvdw6;
                        FscalV[i]        = (12.0*Vvdw12-6.0*Vvdw6)*rinvV*rinvV;                    
                    }
                    else if(ivdw==enbvdwTAB)
                    {
                        /* Table LJ */
                        nnn = n1V+4;
                        
                        /* dispersion */
                        Y          = VFtab[nnn];
                        F          = VFtab[nnn+1];
                        Geps       = epsV*VFtab[nnn+2];
                        Heps2      = eps2V*VFtab[nnn+3];
                        Fp         = F+Geps+Heps2;
                        VV         = Y+epsV*Fp;
                        FF         = Fp+Geps+2.0*Heps2;
                        Vvdw[i]   += c6[i]*VV;
                        FscalV[i] -= c6[i]*tabscale*FF*rinvV;                    
                    
                        /* repulsion */
                        Y          = VFtab[nnn+4];
                        F          = VFtab[nnn+5];
                        Geps       = epsV*VFtab[nnn+6];
                        Heps2      = eps2V*VFtab[nnn+7];
                        Fp         = F+Geps+Heps2;
                        VV         = Y+epsV*Fp;
                        FF         = Fp+Geps+2.0*Heps2;
                        Vvdw[i]   += c12[i]*VV;
                        FscalV[i] -= c12[i]*tabscale*FF*rinvV;
                    }           
                    /* Buckingham vdw free energy not supported for now */
                }
            }

            Fscal = 0;
            
            if (icoul==enbcoulFEWALD) {
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
                
                for (i=0;i<NSTATES;i++) 
                {
                    vctot      -= LFC[i]*qq[i]*VV;
                    Fscal      -= LFC[i]*qq[i]*FF;
                    dvdl_coul  -= (DLF[i]*qq[i])*VV;
                }
            }
                
            /* Assemble A and B states */
            for (i=0;i<NSTATES;i++) 
            {
                vctot         += LFC[i]*Vcoul[i];
                vvtot         += LFV[i]*Vvdw[i];
                
                Fscal         += (LFC[i]*FscalC[i]*rinv4C[i])*r4;
                Fscal         += (LFV[i]*FscalV[i]*rinv4V[i])*r4;
                
                dvdl_coul     += Vcoul[i]*DLF[i];
                dvdl_coul     += DLF[i]*LFC[i]*dalf_coul[i]*FscalC[i]*sigma6[i]*rinv4C[i];
                
                dvdl_vdw      += Vvdw[i]*DLF[i];
                dvdl_vdw      += DLF[i]*LFV[i]*dalf_vdw[i]*FscalV[i]*sigma6[i]*rinv4V[i];
            }
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
        Vv[ggid]           = Vv[ggid] + vvtot;
    }
    
    dvdl[efptCOUL]     += dvdl_coul;
    dvdl[efptVDW]      += dvdl_vdw;
    *outeriter       = nri;            
    *inneriter       = nj1;            
}
