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
                          real                 sc_r_power,
                          real                 sigma6_def,
                          real                 sigma6_min,
                          gmx_bool             bDoForces,
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
    real          dx,dy,dz,rsq,rinv;
    real          c6[NSTATES],c12[NSTATES];
    real          LFC[NSTATES],LFV[NSTATES],DLF[NSTATES];
    double        dvdl_coul,dvdl_vdw;
    real          lfac_coul[NSTATES],dlfac_coul[NSTATES],lfac_vdw[NSTATES],dlfac_vdw[NSTATES];
    real          sigma6[NSTATES],alpha_vdw_eff,alpha_coul_eff,sigma2_def,sigma2_min;
    real          rp,rpm2,rC,rV,rinvC,rpinvC,rinvV,rpinvV;
    real          sigma2[NSTATES],sigma_pow[NSTATES],sigma_powm2[NSTATES],rs,rs2;
    int           do_coultab,do_vdwtab,do_tab,tab_elemsize;
    int           n0,n1C,n1V,nnn;
    real          Y,F,G,H,Fp,Geps,Heps2,epsC,eps2C,epsV,eps2V,VV,FF;
    double        isp=0.564189583547756;
    real          dvdl_part;

    /* fix compiler warnings */
    nj1   = 0;
    n1C   = n1V   = 0;
    epsC  = epsV  = 0;
    eps2C = eps2V = 0;

    dvdl_coul  = 0;
    dvdl_vdw   = 0;

    /* Lambda factor for state A, 1-lambda*/
    LFC[STATE_A] = 1.0 - lambda_coul;
    LFV[STATE_A] = 1.0 - lambda_vdw;

    /* Lambda factor for state B, lambda*/
    LFC[STATE_B] = lambda_coul;
    LFV[STATE_B] = lambda_vdw;

    /*derivative of the lambda factor for state A and B */
    DLF[STATE_A] = -1;
    DLF[STATE_B] = 1;

    for (i=0;i<NSTATES;i++)
    {
        lfac_coul[i]  = (lam_power==2 ? (1-LFC[i])*(1-LFC[i]) : (1-LFC[i]));
        dlfac_coul[i] = DLF[i]*lam_power/sc_r_power*(lam_power==2 ? (1-LFC[i]) : 1);
        lfac_vdw[i]   = (lam_power==2 ? (1-LFV[i])*(1-LFV[i]) : (1-LFV[i]));
        dlfac_vdw[i]  = DLF[i]*lam_power/sc_r_power*(lam_power==2 ? (1-LFV[i]) : 1);
    }
    /* precalculate */
    sigma2_def = pow(sigma6_def,1.0/3.0);
    sigma2_min = pow(sigma6_min,1.0/3.0);

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
            if (sc_r_power == 6.0)
            {
                rpm2             = rsq*rsq; /* r4 */
                rp               = rpm2*rsq; /* r6 */
            }
            else if (sc_r_power == 48.0)
            {
                rp               = rsq*rsq*rsq;  /* r6 */
                rp               = rp*rp; /* r12 */
                rp               = rp*rp; /* r24 */
                rp               = rp*rp; /* r48 */
                rpm2             = rp/rsq; /* r46 */
            }
            else
            {
                rp             = pow(r,sc_r_power);  /* not currently supported as input, but can handle it */
                rpm2           = rp/rsq;
            }

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
                    sigma6[i]       = c12[i]/c6[i];
                    sigma2[i]       = pow(c12[i]/c6[i],1.0/3.0);
                    /* should be able to get rid of this ^^^ internal pow call eventually.  Will require agreement on
                       what data to store externally.  Can't be fixed without larger scale changes, so not 4.6 */
                    if (sigma6[i] < sigma6_min) { /* for disappearing coul and vdw with soft core at the same time */
                        sigma6[i] = sigma6_min;
                        sigma2[i] = sigma2_min;
                    }
                }
                else
                {
                    sigma6[i]       = sigma6_def;
                    sigma2[i]       = sigma2_def;
                }
                if (sc_r_power == 6.0)
                {
                    sigma_pow[i]    = sigma6[i];
                    sigma_powm2[i]  = sigma6[i]/sigma2[i];
                }
                else if (sc_r_power == 48.0) 
                {
                    sigma_pow[i]    = sigma6[i]*sigma6[i];   /* sigma^12 */
                    sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^24 */
                    sigma_pow[i]    = sigma_pow[i]*sigma_pow[i]; /* sigma^48 */
                    sigma_powm2[i]  = sigma_pow[i]/sigma2[i];                    
                }
                else 
                {    /* not really supported as input, but in here for testing the general case*/
                    sigma_pow[i]    = pow(sigma2[i],sc_r_power/2);
                    sigma_powm2[i]  = sigma_pow[i]/(sigma2[i]);
                }
            }

            /* only use softcore if one of the states has a zero endstate - softcore is for avoiding infinities!*/
            if ((c12[STATE_A] > 0) && (c12[STATE_B] > 0)) {
                alpha_vdw_eff    = 0;
                alpha_coul_eff   = 0;
            } else {
                alpha_vdw_eff    = alpha_vdw;
                alpha_coul_eff   = alpha_coul;
            }

            for (i=0;i<NSTATES;i++) 
            {
                FscalC[i]    = 0;
                FscalV[i]    = 0;
                Vcoul[i]     = 0;
                Vvdw[i]      = 0;

                /* Only spend time on A or B state if it is non-zero */
                if( (qq[i] != 0) || (c6[i] != 0) || (c12[i] != 0) )
                {

                    /* this section has to be inside the loop becaue of the dependence on sigma_pow */
                    rpinvC         = 1.0/(alpha_coul_eff*lfac_coul[i]*sigma_pow[i]+rp);
                    rinvC          = pow(rpinvC,1.0/sc_r_power);
                    rC             = 1.0/rinvC;
                    
                    rpinvV         = 1.0/(alpha_vdw_eff*lfac_vdw[i]*sigma_pow[i]+rp);
                    rinvV          = pow(rpinvV,1.0/sc_r_power);
                    rV             = 1.0/rinvV;

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
                        FscalC[i]  = Vcoul[i]*rpinvC;
                    }
                    else if(icoul==enbcoulRF)
                    {
                        /* reaction-field */
                        krsq       = krf*rC*rC;
                        Vcoul[i]   = qq[i]*(rinvC+krsq-crf);
                        FscalC[i]  = qq[i]*(rinvC-2.0*krsq)*rpinvC;
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
                        FscalC[i]  = -qq[i]*tabscale*FF*rC*rpinvC;
                    }

                    if(ivdw==enbvdwLJ)
                    {
                        /* cutoff LJ */
                        if (sc_r_power == 6.0)
                        {
                            rinv6            = rpinvV;
                        }
                        else
                        {
                            rinv6            = pow(rinvV,6.0);
                        }
                        Vvdw6            = c6[i]*rinv6;
                        Vvdw12           = c12[i]*rinv6*rinv6;
                        Vvdw[i]          = Vvdw12-Vvdw6;
                        FscalV[i]        = (12.0*Vvdw12-6.0*Vvdw6)*rpinvV;
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
                        FscalV[i] -= c6[i]*tabscale*FF*rV*rpinvV;

                        /* repulsion */
                        Y          = VFtab[nnn+4];
                        F          = VFtab[nnn+5];
                        Geps       = epsV*VFtab[nnn+6];
                        Heps2      = eps2V*VFtab[nnn+7];
                        Fp         = F+Geps+Heps2;
                        VV         = Y+epsV*Fp;
                        FF         = Fp+Geps+2.0*Heps2;
                        Vvdw[i]   += c12[i]*VV;
                        FscalV[i] -= c12[i]*tabscale*FF*rV*rpinvV;
                    }           
                    /* Buckingham vdw free energy not supported for now */
                }
            }

            Fscal = 0;

            if (icoul==enbcoulFEWALD) {
                /* because we compute the softcore normally,
                   we have to remove the ewald short range portion. Done outside of
                   the states loop because this part doesn't depend on the scaled R */

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

                Fscal         += LFC[i]*FscalC[i]*rpm2;
                Fscal         += LFV[i]*FscalV[i]*rpm2;

                dvdl_coul     += Vcoul[i]*DLF[i] + LFC[i]*alpha_coul_eff*dlfac_coul[i]*FscalC[i]*sigma_pow[i];
                dvdl_vdw      += Vvdw[i]*DLF[i] + LFV[i]*alpha_vdw_eff*dlfac_vdw[i]*FscalV[i]*sigma_pow[i];
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
