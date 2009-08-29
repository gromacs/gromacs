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

#include "types/simple.h"
#include "vec.h"
#include "typedefs.h"
#include "nb_generic_cg.h"

void
gmx_nb_generic_cg_kernel(t_nblist *           nlist,
                         t_forcerec *         fr,
                         t_mdatoms *          mdatoms,
                         real *               x,
                         real *               f,
                         real *               fshift,
                         real *               Vc,
                         real *               Vvdw,
                         real                 tabscale,  
                         real *               VFtab,
                         int *                outeriter,
                         int *                inneriter)
{
    int           nri,ntype,table_nelements,icoul,ivdw;
    real          facel,gbtabscale;
    int           n,is3,i3,k,nj0,nj1,j3,ggid,nnn,n0;
    int           ai0,ai1,ai,aj0,aj1,aj;
    real          shX,shY,shZ;
    real          fscal,tx,ty,tz;
    real          rinvsq;
    real          iq;
    real          qq,vcoul,krsq,vctot;
    int           nti,nvdwparam;
    int           tj;
    real          rt,r,eps,eps2,Y,F,Geps,Heps2,VV,FF,Fp,fijD,fijR;
    real          rinvsix;
    real          Vvdwtot;
    real          Vvdw_rep,Vvdw_disp;
    real          ix,iy,iz,fix,fiy,fiz;
    real          jx,jy,jz;
    real          dx,dy,dz,rsq,rinv;
    real          c6,c12,cexp1,cexp2,br;
    real *        charge;
    real *        shiftvec;
    real *        vdwparam;
    int *         shift;
    int *         type;
    
#ifdef ADRESS
    real *     wf;
    real       weight_cg1;
    real       weight_cg2;
    real       weight_product;
    real       weight_product_cg;
    real       hybscal;
    bool       bMixed;
    bool       bIntPres;
    bool       bCG1;
    bool       bCG2;
    wf                  = mdatoms->wf;
    /* Check if we're doing the coarse grained charge group */
    bCG1                = (mdatoms->ptype[nlist->iinr[0]] == eptVSite);
    bCG2                = (mdatoms->ptype[nlist->jjnr[nlist->jindex[0]]] == eptVSite);
    bIntPres            = fr->userint2;
    /* We exclude interactions between cg and explicit atoms here
     * since the GMX_CG_INNERLOOP neighbor search doesn't exclude
     * them (yet). */
    if(bCG1 == bCG2)
    {
        if(bCG1)
        {
            icoul           = 0;
            ivdw            = 3;
        }
        else
        {
            icoul           = nlist->icoul;
            ivdw            = nlist->ivdw;
        }
#else
        icoul               = nlist->icoul;
        ivdw                = nlist->ivdw;
#endif
        /* avoid compiler warnings for cases that cannot happen */
        nnn                 = 0;
        vcoul               = 0.0;
        eps                 = 0.0;
        eps2                = 0.0;
        
        /* 3 VdW parameters for buckingham, otherwise 2 */
        nvdwparam           = (nlist->ivdw==2) ? 3 : 2;
        table_nelements     = (icoul==3) ? 4 : 0;
        table_nelements    += (ivdw==3) ? 8 : 0;
        
        charge              = mdatoms->chargeA;
        type                = mdatoms->typeA;
        facel               = fr->epsfac;
        shiftvec            = fr->shift_vec[0];
        vdwparam            = fr->nbfp;
        ntype               = fr->ntype;
        
        for(n=0; (n<nlist->nri); n++)
        {
            is3              = 3*nlist->shift[n];     
            shX              = shiftvec[is3];  
            shY              = shiftvec[is3+1];
            shZ              = shiftvec[is3+2];
            nj0              = nlist->jindex[n];      
            nj1              = nlist->jindex[n+1];    
            ai0              = nlist->iinr[n];
            ai1              = nlist->iinr_end[n];
            vctot            = 0;              
            Vvdwtot          = 0;              
            fix              = 0;
            fiy              = 0;
            fiz              = 0;
#ifdef ADRESS
            weight_cg1       = wf[ai0];
#endif
            
            for(k=nj0; (k<nj1); k++)
            {
                aj0              = nlist->jjnr[k];
                aj1              = nlist->jjnr_end[k];
                
#ifdef ADRESS
                weight_cg2       = wf[aj0];
                weight_product   = weight_cg1*weight_cg2;
                weight_product_cg   = (1.0-weight_cg1)*(1.0-weight_cg2);
                /* at least one of the groups is coarse grained */
                if (weight_product == 0)
                {
                    /* if it's a coarse grained loop, include this molecule */
                    if(bCG1)
                    {
                        bMixed = FALSE;
                    }
                    else
                    {
                        aj1 = -1;
                    }
                }
                /* at least one of the groups is explicit */
                else if (weight_product_cg == 0)
                {
                    /* if it's a coarse grained loop, skip this molecule */
                    if(bCG1)
                    {
                        aj1 = -1;
                    }
                    else
                    {
                        bMixed = FALSE;
                    }
                }
                /* both have double identity, get hybrid scaling factor */
                else 
                {
                    /* IMPORTANT: If you change the scaling function, change the return value in src/mdlib/adress.c too */
                    /* this is the old function */
                    /* hybscal = weight_product; */
                    
                    /* this is the unstretched new function */
                    /* hybscal = exp(-0.6931472*weight_product_cg*weight_product_cg/(weight_product*weight_product)); */
                    
                    /* this is the stretched new function to look like the old cos^2 function */
                    hybscal = exp(-6.2383246*weight_product_cg/(weight_product));
                    if(bCG1)
                    {
                        hybscal = 1.0 - hybscal;
                    }
                    bMixed = TRUE;
                }
#endif
                
                for(ai=ai0; (ai<ai1); ai++)
                {
                    i3               = ai*3;
                    ix               = shX + x[i3+0];
                    iy               = shY + x[i3+1];
                    iz               = shZ + x[i3+2];
                    iq               = facel*charge[ai];
                    nti              = nvdwparam*ntype*type[ai];
                    
                    /* Note that this code currently calculates
                     * all LJ and Coulomb interactions,
                     * even if the LJ parameters or charges are zero.
                     * If required, this can be optimized.
                     */
                    
                    for(aj=aj0; (aj<aj1); aj++)
                    {
                        j3               = aj*3;
                        jx               = x[j3+0];
                        jy               = x[j3+1];
                        jz               = x[j3+2];
                        dx               = ix - jx;      
                        dy               = iy - jy;      
                        dz               = iz - jz;      
                        rsq              = dx*dx+dy*dy+dz*dz;
                        rinv             = gmx_invsqrt(rsq);
                        rinvsq           = rinv*rinv;  
                        fscal            = 0;
                        
                        if (icoul==3 || ivdw==3)
                        {
                            r                = rsq*rinv;
                            rt               = r*tabscale;     
                            n0               = rt;             
                            eps              = rt-n0;          
                            eps2             = eps*eps;        
                            nnn              = table_nelements*n0;           				
                        }
                        
                        /* Coulomb interaction. icoul==0 means no interaction */
                        if (icoul > 0)
                        {
                            qq               = iq*charge[aj]; 
                            
                            switch(icoul)
                            {
                            case 1:
                                /* Vanilla cutoff coulomb */
                                vcoul            = qq*rinv;      
                                fscal            = vcoul*rinvsq; 
                                break;
                                
                            case 2:
                                /* Reaction-field */
                                krsq             = fr->k_rf*rsq;      
                                vcoul            = qq*(rinv+krsq-fr->c_rf);
                                fscal            = qq*(rinv-2.0*krsq)*rinvsq;
                                break;
                                
                            case 3:
                                /* Tabulated coulomb */
                                Y                = VFtab[nnn];     
                                F                = VFtab[nnn+1];   
                                Geps             = eps*VFtab[nnn+2];
                                Heps2            = eps2*VFtab[nnn+3];
                                nnn             += 4;
                                Fp               = F+Geps+Heps2;   
                                VV               = Y+eps*Fp;       
                                FF               = Fp+Geps+2.0*Heps2;
                                vcoul            = qq*VV;          
                                fscal            = -qq*FF*tabscale*rinv;
                                break;
                                
                            case 4:
                                /* GB */
                                gmx_fatal(FARGS,"Death & horror! GB generic interaction not implemented.\n");
                                break;
                                
                            default:
                                gmx_fatal(FARGS,"Death & horror! No generic coulomb interaction for icoul=%d.\n",icoul);
                                break;
                            }
                            vctot            = vctot+vcoul;    
                        } /* End of coulomb interactions */
                        
                        
                        /* VdW interaction. ivdw==0 means no interaction */
                        if (ivdw > 0)
                        {
                            tj               = nti+nvdwparam*type[aj];
                            
                            switch(ivdw)
                            {
                            case 1:
                                /* Vanilla Lennard-Jones cutoff */
                                c6               = vdwparam[tj];   
                                c12              = vdwparam[tj+1]; 
                                
                                rinvsix          = rinvsq*rinvsq*rinvsq;
                                Vvdw_disp        = c6*rinvsix;     
                                Vvdw_rep         = c12*rinvsix*rinvsix;
                                fscal           += (12.0*Vvdw_rep-6.0*Vvdw_disp)*rinvsq;
                                Vvdwtot          = Vvdwtot+Vvdw_rep-Vvdw_disp;
                                break;
                                
                            case 2:
                                /* Buckingham */
                                c6               = vdwparam[tj];   
                                cexp1            = vdwparam[tj+1]; 
                                cexp2            = vdwparam[tj+2]; 
                                
                                rinvsix          = rinvsq*rinvsq*rinvsq;
                                Vvdw_disp        = c6*rinvsix;     
                                br               = cexp2*rsq*rinv;
                                Vvdw_rep         = cexp1*exp(-br); 
                                fscal           += (br*Vvdw_rep-6.0*Vvdw_disp)*rinvsq;
                                Vvdwtot          = Vvdwtot+Vvdw_rep-Vvdw_disp;
                                break;
                                
                            case 3:
                                /* Tabulated VdW */
                                c6               = vdwparam[tj];   
                                c12              = vdwparam[tj+1]; 
#ifdef ADRESS
                                /* Interface pressure correction
                                 * c6 = HYB potential, c12 = CG potential */
                                if(bMixed && bIntPres && bCG1)
                                {
                                    /* could use a temporary variable here, but this should work */
                                    c12 = sqrt(weight_product);
                                    if(c12 < 0.5)
                                    {
                                        c12 = cos(M_PI*c12);
                                        c12*= c12;
                                    }
                                    else
                                    {
                                        c12 = (c12-0.25);
                                        c12*= 4.0*c12;
                                    }
                                    c6 = 1.0-c12;
                                }
#endif
                                
                                Y                = VFtab[nnn];     
                                F                = VFtab[nnn+1];   
                                Geps             = eps*VFtab[nnn+2];
                                Heps2            = eps2*VFtab[nnn+3];
                                Fp               = F+Geps+Heps2;   
                                VV               = Y+eps*Fp;       
                                FF               = Fp+Geps+2.0*Heps2;
                                Vvdw_disp        = c6*VV;          
                                fijD             = c6*FF;          
                                nnn             += 4;          
                                Y                = VFtab[nnn];     
                                F                = VFtab[nnn+1];   
                                Geps             = eps*VFtab[nnn+2];
                                Heps2            = eps2*VFtab[nnn+3];
                                Fp               = F+Geps+Heps2;   
                                VV               = Y+eps*Fp;       
                                FF               = Fp+Geps+2.0*Heps2;
                                Vvdw_rep         = c12*VV;         
                                fijR             = c12*FF;         
                                fscal           += -(fijD+fijR)*tabscale*rinv;
                                Vvdwtot          = Vvdwtot + Vvdw_disp + Vvdw_rep;						
                                break;
                                
                            default:
                                gmx_fatal(FARGS,"Death & horror! No generic VdW interaction for ivdw=%d.\n",ivdw);
                                break;
                            }
                        } /* end VdW interactions */
                        
                        
#ifdef ADRESS
                        /* force weight is one anyway */
                        if (bMixed)
                        {
                            fscal *= hybscal;
                        }
#endif
                        
                        tx               = fscal*dx;     
                        ty               = fscal*dy;     
                        tz               = fscal*dz;     
                        f[i3+0]         += tx;
                        f[i3+1]         += ty;
                        f[i3+2]         += tz;
                        f[j3+0]         -= tx;
                        f[j3+1]         -= ty;
                        f[j3+2]         -= tz;
                        fix             += tx;
                        fiy             += ty;
                        fiz             += tz;
                    }
                }
            }
            
            
            fshift[is3]     += fix;
            fshift[is3+1]   += fiy;
            fshift[is3+2]   += fiz;
            ggid             = nlist->gid[n];         
            Vc[ggid]        += vctot;
            Vvdw[ggid]      += Vvdwtot;
        }
#ifdef ADRESS
/* matches if(bCG1 == bCG2) */
    }
#endif
    
    *outeriter       = nlist->nri;            
    *inneriter       = nlist->jindex[n];          	
}

