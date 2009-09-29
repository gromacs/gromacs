/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.5
 * Written by Christoph Junghans, Brad Lambeth, and possibly others.
 * Copyright (c) 2009 Christoph Junghans, Brad Lambeth.
 * All rights reserved.

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
#include "nb_generic_adress.h"

 void
 gmx_nb_generic_adress_kernel(t_nblist *           nlist,
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
       
     real *     wf;
     bool       bnew_wf;
     real       weight_cg1;
     real       weight_cg2;
     real       weight_product;
     real       weight_product_cg;
     real       hybscal;
     real       c_tmp;
     bool       bHybrid; /*Are we in the hybrid zone ?*/
     bool       bIntPres; /*Is interface pressure correction enabled ?*/
     bool       bCG; /*Are we calulating cg-cg interactions?*/
     wf                  = mdatoms->wf;
     bnew_wf             = fr->badress_new_wf;
     /* Check if we're doing the coarse grained charge group */
     bCG                 = (mdatoms->ptype[nlist->iinr[0]] == eptVSite);
     bIntPres            = fr->badress_pcor;
 
     /* Coarse grained - explicit interaction should be excluded via
      * energygrp_excl in the mdp file !!!!
      * due to that we know that if is bCG is true, we are calculating 
      * CG-CG interactions, else EX-EX interactions */
 
     if(bCG) 
     {
         icoul           = 0;
         ivdw            = 3;
     }
     else
     {
         icoul           = nlist->icoul;
         ivdw            = fr->adress_ivdw;
     }
     
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
         
         weight_cg1       = wf[ai0];
          
         for(k=nj0; (k<nj1); k++)
         {
             aj0              = nlist->jjnr[k];
             aj1              = nlist->jjnr_end[k];
             
            weight_cg2       = wf[aj0];
            weight_product   = weight_cg1*weight_cg2;
            weight_product_cg   = (1.0-weight_cg1)*(1.0-weight_cg2);
            /* at least one of the groups is coarse grained */
	    /* -> we should not calculate ex-ex interactions */
            if (weight_product == 0)
            {
                /* if it's a coarse grained loop, include this molecule */
                if(bCG)
                {
                    bHybrid = FALSE;
	            hybscal = 1.0;
                }
		/* if it's a explicit loop, goto the next molecule */
                else
                {
                    continue;
                }
            }
            else if ((bnew_wf && (weight_product_cg == 0)) || (!bnew_wf && (weight_product == 1)))
            {
                /* if it's a coarse grained loop, skip this molecule */
                if(bCG)
                {
                    continue;
                }
		/* if it's a explicit loop, include this molecule */
                else
                {
                    bHybrid = FALSE;
	            hybscal = 1.0;
                }
            }
            /* both have double identity, get hybrid scaling factor */
            else 
            {
                bHybrid = TRUE;
                if (bnew_wf)
                {
                    hybscal = 1.0 - exp(-0.6931472*weight_product/(weight_product_cg));
                }
                else
                {
                    /* this is the old function */
                    hybscal = weight_product;
                }
                
                /* hybscal is the scaling factor for the forces in the hybrid region
                 * for EX-EX interactions it is just given by the above formular
                 * for CG-CG interactions it is 1- the EX-EX scale */
                if(bCG)
                {
                    hybscal = 1.0 - hybscal;
                }
            }
            
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
                             
                             /* Interface pressure correction (bIntPres)
 			     * it only act in the hybrid zone (bHybrid)
 			     * it is added to the cg-cg interaction (bCG)
                              * c6 = HYB potential, c12 = CG potential
 			     * this implies that in the case of Interface 
 			     * pressure correction and CG-CG interaction
 			     * c6 = 0 otherwise */
                             if(bIntPres && bCG)
                             {
 			        if (bHybrid)
 				{
                                     c_tmp = sqrt(weight_product);
                                     if(c_tmp < 0.5)
                                     {
                                         c_tmp = cos(M_PI*c_tmp);
                                         c_tmp*= c_tmp;
                                     }
                                     else
                                     {
                                         c_tmp = (c_tmp-0.5);
                                         c_tmp*= 4.0*c_tmp;
                                     }
                                     c12 *= c_tmp;
                                     c6 = 1.0-c12;
 				}
 				else
 			        {
 				  c6=0;
 				}
                             }
                            
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
                     
                     
                    /* force weight is one anyway */
                    if (bHybrid)
                    {
                        fscal *= hybscal;
                    }
                    
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
    
     *outeriter       = nlist->nri;            
     *inneriter       = nlist->jindex[n];          	
}

