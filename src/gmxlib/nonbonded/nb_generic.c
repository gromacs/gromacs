/*
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
#include "nb_generic.h"

void
gmx_nb_generic_kernel(t_nblist *           nlist,
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
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid,nnn,n0;
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
	
	icoul               = nlist->icoul;
	ivdw                = nlist->ivdw;

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
        ii               = nlist->iinr[n];        
        ii3              = 3*ii;           
        ix               = shX + x[ii3+0];
        iy               = shY + x[ii3+1];
        iz               = shZ + x[ii3+2];
        iq               = facel*charge[ii];
        nti              = nvdwparam*ntype*type[ii];
        vctot            = 0;              
        Vvdwtot          = 0;              
        fix              = 0;              
        fiy              = 0;              
        fiz              = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = nlist->jjnr[k];        
            j3               = 3*jnr;          
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
			
			if(icoul==3 || ivdw==3)
			{
				r                = rsq*rinv;
				rt               = r*tabscale;     
				n0               = rt;             
				eps              = rt-n0;          
				eps2             = eps*eps;        
				nnn              = table_nelements*n0;           				
			}
			
			/* Coulomb interaction. icoul==0 means no interaction */
			if(icoul>0)
			{
				qq               = iq*charge[jnr]; 

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
			if(ivdw>0)
			{
				tj               = nti+nvdwparam*type[jnr];
				
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
			
			
            tx               = fscal*dx;     
            ty               = fscal*dy;     
            tz               = fscal*dz;     
            fix              = fix + tx;      
            fiy              = fiy + ty;      
            fiz              = fiz + tz;      
            f[j3+0]          = f[j3+0] - tx;
            f[j3+1]          = f[j3+1] - ty;
            f[j3+2]          = f[j3+2] - tz;
        }
        
        f[ii3+0]         = f[ii3+0] + fix;
        f[ii3+1]         = f[ii3+1] + fiy;
        f[ii3+2]         = f[ii3+2] + fiz;
        fshift[is3]      = fshift[is3]+fix;
        fshift[is3+1]    = fshift[is3+1]+fiy;
        fshift[is3+2]    = fshift[is3+2]+fiz;
        ggid             = nlist->gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
        Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
    }
    
    *outeriter       = nlist->nri;            
    *inneriter       = nlist->jindex[n];          	
}

