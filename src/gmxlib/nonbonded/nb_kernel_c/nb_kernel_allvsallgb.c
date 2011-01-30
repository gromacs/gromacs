/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "types/simple.h"

#include "vec.h"
#include "smalloc.h"

#include "nb_kernel_allvsallgb.h"

typedef struct 
{
    real **    pvdwparam;
    int *      jindex;
    int **     exclusion_mask;
} 
gmx_allvsall_data_t;
		
static int 
calc_maxoffset(int i,int natoms)
{
    int maxoffset;
    
    if ((natoms % 2) == 1)
    {
        /* Odd number of atoms, easy */
        maxoffset = natoms/2;
    }
    else if ((natoms % 4) == 0)
    {
        /* Multiple of four is hard */
        if (i < natoms/2)
        {
            if ((i % 2) == 0)
            {
                maxoffset=natoms/2;
            }
            else
            {
                maxoffset=natoms/2-1;
            }
        }
        else
        {
            if ((i % 2) == 1)
            {
                maxoffset=natoms/2;
            }
            else
            {
                maxoffset=natoms/2-1;
            }
        }
    }
    else
    {
        /* natoms/2 = odd */
        if ((i % 2) == 0)
        {
            maxoffset=natoms/2;
        }
        else
        {
            maxoffset=natoms/2-1;
        }
    }
    
    return maxoffset;
}


static void
setup_exclusions_and_indices(gmx_allvsall_data_t *   aadata,
                             t_blocka *              excl, 
                             int                     natoms)
{
    int i,j,k;
    int nj0,nj1;
    int max_offset;
    int max_excl_offset;
    int iexcl;
    int nj;
    
    /* This routine can appear to be a bit complex, but it is mostly book-keeping.
     * To enable the fast all-vs-all kernel we need to be able to stream through all coordinates
     * whether they should interact or not. 
     *
     * To avoid looping over the exclusions, we create a simple mask that is 1 if the interaction
     * should be present, otherwise 0. Since exclusions typically only occur when i & j are close,
     * we create a jindex array with three elements per i atom: the starting point, the point to
     * which we need to check exclusions, and the end point.
     * This way we only have to allocate a short exclusion mask per i atom.
     */
    
    /* Allocate memory for our modified jindex array */
    snew(aadata->jindex,3*natoms);
    
    /* Pointer to lists with exclusion masks */
    snew(aadata->exclusion_mask,natoms);
    
    for(i=0;i<natoms;i++)
    {
        /* Start */
        aadata->jindex[3*i]   = i+1;
        max_offset = calc_maxoffset(i,natoms);

        /* Exclusions */
        nj0   = excl->index[i];
        nj1   = excl->index[i+1];

        /* first check the max range */
        max_excl_offset = -1;
        
        for(j=nj0; j<nj1; j++)
        {
            iexcl = excl->a[j];
            
            k = iexcl - i;
            
            if( k+natoms <= max_offset )
            {
                k+=natoms;
            }
            
            max_excl_offset = (k > max_excl_offset) ? k : max_excl_offset;
        }
        
        max_excl_offset = (max_offset < max_excl_offset) ? max_offset : max_excl_offset;
        
        aadata->jindex[3*i+1] = i+1+max_excl_offset;        
        
        snew(aadata->exclusion_mask[i],max_excl_offset);
        /* Include everything by default */
        for(j=0;j<max_excl_offset;j++)
        {
            /* Use all-ones to mark interactions that should be present, compatible with SSE */
            aadata->exclusion_mask[i][j] = 0xFFFFFFFF;
        }
        
        /* Go through exclusions again */
        for(j=nj0; j<nj1; j++)
        {
            iexcl = excl->a[j];
            
            k = iexcl - i;
            
            if( k+natoms <= max_offset )
            {
                k+=natoms;
            }
            
            if(k>0 && k<=max_excl_offset)
            {
                /* Excluded, kill it! */
                aadata->exclusion_mask[i][k-1] = 0;
            }
        }
        
        /* End */
        aadata->jindex[3*i+2] = i+1+max_offset;        
    }
}


static void
setup_aadata(gmx_allvsall_data_t **  p_aadata,
			 t_blocka *              excl, 
             int                     natoms,
             int *                   type,
             int                     ntype,
             real *                  pvdwparam)
{
	int i,j,idx;
	gmx_allvsall_data_t *aadata;
    real *p;
        
	snew(aadata,1);
	*p_aadata = aadata;
    
    /* Generate vdw params */
    snew(aadata->pvdwparam,ntype);
        
    for(i=0;i<ntype;i++)
    {
    	snew(aadata->pvdwparam[i],2*natoms);
        p = aadata->pvdwparam[i];

        /* Lets keep it simple and use multiple steps - first create temp. c6/c12 arrays */
        for(j=0;j<natoms;j++)
        {
            idx             = i*ntype+type[j];
            p[2*j]          = pvdwparam[2*idx];
            p[2*j+1]        = pvdwparam[2*idx+1];
        }        
    }
    
    setup_exclusions_and_indices(aadata,excl,natoms);
}



void
nb_kernel_allvsallgb(t_forcerec *           fr,
                     t_mdatoms *            mdatoms,
                     t_blocka *             excl,    
                     real *                 x,
                     real *                 f,
                     real *                 Vc,
                     real *                 Vvdw,
                     real *                 vpol,
                     int *                  outeriter,
                     int *                  inneriter,
                     void *                 work)
{
	gmx_allvsall_data_t *aadata;
	int        natoms;
	int        ni0,ni1;
	int        nj0,nj1,nj2;
	int        i,j,k;
	real *     charge;
	int *      type;
    real       facel;
	real *     pvdw;
	int        ggid;
    int *      mask;
    real *     GBtab;
    real       gbfactor;
    real *     invsqrta;
    real *     dvda;
    real       vgbtot,dvdasum;
    int        nnn,n0;
    
    real       ix,iy,iz,iq;
    real       fix,fiy,fiz;
    real       jx,jy,jz,qq;
    real       dx,dy,dz;
    real       tx,ty,tz;
    real       rsq,rinv,rinvsq,rinvsix;
    real       vcoul,vctot;
    real       c6,c12,Vvdw6,Vvdw12,Vvdwtot;
    real       fscal,dvdatmp,fijC,vgb;
    real       Y,F,Fp,Geps,Heps2,VV,FF,eps,eps2,r,rt;
    real       dvdaj,gbscale,isaprod,isai,isaj,gbtabscale;
    
	charge              = mdatoms->chargeA;
	type                = mdatoms->typeA;
    gbfactor            = ((1.0/fr->epsilon_r) - (1.0/fr->gb_epsilon_solvent));
	facel               = fr->epsfac;
    GBtab               = fr->gbtab.tab;
    gbtabscale          = fr->gbtab.scale;
    invsqrta            = fr->invsqrta;
    dvda                = fr->dvda;
    
    natoms              = mdatoms->nr;
	ni0                 = mdatoms->start;
	ni1                 = mdatoms->start+mdatoms->homenr;
    
    aadata = *((gmx_allvsall_data_t **)work);

	if(aadata==NULL)
	{
		setup_aadata(&aadata,excl,natoms,type,fr->ntype,fr->nbfp);
        *((gmx_allvsall_data_t **)work) = aadata;
	}

	for(i=ni0; i<ni1; i++)
	{
		/* We assume shifts are NOT used for all-vs-all interactions */
		
		/* Load i atom data */
        ix                = x[3*i];
        iy                = x[3*i+1];
        iz                = x[3*i+2];
        iq                = facel*charge[i];
        
        isai              = invsqrta[i];

        pvdw              = aadata->pvdwparam[type[i]];
        
		/* Zero the potential energy for this list */
		Vvdwtot           = 0.0;
        vctot             = 0.0;
        vgbtot            = 0.0;
        dvdasum           = 0.0;              

		/* Clear i atom forces */
        fix               = 0.0;
        fiy               = 0.0;
        fiz               = 0.0;
        
		/* Load limits for loop over neighbors */
		nj0              = aadata->jindex[3*i];
		nj1              = aadata->jindex[3*i+1];
		nj2              = aadata->jindex[3*i+2];

        mask             = aadata->exclusion_mask[i];
                
        /* Prologue part, including exclusion mask */
        for(j=nj0; j<nj1; j++,mask++)
        {          
            if(*mask!=0)
            {
                k = j%natoms;
                
                /* load j atom coordinates */
                jx                = x[3*k];
                jy                = x[3*k+1];
                jz                = x[3*k+2];
                
                /* Calculate distance */
                dx                = ix - jx;      
                dy                = iy - jy;      
                dz                = iz - jz;      
                rsq               = dx*dx+dy*dy+dz*dz;
                
                /* Calculate 1/r and 1/r2 */
                rinv             = gmx_invsqrt(rsq);
                
                /* Load parameters for j atom */
                isaj             = invsqrta[k];  
                isaprod          = isai*isaj;      
                qq               = iq*charge[k]; 
                vcoul            = qq*rinv;      
                fscal            = vcoul*rinv;   
                qq               = isaprod*(-qq)*gbfactor;  
                gbscale          = isaprod*gbtabscale;
                c6                = pvdw[2*k];
                c12               = pvdw[2*k+1];
                rinvsq           = rinv*rinv;  
                
                /* Tabulated Generalized-Born interaction */
                dvdaj            = dvda[k];      
                r                = rsq*rinv;   
                
                /* Calculate table index */
                rt               = r*gbscale;      
                n0               = rt;             
                eps              = rt-n0;          
                eps2             = eps*eps;        
                nnn              = 4*n0;           
                Y                = GBtab[nnn];     
                F                = GBtab[nnn+1];   
                Geps             = eps*GBtab[nnn+2];
                Heps2            = eps2*GBtab[nnn+3];
                Fp               = F+Geps+Heps2;   
                VV               = Y+eps*Fp;       
                FF               = Fp+Geps+2.0*Heps2;
                vgb              = qq*VV;          
                fijC             = qq*FF*gbscale;  
                dvdatmp          = -0.5*(vgb+fijC*r);
                dvdasum          = dvdasum + dvdatmp;
                dvda[k]          = dvdaj+dvdatmp*isaj*isaj;
                vctot            = vctot + vcoul;  
                vgbtot           = vgbtot + vgb;
                
                /* Lennard-Jones interaction */
                rinvsix          = rinvsq*rinvsq*rinvsq;
                Vvdw6            = c6*rinvsix;     
                Vvdw12           = c12*rinvsix*rinvsix;
                Vvdwtot          = Vvdwtot+Vvdw12-Vvdw6;
                fscal            = (12.0*Vvdw12-6.0*Vvdw6)*rinvsq-(fijC-fscal)*rinv;
                                
                /* Calculate temporary vectorial force */
                tx                = fscal*dx;     
                ty                = fscal*dy;     
                tz                = fscal*dz;     
                
                /* Increment i atom force */
                fix               = fix + tx;      
                fiy               = fiy + ty;      
                fiz               = fiz + tz;      
            
                /* Decrement j atom force */
                f[3*k]            = f[3*k]   - tx;
                f[3*k+1]          = f[3*k+1] - ty;
                f[3*k+2]          = f[3*k+2] - tz;
            }
            /* Inner loop uses 38 flops/iteration */
        }

        /* Main part, no exclusions */
        for(j=nj1; j<nj2; j++)
        {       
            k = j%natoms;

            /* load j atom coordinates */
            jx                = x[3*k];
            jy                = x[3*k+1];
            jz                = x[3*k+2];
            
            /* Calculate distance */
            dx                = ix - jx;      
            dy                = iy - jy;      
            dz                = iz - jz;      
            rsq               = dx*dx+dy*dy+dz*dz;
            
            /* Calculate 1/r and 1/r2 */
            rinv             = gmx_invsqrt(rsq);
            
            /* Load parameters for j atom */
            isaj             = invsqrta[k];  
            isaprod          = isai*isaj;      
            qq               = iq*charge[k]; 
            vcoul            = qq*rinv;      
            fscal            = vcoul*rinv;   
            qq               = isaprod*(-qq)*gbfactor;  
            gbscale          = isaprod*gbtabscale;
            c6                = pvdw[2*k];
            c12               = pvdw[2*k+1];
            rinvsq           = rinv*rinv;  
            
            /* Tabulated Generalized-Born interaction */
            dvdaj            = dvda[k];      
            r                = rsq*rinv;   
            
            /* Calculate table index */
            rt               = r*gbscale;      
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 4*n0;           
            Y                = GBtab[nnn];     
            F                = GBtab[nnn+1];   
            Geps             = eps*GBtab[nnn+2];
            Heps2            = eps2*GBtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            FF               = Fp+Geps+2.0*Heps2;
            vgb              = qq*VV;          
            fijC             = qq*FF*gbscale;  
            dvdatmp          = -0.5*(vgb+fijC*r);
            dvdasum          = dvdasum + dvdatmp;
            dvda[k]          = dvdaj+dvdatmp*isaj*isaj;
            vctot            = vctot + vcoul;  
            vgbtot           = vgbtot + vgb;

            /* Lennard-Jones interaction */
            rinvsix          = rinvsq*rinvsq*rinvsq;
            Vvdw6            = c6*rinvsix;     
            Vvdw12           = c12*rinvsix*rinvsix;
            Vvdwtot          = Vvdwtot+Vvdw12-Vvdw6;
            fscal            = (12.0*Vvdw12-6.0*Vvdw6)*rinvsq-(fijC-fscal)*rinv;
            
            /* Calculate temporary vectorial force */
            tx                = fscal*dx;     
            ty                = fscal*dy;     
            tz                = fscal*dz;     
            
            /* Increment i atom force */
            fix               = fix + tx;      
            fiy               = fiy + ty;      
            fiz               = fiz + tz;      
            
            /* Decrement j atom force */
            f[3*k]            = f[3*k]   - tx;
            f[3*k+1]          = f[3*k+1] - ty;
            f[3*k+2]          = f[3*k+2] - tz;
            
            /* Inner loop uses 38 flops/iteration */
        }
        
        f[3*i]   += fix;
        f[3*i+1] += fiy;
        f[3*i+2] += fiz;
		
		/* Add potential energies to the group for this list */
		ggid             = 0;         
        
		Vc[ggid]         = Vc[ggid] + vctot;
        Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
        vpol[ggid]       = vpol[ggid] + vgbtot;
        dvda[i]          = dvda[i] + dvdasum*isai*isai;

		/* Outer loop uses 6 flops/iteration */
	}    

    /* Write outer/inner iteration count to pointers */
    *outeriter       = ni1-ni0;         
    *inneriter       = (ni1-ni0)*natoms/2;         
}


