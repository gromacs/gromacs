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

#include "vec.h"
#include "smalloc.h"

#include "nb_kernel_allvsallgb_sse2_single.h"
#include "gmx_sse2_single.h"

#include <xmmintrin.h>
#include <emmintrin.h>

#define SIMD_WIDTH 4
#define UNROLLI    4
#define UNROLLJ    4



typedef struct 
{
	real *     x_align;
	real *     y_align;
	real *     z_align;
	real *     q_align;
	real *     fx_align;
	real *     fy_align;
	real *     fz_align;	
    real *     invsqrta_align;
    real *     dvda_align;
    real **    pvdwaram_align;
    real **    ppvdw;
    int *      jindex;
    int *      imask;
    int **     prologue_mask;
    int **     epilogue_mask;
    int        ntype;
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
setup_exclusions_and_indices_float(gmx_allvsall_data_t *   aadata,
                                   t_blocka *              excl,  
                                   int                     start,
                                   int                     end,
                                   int                     natoms)
{
    int i,j,k;
    int ni0,ni1,nj0,nj1,nj;
    int imin,imax;
    int ibase;
    int firstinteraction;
    int max_offset;
    int max_excl_offset;
    int iexcl;
    int  *pi;
    
    /* This routine can appear to be a bit complex, but it is mostly book-keeping.
     * To enable the fast all-vs-all kernel we need to be able to stream through all coordinates
     * whether they should interact or not. 
     *
     * Since we have already made an extra copy of the coordinates from natoms to 2*natoms-1,
     * atom i should interact with j-atoms i+1 <= j < i+1+maxoffset , unless they are excluded.
     * maxoffset is typically natoms/2, or natoms/2-1 when natoms is even and i>=natoms/2 
     * (the last is to make sure we only include i-j, and not also j-i interactions).
     *
     * The exclusions/inclusions are handled by a mask, set to 0xFFFFFFFF when the interaction is
     * included, and 0 otherwise.
     * 
     * Thus, we will have exclusion masks for:
     *
     * 1) j<=i
     * 2) Explicitly excluded atoms
     * 3) j>i+maxoffset
     *
     * For any normal/sane molecule, there will only be explicit exclusions for j=i+n, with n quite small (1..10).
     * By calculating this range, we can split the kernel into three parts:
     *
     * A) A prologue, where we have exclusion masks to check for all three types of exclusions
     *    (for very small systems this is the only executed part)
     * B) A main part, where no atoms are excluded
     * C) An epilogue, where we only take the last type of exclusions into account
     *
     * 
     * So, this means we need to exclusion mask lists:
     *
     * - One for the first few atoms above i (up to the max explicit exclusion range)
     * - One for the last few atoms around i+maxoffset
     *
     * 
     * Since the SIMD code needs to load aligned coordinates, the masks too will be aligned, i.e. the mask for atom
     * 5 will really start on atom 4. In addition, the kernel will be unrolled both in I and J (think 4*4 tiles),
     * so we should have groups of 4 i-atoms whose exclusion masks start at the same offset, and have the same length.
     */
    
    /* First create a simple mask to say which i atoms should be included. This is useful when the start/end positions
     * are not multiples of UNROLLI.
     */
    
    /* Example: if start=5 and end=17, we will get ni=4 and ni1=20. 
     * The i loop will this go over atoms 4, 17, 18, and 19 and addition to the ones we want to include.
     */
    ni0 = (start/UNROLLI)*UNROLLI;
    ni1 = ((end+UNROLLI-1)/UNROLLI)*UNROLLI;
    
    /* Set the interaction mask to only enable the i atoms we want to include */
    snew(pi,natoms+UNROLLI+2*SIMD_WIDTH);
    aadata->imask = (int *) (((size_t) pi + 16) & (~((size_t) 15)));
    for(i=0;i<natoms+UNROLLI;i++)
    {
        aadata->imask[i] = (i>=start && i<end) ? 0xFFFFFFFF : 0;
    }
    
    /* Allocate memory for our modified jindex array */
    snew(aadata->jindex,4*(natoms+UNROLLI));
    for(i=0;i<4*(natoms+UNROLLI);i++)
    {
        aadata->jindex[i] = 0;
    }
    
    /* Create the exclusion masks for the prologue part */
	snew(aadata->prologue_mask,natoms+UNROLLI); /* list of pointers */
	
    for(i=0;i<natoms+UNROLLI;i++)
    {
        aadata->prologue_mask[i] = NULL;
    }
    
    /* Calculate the largest exclusion range we need for each UNROLLI-tuplet of i atoms. */
    for(ibase=ni0;ibase<ni1;ibase+=UNROLLI)
	{
        max_excl_offset = -1;
        
        /* First find maxoffset for the next 4 atoms (or fewer if we are close to end) */
        imax = ((ibase+UNROLLI) < end) ? (ibase+UNROLLI) : end;
        
        /* Which atom is the first we (might) interact with? */
        imin = natoms; /* Guaranteed to be overwritten by one of 'firstinteraction' */
        for(i=ibase;i<imax;i++)
        {
            firstinteraction = i+1;
            max_offset = calc_maxoffset(i,natoms);
            
            nj0   = excl->index[i];
            nj1   = excl->index[i+1];
            for(j=nj0; j<nj1; j++)
            {
                if(excl->a[j] == firstinteraction)
                {
                    firstinteraction++;
                }
            }
            imin = (firstinteraction < imin) ? firstinteraction : imin;
        }
        /* round down to j unrolling factor */
        imin = (imin/UNROLLJ)*UNROLLJ;
        
        for(i=ibase;i<imax;i++)
        {
            max_offset = calc_maxoffset(i,natoms);
            
            nj0   = excl->index[i];
            nj1   = excl->index[i+1];
            for(j=nj0; j<nj1; j++)
            {                
                k = excl->a[j];
                
                if(k<imin)
                {
                    k += natoms;
                }
                
                if(k>i+max_offset)
                {
                    continue;
                }
                
                k = k - imin;
                
                if( k+natoms <= max_offset )
                {
                    k+=natoms;
                }
                
                max_excl_offset = (k > max_excl_offset) ? k : max_excl_offset;
            }
        }
        
        /* The offset specifies the last atom to be excluded, so add one unit to get an upper loop limit */
        max_excl_offset++;
        /* round up to j unrolling factor */
        max_excl_offset = (max_excl_offset/UNROLLJ+1)*UNROLLJ;
        
        /* Set all the prologue masks length to this value (even for i>end) */
        for(i=ibase;i<ibase+UNROLLI;i++)
        {
            aadata->jindex[4*i]   = imin;
            aadata->jindex[4*i+1] = imin+max_excl_offset;
        }        
    }
    
    /* Now the hard part, loop over it all again to calculate the actual contents of the prologue masks */
    for(ibase=ni0;ibase<ni1;ibase+=UNROLLI)
    {      
        for(i=ibase;i<ibase+UNROLLI;i++)
        {
            nj = aadata->jindex[4*i+1] - aadata->jindex[4*i];
            imin = aadata->jindex[4*i];
            
            /* Allocate aligned memory */
            snew(pi,nj+2*SIMD_WIDTH);
            aadata->prologue_mask[i] = (int *) (((size_t) pi + 16) & (~((size_t) 15)));
            
            /* If natoms is odd, maxoffset=natoms/2 
             * If natoms is even, maxoffset=natoms/2 for even atoms, natoms/2-1 for odd atoms.
             */
            max_offset = calc_maxoffset(i,natoms);
            
            /* Include interactions i+1 <= j < i+maxoffset */
            for(k=0;k<nj;k++)
            {
                j = imin + k;
                
                if( (j>i) && (j<=i+max_offset) )
                {
                    aadata->prologue_mask[i][k] = 0xFFFFFFFF;
                }
                else
                {
                    aadata->prologue_mask[i][k] = 0;
                }
            }
            
            /* Clear out the explicit exclusions */
            if(i<end)
            {
                nj0   = excl->index[i];
                nj1   = excl->index[i+1];
                for(j=nj0; j<nj1; j++)
                {                    
                    if(excl->a[j]>i+max_offset)
                    {
                        continue;
                    }
                    
                    k = excl->a[j] - i;
                    
                    if( k+natoms <= max_offset )
                    {
                        k+=natoms;
                    }
                    
                    k = k+i-imin;
                    if(k>=0)
                    {                        
                        aadata->prologue_mask[i][k] = 0;
                    }
                }
            }
        }
    }
    
    /* Construct the epilogue mask - this just contains the check for maxoffset */
    snew(aadata->epilogue_mask,natoms+UNROLLI);
    
    /* First zero everything to avoid uninitialized data */
    for(i=0;i<natoms+UNROLLI;i++)
    {
        aadata->jindex[4*i+2]    = aadata->jindex[4*i+1];
        aadata->jindex[4*i+3]    = aadata->jindex[4*i+1];
        aadata->epilogue_mask[i] = NULL;
    }
    
    for(ibase=ni0;ibase<ni1;ibase+=UNROLLI)
    {      
        /* Find the lowest index for which we need to use the epilogue */
        imin = ibase;        
        max_offset = calc_maxoffset(imin,natoms);
        
        imin = imin + 1 + max_offset;
        
        /* Find largest index for which we need to use the epilogue */
        imax = ibase + UNROLLI-1;
        imax = (imax < end) ? imax : end; 
        
        /* imax can be either odd or even */
        max_offset = calc_maxoffset(imax,natoms);
        
        imax = imax + 1 + max_offset + UNROLLJ - 1;
        
        for(i=ibase;i<ibase+UNROLLI;i++)
        {
            /* Start of epilogue - round down to j tile limit */
            aadata->jindex[4*i+2] = (imin/UNROLLJ)*UNROLLJ;
            /* Make sure we dont overlap - for small systems everything is done in the prologue */
            aadata->jindex[4*i+2] = (aadata->jindex[4*i+1] > aadata->jindex[4*i+2]) ? aadata->jindex[4*i+1] : aadata->jindex[4*i+2];
            /* Round upwards to j tile limit */
            aadata->jindex[4*i+3] = (imax/UNROLLJ)*UNROLLJ;
            /* Make sure we dont have a negative range for the epilogue */
            aadata->jindex[4*i+3] = (aadata->jindex[4*i+2] > aadata->jindex[4*i+3]) ? aadata->jindex[4*i+2] : aadata->jindex[4*i+3];
        }
    }
    
    /* And fill it with data... */
    
    for(ibase=ni0;ibase<ni1;ibase+=UNROLLI)
    {
        for(i=ibase;i<ibase+UNROLLI;i++)
        {
            
            nj = aadata->jindex[4*i+3] - aadata->jindex[4*i+2];
            
            /* Allocate aligned memory */
            snew(pi,nj+2*SIMD_WIDTH);
            aadata->epilogue_mask[i] = (int *) (((size_t) pi + 16) & (~((size_t) 15)));
            
            max_offset = calc_maxoffset(i,natoms);
            
            for(k=0;k<nj;k++)
            {
                j = aadata->jindex[4*i+2] + k;
                aadata->epilogue_mask[i][k] = (j <= i+max_offset) ? 0xFFFFFFFF : 0;
            }
        }
    }
}

static void
setup_aadata(gmx_allvsall_data_t **  p_aadata,
			 t_blocka *              excl, 
			 int                     start,
			 int                     end,
             int                     natoms,
             int *                   type,
             int                     ntype,
             real *                  pvdwaram)
{
	int i,j,k,idx;
	gmx_allvsall_data_t *aadata;
	real *pr;
    real *p;
	real *c6tmp,*c12tmp;
        
	snew(aadata,1);
	*p_aadata = aadata;
    
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->x_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->y_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->z_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->fx_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->fy_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->fz_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->q_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));	
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->invsqrta_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));	
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->dvda_align = (real *) (((size_t) pr + 16) & (~((size_t) 15)));	
    
    for(i=0;i<2*natoms+SIMD_WIDTH;i++)
	{
		aadata->x_align[i] = 0.0;
		aadata->y_align[i] = 0.0;
		aadata->z_align[i] = 0.0;
		aadata->q_align[i] = 0.0;
		aadata->fx_align[i] = 0.0;
		aadata->fy_align[i] = 0.0;
		aadata->fz_align[i] = 0.0;
		aadata->dvda_align[i] = 0.0;
	}
    
    /* Generate vdw params */
    snew(aadata->pvdwaram_align,ntype);
    snew(c6tmp,2*natoms+SIMD_WIDTH);
    snew(c12tmp,2*natoms+SIMD_WIDTH);
        
    for(i=0;i<ntype;i++)
    {
        /* Note that this 4 does NOT refer to SIMD_WIDTH, but to c6 & c12 params for 2*natoms! */
    	snew(pr,4*natoms+4*SIMD_WIDTH);
        aadata->pvdwaram_align[i] = (real *) (((size_t) pr + 16) & (~((size_t) 15)));
        p=aadata->pvdwaram_align[i];

        /* Lets keep it simple and use multiple steps - first create temp. c6/c12 arrays */
        for(j=0;j<natoms;j++)
        {
            idx             = i*ntype+type[j];
            c6tmp[j]         = pvdwaram[2*idx];
            c12tmp[j]        = pvdwaram[2*idx+1];
            c6tmp[natoms+j]  = c6tmp[j];
            c12tmp[natoms+j] = c12tmp[j];
        }
        for(j=2*natoms;j<2*natoms+SIMD_WIDTH;j++)
        {
            c6tmp[j]  = 0.0;
            c12tmp[j] = 0.0;
        }
        
        /* Merge into a single array: c6,c6,c6,c6,c12,c12,c12,c12,c6,c6,c6,c6,c12,c12,c12,c12,etc. */

         for(j=0;j<2*natoms;j+=UNROLLJ)
        {
            for(k=0;k<UNROLLJ;k++)
            {
                idx = j+k;
                p[2*j+k]         = c6tmp[idx];
                p[2*j+UNROLLJ+k] = c12tmp[idx];
            }
        }
    }
    sfree(c6tmp);
    sfree(c12tmp);
  
    snew(aadata->ppvdw,natoms+UNROLLI);
    for(i=0;i<natoms;i++)
    {
        aadata->ppvdw[i] = aadata->pvdwaram_align[type[i]];
    }
    for(i=natoms;i<natoms+UNROLLI;i++)
    {
        aadata->ppvdw[i] = aadata->pvdwaram_align[0];
    }
    
    setup_exclusions_and_indices_float(aadata,excl,start,end,natoms);
}




void
nb_kernel_allvsallgb_sse2_single(t_forcerec *           fr,
                                 t_mdatoms *            mdatoms,
                                 t_blocka *             excl,    
                                 real *                 x,
                                 real *                 f,
                                 real *                 vc,
                                 real *                 vvdw,
                                 real *                 vpol,
                                 int *                  outeriter,
                                 int *                  inneriter,
                                 void *                 work)
{
	int        natoms;
	int        ni0,ni1;
	int        nj0,nj1,nj2,nj3;
	int        i,j,k;
    real       gbfactor;
	real *     charge;
	int *      type;
    real       facel;
	real **    pvdwaram_align;
    real **    ppvdw;
    real *     pvdw0;
    real *     pvdw1;
    real *     pvdw2;
    real *     pvdw3;
    real *     GBtab;
	int        ggid;
	gmx_allvsall_data_t *aadata;
	real *     x_align;
	real *     y_align;
	real *     z_align;
	real *     fx_align;
	real *     fy_align;
	real *     fz_align;
	real *     q_align;
    real *     invsqrta_align;
    real *     dvda_align;
	int        maxoffset;
    int *      jindex;
    int **     prologue_mask;
    int **     epilogue_mask;
    int *      pmask0;
    int *      pmask1;
    int *      pmask2;
    int *      pmask3;
    int *      emask0;
    int *      emask1;
    int *      emask2;
    int *      emask3;
    int *      imask;
    int        nnn;
    real       tmpsum[4];
    
    __m128     ix_SSE0,iy_SSE0,iz_SSE0;
    __m128     ix_SSE1,iy_SSE1,iz_SSE1;
    __m128     ix_SSE2,iy_SSE2,iz_SSE2;
    __m128     ix_SSE3,iy_SSE3,iz_SSE3;
	__m128     fix_SSE0,fiy_SSE0,fiz_SSE0;
	__m128     fix_SSE1,fiy_SSE1,fiz_SSE1;
	__m128     fix_SSE2,fiy_SSE2,fiz_SSE2;
	__m128     fix_SSE3,fiy_SSE3,fiz_SSE3;
	__m128     fjx_SSE,fjy_SSE,fjz_SSE;
	__m128     jx_SSE,jy_SSE,jz_SSE,jq_SSE;
	__m128     dx_SSE0,dy_SSE0,dz_SSE0;
	__m128     dx_SSE1,dy_SSE1,dz_SSE1;
	__m128     dx_SSE2,dy_SSE2,dz_SSE2;
	__m128     dx_SSE3,dy_SSE3,dz_SSE3;
	__m128     tx_SSE0,ty_SSE0,tz_SSE0;
	__m128     tx_SSE1,ty_SSE1,tz_SSE1;
	__m128     tx_SSE2,ty_SSE2,tz_SSE2;
	__m128     tx_SSE3,ty_SSE3,tz_SSE3;
	__m128     rsq_SSE0,rinv_SSE0,rinvsq_SSE0,rinvsix_SSE0;
	__m128     rsq_SSE1,rinv_SSE1,rinvsq_SSE1,rinvsix_SSE1;
	__m128     rsq_SSE2,rinv_SSE2,rinvsq_SSE2,rinvsix_SSE2;
	__m128     rsq_SSE3,rinv_SSE3,rinvsq_SSE3,rinvsix_SSE3;
	__m128     qq_SSE0,iq_SSE0;
	__m128     qq_SSE1,iq_SSE1;
	__m128     qq_SSE2,iq_SSE2;
	__m128     qq_SSE3,iq_SSE3;
	__m128     vcoul_SSE0,vvdw6_SSE0,vvdw12_SSE0,fscal_SSE0;
	__m128     vcoul_SSE1,vvdw6_SSE1,vvdw12_SSE1,fscal_SSE1;
	__m128     vcoul_SSE2,vvdw6_SSE2,vvdw12_SSE2,fscal_SSE2;
	__m128     vcoul_SSE3,vvdw6_SSE3,vvdw12_SSE3,fscal_SSE3;
    __m128     c6_SSE0,c12_SSE0;
    __m128     c6_SSE1,c12_SSE1;
    __m128     c6_SSE2,c12_SSE2;
    __m128     c6_SSE3,c12_SSE3;
    __m128     Y_SSE0,F_SSE0,G_SSE0,H_SSE0;
    __m128     Y_SSE1,F_SSE1,G_SSE1,H_SSE1;
    __m128     Y_SSE2,F_SSE2,G_SSE2,H_SSE2;
    __m128     Y_SSE3,F_SSE3,G_SSE3,H_SSE3;
    __m128     r_SSE0,r_SSE1,r_SSE2,r_SSE3;
    __m128     rtab_SSE0,rtab_SSE1,rtab_SSE2,rtab_SSE3;
    __m128     eps_SSE0,eps_SSE1,eps_SSE2,eps_SSE3;
    __m128     vgb_SSE0,vgb_SSE1,vgb_SSE2,vgb_SSE3;
    __m128     fijGB_SSE0,fijGB_SSE1,fijGB_SSE2,fijGB_SSE3;
    __m128     isaprod_SSE0,isaprod_SSE1,isaprod_SSE2,isaprod_SSE3;
    __m128     isai_SSE,isai_SSE0,isai_SSE1,isai_SSE2,isai_SSE3;
    __m128     dvdatmp_SSE0,dvdatmp_SSE1,dvdatmp_SSE2,dvdatmp_SSE3;
    __m128     gbscale_SSE0,gbscale_SSE1,gbscale_SSE2,gbscale_SSE3;
    __m128     dvdasum_SSE0,dvdasum_SSE1,dvdasum_SSE2,dvdasum_SSE3;
    __m128     gbtabscale_SSE,gbfactor_SSE;
	__m128     vctot_SSE,vvdwtot_SSE,vgbtot_SSE,isaj_SSE;
	__m128     half_SSE,two_SSE,six_SSE,twelve_SSE;
	__m128     imask_SSE0,imask_SSE1,imask_SSE2,imask_SSE3;
    __m128     jmask_SSE0,jmask_SSE1,jmask_SSE2,jmask_SSE3;
    __m128i    n0_SSE0,n0_SSE1,n0_SSE2,n0_SSE3;
    __m128i    nnn_SSE0,nnn_SSE1,nnn_SSE2,nnn_SSE3;
    
	charge              = mdatoms->chargeA;
	type                = mdatoms->typeA;
    gbfactor            = ((1.0/fr->epsilon_r) - (1.0/fr->gb_epsilon_solvent));
	facel               = fr->epsfac;
    GBtab               = fr->gbtab.tab;
    
	natoms              = mdatoms->nr;
	ni0                 = (mdatoms->start/SIMD_WIDTH)*SIMD_WIDTH;
	ni1                 = mdatoms->start+mdatoms->homenr;
    
    half_SSE       = _mm_set1_ps(0.5);
    two_SSE        = _mm_set1_ps(2.0);
    six_SSE        = _mm_set1_ps(6.0);
	twelve_SSE     = _mm_set1_ps(12.0);
    gbfactor_SSE   = _mm_set1_ps( -gbfactor ); /* Yes, sign should be negative (saves 1 SSE instruction) */
    gbtabscale_SSE = _mm_set1_ps(fr->gbtab.scale);
    
    aadata = *((gmx_allvsall_data_t **)work);

	if(aadata==NULL)
	{
		setup_aadata(&aadata,excl,mdatoms->start,mdatoms->start+mdatoms->homenr,natoms,type,fr->ntype,fr->nbfp);
        *((gmx_allvsall_data_t **)work) = aadata;
	}
    
	x_align = aadata->x_align;
	y_align = aadata->y_align;
	z_align = aadata->z_align;
	fx_align = aadata->fx_align;
	fy_align = aadata->fy_align;
	fz_align = aadata->fz_align; 
	q_align = aadata->q_align;
	pvdwaram_align = aadata->pvdwaram_align;
    ppvdw   = aadata->ppvdw;
    dvda_align = aadata->dvda_align;
    invsqrta_align = aadata->invsqrta_align;
    prologue_mask = aadata->prologue_mask;
    epilogue_mask = aadata->epilogue_mask;
    jindex        = aadata->jindex;
    imask         = aadata->imask;
    
	for(i=ni0;i<ni1+1+natoms/2;i++)
    {
        k = i%natoms;
		x_align[i]       = x[3*k];
		y_align[i]       = x[3*k+1];
		z_align[i]       = x[3*k+2];
		q_align[i]       = charge[k];
        invsqrta_align[i] = fr->invsqrta[k];
		fx_align[i]       = 0;
		fy_align[i]       = 0;
		fz_align[i]       = 0;
        dvda_align[i]     = 0;
	}
    
    
	for(i=ni0; i<ni1; i+=UNROLLI)
	{
		/* We assume shifts are NOT used for all-vs-all interactions */
		
		/* Load i atom data */
		ix_SSE0          = _mm_load1_ps(x_align+i);
		ix_SSE1          = _mm_load1_ps(x_align+i+1);
		ix_SSE2          = _mm_load1_ps(x_align+i+2);
		ix_SSE3          = _mm_load1_ps(x_align+i+3);
		iy_SSE0          = _mm_load1_ps(y_align+i);
		iy_SSE1          = _mm_load1_ps(y_align+i+1);
		iy_SSE2          = _mm_load1_ps(y_align+i+2);
		iy_SSE3          = _mm_load1_ps(y_align+i+3);
		iz_SSE0          = _mm_load1_ps(z_align+i);
		iz_SSE1          = _mm_load1_ps(z_align+i+1);
		iz_SSE2          = _mm_load1_ps(z_align+i+2);
		iz_SSE3          = _mm_load1_ps(z_align+i+3);
		iq_SSE0          = _mm_set1_ps(facel*q_align[i]);
		iq_SSE1          = _mm_set1_ps(facel*q_align[i+1]);
		iq_SSE2          = _mm_set1_ps(facel*q_align[i+2]);
		iq_SSE3          = _mm_set1_ps(facel*q_align[i+3]);
        
        isai_SSE         = _mm_load_ps(invsqrta_align+i);
        isai_SSE0        = _mm_shuffle_ps(isai_SSE,isai_SSE,_MM_SHUFFLE(0,0,0,0));
        isai_SSE1        = _mm_shuffle_ps(isai_SSE,isai_SSE,_MM_SHUFFLE(1,1,1,1));
        isai_SSE2        = _mm_shuffle_ps(isai_SSE,isai_SSE,_MM_SHUFFLE(2,2,2,2));
        isai_SSE3        = _mm_shuffle_ps(isai_SSE,isai_SSE,_MM_SHUFFLE(3,3,3,3));
        
        pvdw0            = ppvdw[i];
        pvdw1            = ppvdw[i+1];
        pvdw2            = ppvdw[i+2];
        pvdw3            = ppvdw[i+3];

		/* Zero the potential energy for this list */
		vvdwtot_SSE      = _mm_setzero_ps();
		vctot_SSE        = _mm_setzero_ps();
        vgbtot_SSE       = _mm_setzero_ps();
        
		/* Clear i atom forces */
		fix_SSE0           = _mm_setzero_ps();
		fix_SSE1           = _mm_setzero_ps();
		fix_SSE2           = _mm_setzero_ps();
		fix_SSE3           = _mm_setzero_ps();
		fiy_SSE0           = _mm_setzero_ps();
		fiy_SSE1           = _mm_setzero_ps();
		fiy_SSE2           = _mm_setzero_ps();
		fiy_SSE3           = _mm_setzero_ps();
		fiz_SSE0           = _mm_setzero_ps();
		fiz_SSE1           = _mm_setzero_ps();
		fiz_SSE2           = _mm_setzero_ps();
		fiz_SSE3           = _mm_setzero_ps();
        dvdasum_SSE0       = _mm_setzero_ps();
        dvdasum_SSE1       = _mm_setzero_ps();
        dvdasum_SSE2       = _mm_setzero_ps();
        dvdasum_SSE3       = _mm_setzero_ps();
        
		/* Load limits for loop over neighbors */
		nj0              = jindex[4*i];
		nj1              = jindex[4*i+1];
		nj2              = jindex[4*i+2];
		nj3              = jindex[4*i+3];

        pmask0           = prologue_mask[i];
        pmask1           = prologue_mask[i+1];
        pmask2           = prologue_mask[i+2];
        pmask3           = prologue_mask[i+3];
        emask0           = epilogue_mask[i];
        emask1           = epilogue_mask[i+1];
        emask2           = epilogue_mask[i+2];
        emask3           = epilogue_mask[i+3];
        imask_SSE0        = _mm_load1_ps((real *)(imask+i));
        imask_SSE1        = _mm_load1_ps((real *)(imask+i+1));
        imask_SSE2        = _mm_load1_ps((real *)(imask+i+2));
        imask_SSE3        = _mm_load1_ps((real *)(imask+i+3));
         
        for(j=nj0; j<nj1; j+=UNROLLJ)
        {            
            jmask_SSE0 = _mm_load_ps((real *)pmask0);
            jmask_SSE1 = _mm_load_ps((real *)pmask1);
            jmask_SSE2 = _mm_load_ps((real *)pmask2);
            jmask_SSE3 = _mm_load_ps((real *)pmask3);
            pmask0 += UNROLLJ;
            pmask1 += UNROLLJ;
            pmask2 += UNROLLJ;
            pmask3 += UNROLLJ;
         
            /* load j atom coordinates */
            jx_SSE           = _mm_load_ps(x_align+j);
            jy_SSE           = _mm_load_ps(y_align+j);
            jz_SSE           = _mm_load_ps(z_align+j);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0,jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0,jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0,jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1,jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1,jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1,jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2,jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2,jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2,jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3,jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3,jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3,jz_SSE);
            
            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0,dy_SSE0,dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1,dy_SSE1,dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2,dy_SSE2,dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3,dy_SSE3,dz_SSE3);

            /* Combine masks */
            jmask_SSE0         = _mm_and_ps(jmask_SSE0,imask_SSE0);
            jmask_SSE1         = _mm_and_ps(jmask_SSE1,imask_SSE1);
            jmask_SSE2         = _mm_and_ps(jmask_SSE2,imask_SSE2);
            jmask_SSE3         = _mm_and_ps(jmask_SSE3,imask_SSE3);
            
            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);
            
            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0,jmask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1,jmask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2,jmask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3,jmask_SSE3);
                        
            /* Load parameters for j atom */
            jq_SSE             = _mm_load_ps(q_align+j);
            qq_SSE0            = _mm_mul_ps(iq_SSE0,jq_SSE);
            qq_SSE1            = _mm_mul_ps(iq_SSE1,jq_SSE);
            qq_SSE2            = _mm_mul_ps(iq_SSE2,jq_SSE);
            qq_SSE3            = _mm_mul_ps(iq_SSE3,jq_SSE);
            qq_SSE0            = _mm_and_ps(qq_SSE0,jmask_SSE0);
            qq_SSE1            = _mm_and_ps(qq_SSE1,jmask_SSE1);
            qq_SSE2            = _mm_and_ps(qq_SSE2,jmask_SSE2);
            qq_SSE3            = _mm_and_ps(qq_SSE3,jmask_SSE3);
            
            c6_SSE0            = _mm_load_ps(pvdw0+2*j);
            c6_SSE1            = _mm_load_ps(pvdw1+2*j);
            c6_SSE2            = _mm_load_ps(pvdw2+2*j);
            c6_SSE3            = _mm_load_ps(pvdw3+2*j);
            c12_SSE0           = _mm_load_ps(pvdw0+2*j+UNROLLJ);
            c12_SSE1           = _mm_load_ps(pvdw1+2*j+UNROLLJ);
            c12_SSE2           = _mm_load_ps(pvdw2+2*j+UNROLLJ);
            c12_SSE3           = _mm_load_ps(pvdw3+2*j+UNROLLJ);
            
            isaj_SSE           = _mm_load_ps(invsqrta_align+j);
            isaprod_SSE0       = _mm_mul_ps(isai_SSE0,isaj_SSE);
            isaprod_SSE1       = _mm_mul_ps(isai_SSE1,isaj_SSE);
            isaprod_SSE2       = _mm_mul_ps(isai_SSE2,isaj_SSE);
            isaprod_SSE3       = _mm_mul_ps(isai_SSE3,isaj_SSE);
            
            rinvsq_SSE0        = _mm_mul_ps(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1        = _mm_mul_ps(rinv_SSE1,rinv_SSE1);
            rinvsq_SSE2        = _mm_mul_ps(rinv_SSE2,rinv_SSE2);
            rinvsq_SSE3        = _mm_mul_ps(rinv_SSE3,rinv_SSE3);
            
            /* Coulomb interaction */
            vcoul_SSE0         = _mm_mul_ps(qq_SSE0,rinv_SSE0);
            vcoul_SSE1         = _mm_mul_ps(qq_SSE1,rinv_SSE1);
            vcoul_SSE2         = _mm_mul_ps(qq_SSE2,rinv_SSE2);
            vcoul_SSE3         = _mm_mul_ps(qq_SSE3,rinv_SSE3);

            fscal_SSE0         = _mm_mul_ps(vcoul_SSE0,rinv_SSE0);
            fscal_SSE1         = _mm_mul_ps(vcoul_SSE1,rinv_SSE1);
            fscal_SSE2         = _mm_mul_ps(vcoul_SSE2,rinv_SSE2);
            fscal_SSE3         = _mm_mul_ps(vcoul_SSE3,rinv_SSE3);
                        
            vctot_SSE          = _mm_add_ps(vctot_SSE, gmx_mm_sum4_ps(vcoul_SSE0,vcoul_SSE1,vcoul_SSE2,vcoul_SSE3));
            
            /* Polarization interaction */
            qq_SSE0            = _mm_mul_ps(qq_SSE0,isaprod_SSE0);
            qq_SSE1            = _mm_mul_ps(qq_SSE1,isaprod_SSE1);
            qq_SSE2            = _mm_mul_ps(qq_SSE2,isaprod_SSE2);
            qq_SSE3            = _mm_mul_ps(qq_SSE3,isaprod_SSE3);
            qq_SSE0            = _mm_mul_ps(qq_SSE0,gbfactor_SSE);
            qq_SSE1            = _mm_mul_ps(qq_SSE1,gbfactor_SSE);
            qq_SSE2            = _mm_mul_ps(qq_SSE2,gbfactor_SSE);
            qq_SSE3            = _mm_mul_ps(qq_SSE3,gbfactor_SSE);
            
            gbscale_SSE0       = _mm_mul_ps(isaprod_SSE0,gbtabscale_SSE);
            gbscale_SSE1       = _mm_mul_ps(isaprod_SSE1,gbtabscale_SSE);
            gbscale_SSE2       = _mm_mul_ps(isaprod_SSE2,gbtabscale_SSE);
            gbscale_SSE3       = _mm_mul_ps(isaprod_SSE3,gbtabscale_SSE);
            
            /* Use tables for (rescaled) GB potential */
            r_SSE0             = _mm_mul_ps(rsq_SSE0,rinv_SSE0);
            r_SSE1             = _mm_mul_ps(rsq_SSE1,rinv_SSE1);
            r_SSE2             = _mm_mul_ps(rsq_SSE2,rinv_SSE2);
            r_SSE3             = _mm_mul_ps(rsq_SSE3,rinv_SSE3);
            rtab_SSE0          = _mm_mul_ps(r_SSE0,gbscale_SSE0);
            rtab_SSE1          = _mm_mul_ps(r_SSE1,gbscale_SSE1);
            rtab_SSE2          = _mm_mul_ps(r_SSE2,gbscale_SSE2);
            rtab_SSE3          = _mm_mul_ps(r_SSE3,gbscale_SSE3);
            
            n0_SSE0            = _mm_cvttps_epi32(rtab_SSE0);
            n0_SSE1            = _mm_cvttps_epi32(rtab_SSE1);
            n0_SSE2            = _mm_cvttps_epi32(rtab_SSE2);
            n0_SSE3            = _mm_cvttps_epi32(rtab_SSE3);
            eps_SSE0           = _mm_sub_ps(rtab_SSE0 , _mm_cvtepi32_ps(n0_SSE0) );
            eps_SSE1           = _mm_sub_ps(rtab_SSE1 , _mm_cvtepi32_ps(n0_SSE1) );
            eps_SSE2           = _mm_sub_ps(rtab_SSE2 , _mm_cvtepi32_ps(n0_SSE2) );
            eps_SSE3           = _mm_sub_ps(rtab_SSE3 , _mm_cvtepi32_ps(n0_SSE3) );
            nnn_SSE0           = _mm_slli_epi32(n0_SSE0,2);
            nnn_SSE1           = _mm_slli_epi32(n0_SSE1,2);
            nnn_SSE2           = _mm_slli_epi32(n0_SSE2,2);
            nnn_SSE3           = _mm_slli_epi32(n0_SSE3,2);
            
            /* the tables are 16-byte aligned, so we can use _mm_load_ps */			
            Y_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,0)); /* YFGH */
            F_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,1)); /* YFGH */
            G_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,2)); /* YFGH */
            H_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE0,F_SSE0,G_SSE0,H_SSE0);
            
            Y_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,0)); /* YFGH */
            F_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,1)); /* YFGH */
            G_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,2)); /* YFGH */
            H_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE1,F_SSE1,G_SSE1,H_SSE1);
            
            Y_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,0)); /* YFGH */
            F_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,1)); /* YFGH */
            G_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,2)); /* YFGH */
            H_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE2,F_SSE2,G_SSE2,H_SSE2);
            
            Y_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,0)); /* YFGH */
            F_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,1)); /* YFGH */
            G_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,2)); /* YFGH */
            H_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE3,F_SSE3,G_SSE3,H_SSE3);
                        
            /* G*eps */
            G_SSE0             = _mm_mul_ps(G_SSE0,eps_SSE0);
            G_SSE1             = _mm_mul_ps(G_SSE1,eps_SSE1);
            G_SSE2             = _mm_mul_ps(G_SSE2,eps_SSE2);
            G_SSE3             = _mm_mul_ps(G_SSE3,eps_SSE3);
            
            /* H*eps2 */
            H_SSE0             = _mm_mul_ps(H_SSE0, _mm_mul_ps(eps_SSE0,eps_SSE0) );
            H_SSE1             = _mm_mul_ps(H_SSE1, _mm_mul_ps(eps_SSE1,eps_SSE1) );
            H_SSE2             = _mm_mul_ps(H_SSE2, _mm_mul_ps(eps_SSE2,eps_SSE2) );
            H_SSE3             = _mm_mul_ps(H_SSE3, _mm_mul_ps(eps_SSE3,eps_SSE3) );

            /* Fp=F+G*eps+H*eps2 */
            F_SSE0             = _mm_add_ps(F_SSE0, _mm_add_ps( G_SSE0 , H_SSE0 ) );
            F_SSE1             = _mm_add_ps(F_SSE1, _mm_add_ps( G_SSE1 , H_SSE1 ) );
            F_SSE2             = _mm_add_ps(F_SSE2, _mm_add_ps( G_SSE2 , H_SSE2 ) );
            F_SSE3             = _mm_add_ps(F_SSE3, _mm_add_ps( G_SSE3 , H_SSE3 ) );
            
            /* VV=Y+eps*Fp */
            Y_SSE0             = _mm_add_ps(Y_SSE0, _mm_mul_ps(F_SSE0, eps_SSE0));
            Y_SSE1             = _mm_add_ps(Y_SSE1, _mm_mul_ps(F_SSE1, eps_SSE1));
            Y_SSE2             = _mm_add_ps(Y_SSE2, _mm_mul_ps(F_SSE2, eps_SSE2));
            Y_SSE3             = _mm_add_ps(Y_SSE3, _mm_mul_ps(F_SSE3, eps_SSE3));

            /* FF = Fp+Geps+2*Heps2 */
            F_SSE0             = _mm_add_ps(F_SSE0, _mm_add_ps(G_SSE0 , _mm_mul_ps(H_SSE0,two_SSE)));
            F_SSE1             = _mm_add_ps(F_SSE1, _mm_add_ps(G_SSE1 , _mm_mul_ps(H_SSE1,two_SSE)));
            F_SSE2             = _mm_add_ps(F_SSE2, _mm_add_ps(G_SSE2 , _mm_mul_ps(H_SSE2,two_SSE)));
            F_SSE3             = _mm_add_ps(F_SSE3, _mm_add_ps(G_SSE3 , _mm_mul_ps(H_SSE3,two_SSE)));

            vgb_SSE0           = _mm_mul_ps(Y_SSE0, qq_SSE0);
            vgb_SSE1           = _mm_mul_ps(Y_SSE1, qq_SSE1);
            vgb_SSE2           = _mm_mul_ps(Y_SSE2, qq_SSE2);
            vgb_SSE3           = _mm_mul_ps(Y_SSE3, qq_SSE3);
            
            fijGB_SSE0         = _mm_mul_ps(F_SSE0, _mm_mul_ps(qq_SSE0,gbscale_SSE0));
            fijGB_SSE1         = _mm_mul_ps(F_SSE1, _mm_mul_ps(qq_SSE1,gbscale_SSE1));
            fijGB_SSE2         = _mm_mul_ps(F_SSE2, _mm_mul_ps(qq_SSE2,gbscale_SSE2));
            fijGB_SSE3         = _mm_mul_ps(F_SSE3, _mm_mul_ps(qq_SSE3,gbscale_SSE3));
            
            /* Note: this dvdatmp has different sign from the usual c code, saves 1 instruction */
            dvdatmp_SSE0       = _mm_mul_ps(_mm_add_ps(vgb_SSE0, _mm_mul_ps(fijGB_SSE0,r_SSE0)) , half_SSE);
            dvdatmp_SSE1       = _mm_mul_ps(_mm_add_ps(vgb_SSE1, _mm_mul_ps(fijGB_SSE1,r_SSE1)) , half_SSE);
            dvdatmp_SSE2       = _mm_mul_ps(_mm_add_ps(vgb_SSE2, _mm_mul_ps(fijGB_SSE2,r_SSE2)) , half_SSE);
            dvdatmp_SSE3       = _mm_mul_ps(_mm_add_ps(vgb_SSE3, _mm_mul_ps(fijGB_SSE3,r_SSE3)) , half_SSE);
            
            vgbtot_SSE         = _mm_add_ps(vgbtot_SSE, gmx_mm_sum4_ps(vgb_SSE0,vgb_SSE1,vgb_SSE2,vgb_SSE3));

            dvdasum_SSE0       = _mm_sub_ps(dvdasum_SSE0, dvdatmp_SSE0);
            dvdasum_SSE1       = _mm_sub_ps(dvdasum_SSE1, dvdatmp_SSE1);
            dvdasum_SSE2       = _mm_sub_ps(dvdasum_SSE2, dvdatmp_SSE2);
            dvdasum_SSE3       = _mm_sub_ps(dvdasum_SSE3, dvdatmp_SSE3);

            dvdatmp_SSE0       = gmx_mm_sum4_ps(dvdatmp_SSE0,dvdatmp_SSE1,dvdatmp_SSE2,dvdatmp_SSE3);
            dvdatmp_SSE0       = _mm_mul_ps(dvdatmp_SSE0, _mm_mul_ps(isaj_SSE,isaj_SSE));

            /* update derivative wrt j atom born radius */
            _mm_store_ps(dvda_align+j,
                         _mm_sub_ps( _mm_load_ps(dvda_align+j) , dvdatmp_SSE0 ));
             
            /* Lennard-Jones interaction */
            rinvsix_SSE0       = _mm_mul_ps(rinvsq_SSE0,_mm_mul_ps(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1       = _mm_mul_ps(rinvsq_SSE1,_mm_mul_ps(rinvsq_SSE1,rinvsq_SSE1));
            rinvsix_SSE2       = _mm_mul_ps(rinvsq_SSE2,_mm_mul_ps(rinvsq_SSE2,rinvsq_SSE2));
            rinvsix_SSE3       = _mm_mul_ps(rinvsq_SSE3,_mm_mul_ps(rinvsq_SSE3,rinvsq_SSE3));
            vvdw6_SSE0         = _mm_mul_ps(c6_SSE0,rinvsix_SSE0);
            vvdw6_SSE1         = _mm_mul_ps(c6_SSE1,rinvsix_SSE1);
            vvdw6_SSE2         = _mm_mul_ps(c6_SSE2,rinvsix_SSE2);
            vvdw6_SSE3         = _mm_mul_ps(c6_SSE3,rinvsix_SSE3);
            vvdw12_SSE0        = _mm_mul_ps(c12_SSE0,_mm_mul_ps(rinvsix_SSE0,rinvsix_SSE0));
            vvdw12_SSE1        = _mm_mul_ps(c12_SSE1,_mm_mul_ps(rinvsix_SSE1,rinvsix_SSE1));
            vvdw12_SSE2        = _mm_mul_ps(c12_SSE2,_mm_mul_ps(rinvsix_SSE2,rinvsix_SSE2));
            vvdw12_SSE3        = _mm_mul_ps(c12_SSE3,_mm_mul_ps(rinvsix_SSE3,rinvsix_SSE3));

            vvdwtot_SSE        = _mm_add_ps(vvdwtot_SSE, gmx_mm_sum4_ps(_mm_sub_ps(vvdw12_SSE0,vvdw6_SSE0),
                                                                     _mm_sub_ps(vvdw12_SSE1,vvdw6_SSE1),
                                                                     _mm_sub_ps(vvdw12_SSE2,vvdw6_SSE2),
                                                                     _mm_sub_ps(vvdw12_SSE3,vvdw6_SSE3)));
            

            fscal_SSE0         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE0, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE0),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE0))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE0,fscal_SSE0),rinv_SSE0 ));
            
            fscal_SSE1         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE1, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE1),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE1))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE1,fscal_SSE1),rinv_SSE1 ));

            fscal_SSE2         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE2, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE2),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE2))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE2,fscal_SSE2),rinv_SSE2 ));

            fscal_SSE3         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE3, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE3),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE3))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE3,fscal_SSE3),rinv_SSE3 ));
            
            /* Calculate temporary vectorial force */
            tx_SSE0            = _mm_mul_ps(fscal_SSE0,dx_SSE0);
            tx_SSE1            = _mm_mul_ps(fscal_SSE1,dx_SSE1);
            tx_SSE2            = _mm_mul_ps(fscal_SSE2,dx_SSE2);
            tx_SSE3            = _mm_mul_ps(fscal_SSE3,dx_SSE3);
            ty_SSE0            = _mm_mul_ps(fscal_SSE0,dy_SSE0);
            ty_SSE1            = _mm_mul_ps(fscal_SSE1,dy_SSE1);
            ty_SSE2            = _mm_mul_ps(fscal_SSE2,dy_SSE2);
            ty_SSE3            = _mm_mul_ps(fscal_SSE3,dy_SSE3);
            tz_SSE0            = _mm_mul_ps(fscal_SSE0,dz_SSE0);
            tz_SSE1            = _mm_mul_ps(fscal_SSE1,dz_SSE1);
            tz_SSE2            = _mm_mul_ps(fscal_SSE2,dz_SSE2);
            tz_SSE3            = _mm_mul_ps(fscal_SSE3,dz_SSE3);
            
            /* Increment i atom force */
            fix_SSE0          = _mm_add_ps(fix_SSE0,tx_SSE0);
            fix_SSE1          = _mm_add_ps(fix_SSE1,tx_SSE1);
            fix_SSE2          = _mm_add_ps(fix_SSE2,tx_SSE2);
            fix_SSE3          = _mm_add_ps(fix_SSE3,tx_SSE3);
            fiy_SSE0          = _mm_add_ps(fiy_SSE0,ty_SSE0);
            fiy_SSE1          = _mm_add_ps(fiy_SSE1,ty_SSE1);
            fiy_SSE2          = _mm_add_ps(fiy_SSE2,ty_SSE2);
            fiy_SSE3          = _mm_add_ps(fiy_SSE3,ty_SSE3);
            fiz_SSE0          = _mm_add_ps(fiz_SSE0,tz_SSE0);
            fiz_SSE1          = _mm_add_ps(fiz_SSE1,tz_SSE1);
            fiz_SSE2          = _mm_add_ps(fiz_SSE2,tz_SSE2);
            fiz_SSE3          = _mm_add_ps(fiz_SSE3,tz_SSE3);
            
            /* Decrement j atom force */
            _mm_store_ps(fx_align+j,
                         _mm_sub_ps( _mm_load_ps(fx_align+j) , gmx_mm_sum4_ps(tx_SSE0,tx_SSE1,tx_SSE2,tx_SSE3) ));
            _mm_store_ps(fy_align+j,
                         _mm_sub_ps( _mm_load_ps(fy_align+j) , gmx_mm_sum4_ps(ty_SSE0,ty_SSE1,ty_SSE2,ty_SSE3) ));
            _mm_store_ps(fz_align+j,
                         _mm_sub_ps( _mm_load_ps(fz_align+j) , gmx_mm_sum4_ps(tz_SSE0,tz_SSE1,tz_SSE2,tz_SSE3) ));
            
            /* Inner loop uses 38 flops/iteration */
        }

        for(j=nj1; j<nj2; j+=UNROLLJ)
        {                      
            /* load j atom coordinates */
            jx_SSE           = _mm_load_ps(x_align+j);
            jy_SSE           = _mm_load_ps(y_align+j);
            jz_SSE           = _mm_load_ps(z_align+j);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0,jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0,jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0,jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1,jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1,jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1,jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2,jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2,jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2,jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3,jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3,jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3,jz_SSE);
            
            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0,dy_SSE0,dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1,dy_SSE1,dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2,dy_SSE2,dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3,dy_SSE3,dz_SSE3);
                        
            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);
            
            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0,imask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1,imask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2,imask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3,imask_SSE3);
            
            /* Load parameters for j atom */
            jq_SSE             = _mm_load_ps(q_align+j);
            qq_SSE0            = _mm_mul_ps(iq_SSE0,jq_SSE);
            qq_SSE1            = _mm_mul_ps(iq_SSE1,jq_SSE);
            qq_SSE2            = _mm_mul_ps(iq_SSE2,jq_SSE);
            qq_SSE3            = _mm_mul_ps(iq_SSE3,jq_SSE);
            qq_SSE0            = _mm_and_ps(qq_SSE0,imask_SSE0);
            qq_SSE1            = _mm_and_ps(qq_SSE1,imask_SSE1);
            qq_SSE2            = _mm_and_ps(qq_SSE2,imask_SSE2);
            qq_SSE3            = _mm_and_ps(qq_SSE3,imask_SSE3);
            
            c6_SSE0            = _mm_load_ps(pvdw0+2*j);
            c6_SSE1            = _mm_load_ps(pvdw1+2*j);
            c6_SSE2            = _mm_load_ps(pvdw2+2*j);
            c6_SSE3            = _mm_load_ps(pvdw3+2*j);
            c12_SSE0           = _mm_load_ps(pvdw0+2*j+UNROLLJ);
            c12_SSE1           = _mm_load_ps(pvdw1+2*j+UNROLLJ);
            c12_SSE2           = _mm_load_ps(pvdw2+2*j+UNROLLJ);
            c12_SSE3           = _mm_load_ps(pvdw3+2*j+UNROLLJ);
            
            isaj_SSE           = _mm_load_ps(invsqrta_align+j);

            isaprod_SSE0       = _mm_mul_ps(isai_SSE0,isaj_SSE);
            isaprod_SSE1       = _mm_mul_ps(isai_SSE1,isaj_SSE);
            isaprod_SSE2       = _mm_mul_ps(isai_SSE2,isaj_SSE);
            isaprod_SSE3       = _mm_mul_ps(isai_SSE3,isaj_SSE);
            
            rinvsq_SSE0        = _mm_mul_ps(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1        = _mm_mul_ps(rinv_SSE1,rinv_SSE1);
            rinvsq_SSE2        = _mm_mul_ps(rinv_SSE2,rinv_SSE2);
            rinvsq_SSE3        = _mm_mul_ps(rinv_SSE3,rinv_SSE3);
            
            /* Coulomb interaction */
            vcoul_SSE0         = _mm_mul_ps(qq_SSE0,rinv_SSE0);
            vcoul_SSE1         = _mm_mul_ps(qq_SSE1,rinv_SSE1);
            vcoul_SSE2         = _mm_mul_ps(qq_SSE2,rinv_SSE2);
            vcoul_SSE3         = _mm_mul_ps(qq_SSE3,rinv_SSE3);
            
            fscal_SSE0         = _mm_mul_ps(vcoul_SSE0,rinv_SSE0);
            fscal_SSE1         = _mm_mul_ps(vcoul_SSE1,rinv_SSE1);
            fscal_SSE2         = _mm_mul_ps(vcoul_SSE2,rinv_SSE2);
            fscal_SSE3         = _mm_mul_ps(vcoul_SSE3,rinv_SSE3);
            
            vctot_SSE          = _mm_add_ps(vctot_SSE, gmx_mm_sum4_ps(vcoul_SSE0,vcoul_SSE1,vcoul_SSE2,vcoul_SSE3));
            
            /* Polarization interaction */
            qq_SSE0            = _mm_mul_ps(qq_SSE0,isaprod_SSE0);
            qq_SSE1            = _mm_mul_ps(qq_SSE1,isaprod_SSE1);
            qq_SSE2            = _mm_mul_ps(qq_SSE2,isaprod_SSE2);
            qq_SSE3            = _mm_mul_ps(qq_SSE3,isaprod_SSE3);
            qq_SSE0            = _mm_mul_ps(qq_SSE0,gbfactor_SSE);
            qq_SSE1            = _mm_mul_ps(qq_SSE1,gbfactor_SSE);
            qq_SSE2            = _mm_mul_ps(qq_SSE2,gbfactor_SSE);
            qq_SSE3            = _mm_mul_ps(qq_SSE3,gbfactor_SSE);
            
            gbscale_SSE0       = _mm_mul_ps(isaprod_SSE0,gbtabscale_SSE);
            gbscale_SSE1       = _mm_mul_ps(isaprod_SSE1,gbtabscale_SSE);
            gbscale_SSE2       = _mm_mul_ps(isaprod_SSE2,gbtabscale_SSE);
            gbscale_SSE3       = _mm_mul_ps(isaprod_SSE3,gbtabscale_SSE);
            
            /* Use tables for (rescaled) GB potential */
            r_SSE0             = _mm_mul_ps(rsq_SSE0,rinv_SSE0);
            r_SSE1             = _mm_mul_ps(rsq_SSE1,rinv_SSE1);
            r_SSE2             = _mm_mul_ps(rsq_SSE2,rinv_SSE2);
            r_SSE3             = _mm_mul_ps(rsq_SSE3,rinv_SSE3);
            rtab_SSE0          = _mm_mul_ps(r_SSE0,gbscale_SSE0);
            rtab_SSE1          = _mm_mul_ps(r_SSE1,gbscale_SSE1);
            rtab_SSE2          = _mm_mul_ps(r_SSE2,gbscale_SSE2);
            rtab_SSE3          = _mm_mul_ps(r_SSE3,gbscale_SSE3);
            
            n0_SSE0            = _mm_cvttps_epi32(rtab_SSE0);
            n0_SSE1            = _mm_cvttps_epi32(rtab_SSE1);
            n0_SSE2            = _mm_cvttps_epi32(rtab_SSE2);
            n0_SSE3            = _mm_cvttps_epi32(rtab_SSE3);
            eps_SSE0           = _mm_sub_ps(rtab_SSE0 , _mm_cvtepi32_ps(n0_SSE0) );
            eps_SSE1           = _mm_sub_ps(rtab_SSE1 , _mm_cvtepi32_ps(n0_SSE1) );
            eps_SSE2           = _mm_sub_ps(rtab_SSE2 , _mm_cvtepi32_ps(n0_SSE2) );
            eps_SSE3           = _mm_sub_ps(rtab_SSE3 , _mm_cvtepi32_ps(n0_SSE3) );
            nnn_SSE0           = _mm_slli_epi32(n0_SSE0,2);
            nnn_SSE1           = _mm_slli_epi32(n0_SSE1,2);
            nnn_SSE2           = _mm_slli_epi32(n0_SSE2,2);
            nnn_SSE3           = _mm_slli_epi32(n0_SSE3,2);
            
            /* the tables are 16-byte aligned, so we can use _mm_load_ps */			
            Y_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,0)); /* YFGH */
            F_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,1)); /* YFGH */
            G_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,2)); /* YFGH */
            H_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE0,F_SSE0,G_SSE0,H_SSE0);
            
            Y_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,0)); /* YFGH */
            F_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,1)); /* YFGH */
            G_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,2)); /* YFGH */
            H_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE1,F_SSE1,G_SSE1,H_SSE1);
            
            Y_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,0)); /* YFGH */
            F_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,1)); /* YFGH */
            G_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,2)); /* YFGH */
            H_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE2,F_SSE2,G_SSE2,H_SSE2);
            
            Y_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,0)); /* YFGH */
            F_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,1)); /* YFGH */
            G_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,2)); /* YFGH */
            H_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE3,F_SSE3,G_SSE3,H_SSE3);
            
            /* G*eps */
            G_SSE0             = _mm_mul_ps(G_SSE0,eps_SSE0);
            G_SSE1             = _mm_mul_ps(G_SSE1,eps_SSE1);
            G_SSE2             = _mm_mul_ps(G_SSE2,eps_SSE2);
            G_SSE3             = _mm_mul_ps(G_SSE3,eps_SSE3);
            
            /* H*eps2 */
            H_SSE0             = _mm_mul_ps(H_SSE0, _mm_mul_ps(eps_SSE0,eps_SSE0) );
            H_SSE1             = _mm_mul_ps(H_SSE1, _mm_mul_ps(eps_SSE1,eps_SSE1) );
            H_SSE2             = _mm_mul_ps(H_SSE2, _mm_mul_ps(eps_SSE2,eps_SSE2) );
            H_SSE3             = _mm_mul_ps(H_SSE3, _mm_mul_ps(eps_SSE3,eps_SSE3) );
            
            /* Fp=F+G*eps+H*eps2 */
            F_SSE0             = _mm_add_ps(F_SSE0, _mm_add_ps( G_SSE0 , H_SSE0 ) );
            F_SSE1             = _mm_add_ps(F_SSE1, _mm_add_ps( G_SSE1 , H_SSE1 ) );
            F_SSE2             = _mm_add_ps(F_SSE2, _mm_add_ps( G_SSE2 , H_SSE2 ) );
            F_SSE3             = _mm_add_ps(F_SSE3, _mm_add_ps( G_SSE3 , H_SSE3 ) );
            
            /* VV=Y+eps*Fp */
            Y_SSE0             = _mm_add_ps(Y_SSE0, _mm_mul_ps(F_SSE0, eps_SSE0));
            Y_SSE1             = _mm_add_ps(Y_SSE1, _mm_mul_ps(F_SSE1, eps_SSE1));
            Y_SSE2             = _mm_add_ps(Y_SSE2, _mm_mul_ps(F_SSE2, eps_SSE2));
            Y_SSE3             = _mm_add_ps(Y_SSE3, _mm_mul_ps(F_SSE3, eps_SSE3));
            
            /* FF = Fp+Geps+2*Heps2 */
            F_SSE0             = _mm_add_ps(F_SSE0, _mm_add_ps(G_SSE0 , _mm_mul_ps(H_SSE0,two_SSE)));
            F_SSE1             = _mm_add_ps(F_SSE1, _mm_add_ps(G_SSE1 , _mm_mul_ps(H_SSE1,two_SSE)));
            F_SSE2             = _mm_add_ps(F_SSE2, _mm_add_ps(G_SSE2 , _mm_mul_ps(H_SSE2,two_SSE)));
            F_SSE3             = _mm_add_ps(F_SSE3, _mm_add_ps(G_SSE3 , _mm_mul_ps(H_SSE3,two_SSE)));
            
            vgb_SSE0           = _mm_mul_ps(Y_SSE0, qq_SSE0);
            vgb_SSE1           = _mm_mul_ps(Y_SSE1, qq_SSE1);
            vgb_SSE2           = _mm_mul_ps(Y_SSE2, qq_SSE2);
            vgb_SSE3           = _mm_mul_ps(Y_SSE3, qq_SSE3);
            
            fijGB_SSE0         = _mm_mul_ps(F_SSE0, _mm_mul_ps(qq_SSE0,gbscale_SSE0));
            fijGB_SSE1         = _mm_mul_ps(F_SSE1, _mm_mul_ps(qq_SSE1,gbscale_SSE1));
            fijGB_SSE2         = _mm_mul_ps(F_SSE2, _mm_mul_ps(qq_SSE2,gbscale_SSE2));
            fijGB_SSE3         = _mm_mul_ps(F_SSE3, _mm_mul_ps(qq_SSE3,gbscale_SSE3));
        
            
            /* Note: this dvdatmp has different sign from the usual c code, saves 1 instruction */
            dvdatmp_SSE0       = _mm_mul_ps(_mm_add_ps(vgb_SSE0, _mm_mul_ps(fijGB_SSE0,r_SSE0)) , half_SSE);
            dvdatmp_SSE1       = _mm_mul_ps(_mm_add_ps(vgb_SSE1, _mm_mul_ps(fijGB_SSE1,r_SSE1)) , half_SSE);
            dvdatmp_SSE2       = _mm_mul_ps(_mm_add_ps(vgb_SSE2, _mm_mul_ps(fijGB_SSE2,r_SSE2)) , half_SSE);
            dvdatmp_SSE3       = _mm_mul_ps(_mm_add_ps(vgb_SSE3, _mm_mul_ps(fijGB_SSE3,r_SSE3)) , half_SSE);

            vgbtot_SSE         = _mm_add_ps(vgbtot_SSE, gmx_mm_sum4_ps(vgb_SSE0,vgb_SSE1,vgb_SSE2,vgb_SSE3));
            
            dvdasum_SSE0       = _mm_sub_ps(dvdasum_SSE0, dvdatmp_SSE0);
            dvdasum_SSE1       = _mm_sub_ps(dvdasum_SSE1, dvdatmp_SSE1);
            dvdasum_SSE2       = _mm_sub_ps(dvdasum_SSE2, dvdatmp_SSE2);
            dvdasum_SSE3       = _mm_sub_ps(dvdasum_SSE3, dvdatmp_SSE3);
            
            dvdatmp_SSE0       = gmx_mm_sum4_ps(dvdatmp_SSE0,dvdatmp_SSE1,dvdatmp_SSE2,dvdatmp_SSE3);
            dvdatmp_SSE0       = _mm_mul_ps(dvdatmp_SSE0, _mm_mul_ps(isaj_SSE,isaj_SSE));
            
            /* update derivative wrt j atom born radius */
            _mm_store_ps(dvda_align+j,
                         _mm_sub_ps( _mm_load_ps(dvda_align+j) , dvdatmp_SSE0 ));
            
            /* Lennard-Jones interaction */
            rinvsix_SSE0       = _mm_mul_ps(rinvsq_SSE0,_mm_mul_ps(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1       = _mm_mul_ps(rinvsq_SSE1,_mm_mul_ps(rinvsq_SSE1,rinvsq_SSE1));
            rinvsix_SSE2       = _mm_mul_ps(rinvsq_SSE2,_mm_mul_ps(rinvsq_SSE2,rinvsq_SSE2));
            rinvsix_SSE3       = _mm_mul_ps(rinvsq_SSE3,_mm_mul_ps(rinvsq_SSE3,rinvsq_SSE3));
            vvdw6_SSE0         = _mm_mul_ps(c6_SSE0,rinvsix_SSE0);
            vvdw6_SSE1         = _mm_mul_ps(c6_SSE1,rinvsix_SSE1);
            vvdw6_SSE2         = _mm_mul_ps(c6_SSE2,rinvsix_SSE2);
            vvdw6_SSE3         = _mm_mul_ps(c6_SSE3,rinvsix_SSE3);
            vvdw12_SSE0        = _mm_mul_ps(c12_SSE0,_mm_mul_ps(rinvsix_SSE0,rinvsix_SSE0));
            vvdw12_SSE1        = _mm_mul_ps(c12_SSE1,_mm_mul_ps(rinvsix_SSE1,rinvsix_SSE1));
            vvdw12_SSE2        = _mm_mul_ps(c12_SSE2,_mm_mul_ps(rinvsix_SSE2,rinvsix_SSE2));
            vvdw12_SSE3        = _mm_mul_ps(c12_SSE3,_mm_mul_ps(rinvsix_SSE3,rinvsix_SSE3));
            
            vvdwtot_SSE        = _mm_add_ps(vvdwtot_SSE, gmx_mm_sum4_ps(_mm_sub_ps(vvdw12_SSE0,vvdw6_SSE0),
                                                                     _mm_sub_ps(vvdw12_SSE1,vvdw6_SSE1),
                                                                     _mm_sub_ps(vvdw12_SSE2,vvdw6_SSE2),
                                                                     _mm_sub_ps(vvdw12_SSE3,vvdw6_SSE3)));
            
            fscal_SSE0         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE0, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE0),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE0))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE0,fscal_SSE0),rinv_SSE0 ));
            
            fscal_SSE1         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE1, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE1),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE1))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE1,fscal_SSE1),rinv_SSE1 ));
            
            fscal_SSE2         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE2, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE2),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE2))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE2,fscal_SSE2),rinv_SSE2 ));
            
            fscal_SSE3         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE3, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE3),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE3))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE3,fscal_SSE3),rinv_SSE3 ));
                        
            /* Calculate temporary vectorial force */
            tx_SSE0            = _mm_mul_ps(fscal_SSE0,dx_SSE0);
            tx_SSE1            = _mm_mul_ps(fscal_SSE1,dx_SSE1);
            tx_SSE2            = _mm_mul_ps(fscal_SSE2,dx_SSE2);
            tx_SSE3            = _mm_mul_ps(fscal_SSE3,dx_SSE3);
            ty_SSE0            = _mm_mul_ps(fscal_SSE0,dy_SSE0);
            ty_SSE1            = _mm_mul_ps(fscal_SSE1,dy_SSE1);
            ty_SSE2            = _mm_mul_ps(fscal_SSE2,dy_SSE2);
            ty_SSE3            = _mm_mul_ps(fscal_SSE3,dy_SSE3);
            tz_SSE0            = _mm_mul_ps(fscal_SSE0,dz_SSE0);
            tz_SSE1            = _mm_mul_ps(fscal_SSE1,dz_SSE1);
            tz_SSE2            = _mm_mul_ps(fscal_SSE2,dz_SSE2);
            tz_SSE3            = _mm_mul_ps(fscal_SSE3,dz_SSE3);
            
            /* Increment i atom force */
            fix_SSE0          = _mm_add_ps(fix_SSE0,tx_SSE0);
            fix_SSE1          = _mm_add_ps(fix_SSE1,tx_SSE1);
            fix_SSE2          = _mm_add_ps(fix_SSE2,tx_SSE2);
            fix_SSE3          = _mm_add_ps(fix_SSE3,tx_SSE3);
            fiy_SSE0          = _mm_add_ps(fiy_SSE0,ty_SSE0);
            fiy_SSE1          = _mm_add_ps(fiy_SSE1,ty_SSE1);
            fiy_SSE2          = _mm_add_ps(fiy_SSE2,ty_SSE2);
            fiy_SSE3          = _mm_add_ps(fiy_SSE3,ty_SSE3);
            fiz_SSE0          = _mm_add_ps(fiz_SSE0,tz_SSE0);
            fiz_SSE1          = _mm_add_ps(fiz_SSE1,tz_SSE1);
            fiz_SSE2          = _mm_add_ps(fiz_SSE2,tz_SSE2);
            fiz_SSE3          = _mm_add_ps(fiz_SSE3,tz_SSE3);
            
            /* Decrement j atom force */
            _mm_store_ps(fx_align+j,
                         _mm_sub_ps( _mm_load_ps(fx_align+j) , gmx_mm_sum4_ps(tx_SSE0,tx_SSE1,tx_SSE2,tx_SSE3) ));
            _mm_store_ps(fy_align+j,
                         _mm_sub_ps( _mm_load_ps(fy_align+j) , gmx_mm_sum4_ps(ty_SSE0,ty_SSE1,ty_SSE2,ty_SSE3) ));
            _mm_store_ps(fz_align+j,
                         _mm_sub_ps( _mm_load_ps(fz_align+j) , gmx_mm_sum4_ps(tz_SSE0,tz_SSE1,tz_SSE2,tz_SSE3) ));
            /* Inner loop uses 38 flops/iteration */
        }

        for(j=nj2; j<nj3; j+=UNROLLJ)
        {
            jmask_SSE0 = _mm_load_ps((real *)emask0);
            jmask_SSE1 = _mm_load_ps((real *)emask1);
            jmask_SSE2 = _mm_load_ps((real *)emask2);
            jmask_SSE3 = _mm_load_ps((real *)emask3);
            emask0 += UNROLLJ;
            emask1 += UNROLLJ;
            emask2 += UNROLLJ;
            emask3 += UNROLLJ;

            /* load j atom coordinates */
            jx_SSE           = _mm_load_ps(x_align+j);
            jy_SSE           = _mm_load_ps(y_align+j);
            jz_SSE           = _mm_load_ps(z_align+j);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_ps(ix_SSE0,jx_SSE);
            dy_SSE0            = _mm_sub_ps(iy_SSE0,jy_SSE);
            dz_SSE0            = _mm_sub_ps(iz_SSE0,jz_SSE);
            dx_SSE1            = _mm_sub_ps(ix_SSE1,jx_SSE);
            dy_SSE1            = _mm_sub_ps(iy_SSE1,jy_SSE);
            dz_SSE1            = _mm_sub_ps(iz_SSE1,jz_SSE);
            dx_SSE2            = _mm_sub_ps(ix_SSE2,jx_SSE);
            dy_SSE2            = _mm_sub_ps(iy_SSE2,jy_SSE);
            dz_SSE2            = _mm_sub_ps(iz_SSE2,jz_SSE);
            dx_SSE3            = _mm_sub_ps(ix_SSE3,jx_SSE);
            dy_SSE3            = _mm_sub_ps(iy_SSE3,jy_SSE);
            dz_SSE3            = _mm_sub_ps(iz_SSE3,jz_SSE);
            
            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_ps(dx_SSE0,dy_SSE0,dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_ps(dx_SSE1,dy_SSE1,dz_SSE1);
            rsq_SSE2           = gmx_mm_calc_rsq_ps(dx_SSE2,dy_SSE2,dz_SSE2);
            rsq_SSE3           = gmx_mm_calc_rsq_ps(dx_SSE3,dy_SSE3,dz_SSE3);
            
            /* Combine masks */
            jmask_SSE0         = _mm_and_ps(jmask_SSE0,imask_SSE0);
            jmask_SSE1         = _mm_and_ps(jmask_SSE1,imask_SSE1);
            jmask_SSE2         = _mm_and_ps(jmask_SSE2,imask_SSE2);
            jmask_SSE3         = _mm_and_ps(jmask_SSE3,imask_SSE3);
            
            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_ps(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_ps(rsq_SSE1);
            rinv_SSE2          = gmx_mm_invsqrt_ps(rsq_SSE2);
            rinv_SSE3          = gmx_mm_invsqrt_ps(rsq_SSE3);
            
            /* Apply mask */
            rinv_SSE0          = _mm_and_ps(rinv_SSE0,jmask_SSE0);
            rinv_SSE1          = _mm_and_ps(rinv_SSE1,jmask_SSE1);
            rinv_SSE2          = _mm_and_ps(rinv_SSE2,jmask_SSE2);
            rinv_SSE3          = _mm_and_ps(rinv_SSE3,jmask_SSE3);
            
            /* Load parameters for j atom */
            jq_SSE             = _mm_load_ps(q_align+j);
            qq_SSE0            = _mm_mul_ps(iq_SSE0,jq_SSE);
            qq_SSE1            = _mm_mul_ps(iq_SSE1,jq_SSE);
            qq_SSE2            = _mm_mul_ps(iq_SSE2,jq_SSE);
            qq_SSE3            = _mm_mul_ps(iq_SSE3,jq_SSE);
            qq_SSE0            = _mm_and_ps(qq_SSE0,jmask_SSE0);
            qq_SSE1            = _mm_and_ps(qq_SSE1,jmask_SSE1);
            qq_SSE2            = _mm_and_ps(qq_SSE2,jmask_SSE2);
            qq_SSE3            = _mm_and_ps(qq_SSE3,jmask_SSE3);
            
            c6_SSE0            = _mm_load_ps(pvdw0+2*j);
            c6_SSE1            = _mm_load_ps(pvdw1+2*j);
            c6_SSE2            = _mm_load_ps(pvdw2+2*j);
            c6_SSE3            = _mm_load_ps(pvdw3+2*j);
            c12_SSE0           = _mm_load_ps(pvdw0+2*j+UNROLLJ);
            c12_SSE1           = _mm_load_ps(pvdw1+2*j+UNROLLJ);
            c12_SSE2           = _mm_load_ps(pvdw2+2*j+UNROLLJ);
            c12_SSE3           = _mm_load_ps(pvdw3+2*j+UNROLLJ);
            
            isaj_SSE           = _mm_load_ps(invsqrta_align+j);

            isaprod_SSE0       = _mm_mul_ps(isai_SSE0,isaj_SSE);
            isaprod_SSE1       = _mm_mul_ps(isai_SSE1,isaj_SSE);
            isaprod_SSE2       = _mm_mul_ps(isai_SSE2,isaj_SSE);
            isaprod_SSE3       = _mm_mul_ps(isai_SSE3,isaj_SSE);
            
            rinvsq_SSE0        = _mm_mul_ps(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1        = _mm_mul_ps(rinv_SSE1,rinv_SSE1);
            rinvsq_SSE2        = _mm_mul_ps(rinv_SSE2,rinv_SSE2);
            rinvsq_SSE3        = _mm_mul_ps(rinv_SSE3,rinv_SSE3);
            
            /* Coulomb interaction */
            vcoul_SSE0         = _mm_mul_ps(qq_SSE0,rinv_SSE0);
            vcoul_SSE1         = _mm_mul_ps(qq_SSE1,rinv_SSE1);
            vcoul_SSE2         = _mm_mul_ps(qq_SSE2,rinv_SSE2);
            vcoul_SSE3         = _mm_mul_ps(qq_SSE3,rinv_SSE3);
            
            fscal_SSE0         = _mm_mul_ps(vcoul_SSE0,rinv_SSE0);
            fscal_SSE1         = _mm_mul_ps(vcoul_SSE1,rinv_SSE1);
            fscal_SSE2         = _mm_mul_ps(vcoul_SSE2,rinv_SSE2);
            fscal_SSE3         = _mm_mul_ps(vcoul_SSE3,rinv_SSE3);
            
            vctot_SSE          = _mm_add_ps(vctot_SSE, gmx_mm_sum4_ps(vcoul_SSE0,vcoul_SSE1,vcoul_SSE2,vcoul_SSE3));
            
            /* Polarization interaction */
            qq_SSE0            = _mm_mul_ps(qq_SSE0,isaprod_SSE0);
            qq_SSE1            = _mm_mul_ps(qq_SSE1,isaprod_SSE1);
            qq_SSE2            = _mm_mul_ps(qq_SSE2,isaprod_SSE2);
            qq_SSE3            = _mm_mul_ps(qq_SSE3,isaprod_SSE3);
            qq_SSE0            = _mm_mul_ps(qq_SSE0,gbfactor_SSE);
            qq_SSE1            = _mm_mul_ps(qq_SSE1,gbfactor_SSE);
            qq_SSE2            = _mm_mul_ps(qq_SSE2,gbfactor_SSE);
            qq_SSE3            = _mm_mul_ps(qq_SSE3,gbfactor_SSE);
            
            gbscale_SSE0       = _mm_mul_ps(isaprod_SSE0,gbtabscale_SSE);
            gbscale_SSE1       = _mm_mul_ps(isaprod_SSE1,gbtabscale_SSE);
            gbscale_SSE2       = _mm_mul_ps(isaprod_SSE2,gbtabscale_SSE);
            gbscale_SSE3       = _mm_mul_ps(isaprod_SSE3,gbtabscale_SSE);
            
            /* Use tables for (rescaled) GB potential */
            r_SSE0             = _mm_mul_ps(rsq_SSE0,rinv_SSE0);
            r_SSE1             = _mm_mul_ps(rsq_SSE1,rinv_SSE1);
            r_SSE2             = _mm_mul_ps(rsq_SSE2,rinv_SSE2);
            r_SSE3             = _mm_mul_ps(rsq_SSE3,rinv_SSE3);
            rtab_SSE0          = _mm_mul_ps(r_SSE0,gbscale_SSE0);
            rtab_SSE1          = _mm_mul_ps(r_SSE1,gbscale_SSE1);
            rtab_SSE2          = _mm_mul_ps(r_SSE2,gbscale_SSE2);
            rtab_SSE3          = _mm_mul_ps(r_SSE3,gbscale_SSE3);
            
            n0_SSE0            = _mm_cvttps_epi32(rtab_SSE0);
            n0_SSE1            = _mm_cvttps_epi32(rtab_SSE1);
            n0_SSE2            = _mm_cvttps_epi32(rtab_SSE2);
            n0_SSE3            = _mm_cvttps_epi32(rtab_SSE3);
            eps_SSE0           = _mm_sub_ps(rtab_SSE0 , _mm_cvtepi32_ps(n0_SSE0) );
            eps_SSE1           = _mm_sub_ps(rtab_SSE1 , _mm_cvtepi32_ps(n0_SSE1) );
            eps_SSE2           = _mm_sub_ps(rtab_SSE2 , _mm_cvtepi32_ps(n0_SSE2) );
            eps_SSE3           = _mm_sub_ps(rtab_SSE3 , _mm_cvtepi32_ps(n0_SSE3) );
            nnn_SSE0           = _mm_slli_epi32(n0_SSE0,2);
            nnn_SSE1           = _mm_slli_epi32(n0_SSE1,2);
            nnn_SSE2           = _mm_slli_epi32(n0_SSE2,2);
            nnn_SSE3           = _mm_slli_epi32(n0_SSE3,2);
            
            /* the tables are 16-byte aligned, so we can use _mm_load_ps */			
            Y_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,0)); /* YFGH */
            F_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,1)); /* YFGH */
            G_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,2)); /* YFGH */
            H_SSE0             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE0,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE0,F_SSE0,G_SSE0,H_SSE0);
            
            Y_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,0)); /* YFGH */
            F_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,1)); /* YFGH */
            G_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,2)); /* YFGH */
            H_SSE1             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE1,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE1,F_SSE1,G_SSE1,H_SSE1);
            
            Y_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,0)); /* YFGH */
            F_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,1)); /* YFGH */
            G_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,2)); /* YFGH */
            H_SSE2             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE2,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE2,F_SSE2,G_SSE2,H_SSE2);
            
            Y_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,0)); /* YFGH */
            F_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,1)); /* YFGH */
            G_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,2)); /* YFGH */
            H_SSE3             = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn_SSE3,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y_SSE3,F_SSE3,G_SSE3,H_SSE3);
            
            /* G*eps */
            G_SSE0             = _mm_mul_ps(G_SSE0,eps_SSE0);
            G_SSE1             = _mm_mul_ps(G_SSE1,eps_SSE1);
            G_SSE2             = _mm_mul_ps(G_SSE2,eps_SSE2);
            G_SSE3             = _mm_mul_ps(G_SSE3,eps_SSE3);
            
            /* H*eps2 */
            H_SSE0             = _mm_mul_ps(H_SSE0, _mm_mul_ps(eps_SSE0,eps_SSE0) );
            H_SSE1             = _mm_mul_ps(H_SSE1, _mm_mul_ps(eps_SSE1,eps_SSE1) );
            H_SSE2             = _mm_mul_ps(H_SSE2, _mm_mul_ps(eps_SSE2,eps_SSE2) );
            H_SSE3             = _mm_mul_ps(H_SSE3, _mm_mul_ps(eps_SSE3,eps_SSE3) );
            
            /* Fp=F+G*eps+H*eps2 */
            F_SSE0             = _mm_add_ps(F_SSE0, _mm_add_ps( G_SSE0 , H_SSE0 ) );
            F_SSE1             = _mm_add_ps(F_SSE1, _mm_add_ps( G_SSE1 , H_SSE1 ) );
            F_SSE2             = _mm_add_ps(F_SSE2, _mm_add_ps( G_SSE2 , H_SSE2 ) );
            F_SSE3             = _mm_add_ps(F_SSE3, _mm_add_ps( G_SSE3 , H_SSE3 ) );
            
            /* VV=Y+eps*Fp */
            Y_SSE0             = _mm_add_ps(Y_SSE0, _mm_mul_ps(F_SSE0, eps_SSE0));
            Y_SSE1             = _mm_add_ps(Y_SSE1, _mm_mul_ps(F_SSE1, eps_SSE1));
            Y_SSE2             = _mm_add_ps(Y_SSE2, _mm_mul_ps(F_SSE2, eps_SSE2));
            Y_SSE3             = _mm_add_ps(Y_SSE3, _mm_mul_ps(F_SSE3, eps_SSE3));
            
            /* FF = Fp+Geps+2*Heps2 */
            F_SSE0             = _mm_add_ps(F_SSE0, _mm_add_ps(G_SSE0 , _mm_mul_ps(H_SSE0,two_SSE)));
            F_SSE1             = _mm_add_ps(F_SSE1, _mm_add_ps(G_SSE1 , _mm_mul_ps(H_SSE1,two_SSE)));
            F_SSE2             = _mm_add_ps(F_SSE2, _mm_add_ps(G_SSE2 , _mm_mul_ps(H_SSE2,two_SSE)));
            F_SSE3             = _mm_add_ps(F_SSE3, _mm_add_ps(G_SSE3 , _mm_mul_ps(H_SSE3,two_SSE)));
            
            vgb_SSE0           = _mm_mul_ps(Y_SSE0, qq_SSE0);
            vgb_SSE1           = _mm_mul_ps(Y_SSE1, qq_SSE1);
            vgb_SSE2           = _mm_mul_ps(Y_SSE2, qq_SSE2);
            vgb_SSE3           = _mm_mul_ps(Y_SSE3, qq_SSE3);
            
            fijGB_SSE0         = _mm_mul_ps(F_SSE0, _mm_mul_ps(qq_SSE0,gbscale_SSE0));
            fijGB_SSE1         = _mm_mul_ps(F_SSE1, _mm_mul_ps(qq_SSE1,gbscale_SSE1));
            fijGB_SSE2         = _mm_mul_ps(F_SSE2, _mm_mul_ps(qq_SSE2,gbscale_SSE2));
            fijGB_SSE3         = _mm_mul_ps(F_SSE3, _mm_mul_ps(qq_SSE3,gbscale_SSE3));
            
            /* Note: this dvdatmp has different sign from the usual c code, saves 1 instruction */
            dvdatmp_SSE0       = _mm_mul_ps(_mm_add_ps(vgb_SSE0, _mm_mul_ps(fijGB_SSE0,r_SSE0)) , half_SSE);
            dvdatmp_SSE1       = _mm_mul_ps(_mm_add_ps(vgb_SSE1, _mm_mul_ps(fijGB_SSE1,r_SSE1)) , half_SSE);
            dvdatmp_SSE2       = _mm_mul_ps(_mm_add_ps(vgb_SSE2, _mm_mul_ps(fijGB_SSE2,r_SSE2)) , half_SSE);
            dvdatmp_SSE3       = _mm_mul_ps(_mm_add_ps(vgb_SSE3, _mm_mul_ps(fijGB_SSE3,r_SSE3)) , half_SSE);
            
            vgbtot_SSE         = _mm_add_ps(vgbtot_SSE, gmx_mm_sum4_ps(vgb_SSE0,vgb_SSE1,vgb_SSE2,vgb_SSE3));
            
            dvdasum_SSE0       = _mm_sub_ps(dvdasum_SSE0, dvdatmp_SSE0);
            dvdasum_SSE1       = _mm_sub_ps(dvdasum_SSE1, dvdatmp_SSE1);
            dvdasum_SSE2       = _mm_sub_ps(dvdasum_SSE2, dvdatmp_SSE2);
            dvdasum_SSE3       = _mm_sub_ps(dvdasum_SSE3, dvdatmp_SSE3);
            
            dvdatmp_SSE0       = gmx_mm_sum4_ps(dvdatmp_SSE0,dvdatmp_SSE1,dvdatmp_SSE2,dvdatmp_SSE3);
            dvdatmp_SSE0       = _mm_mul_ps(dvdatmp_SSE0, _mm_mul_ps(isaj_SSE,isaj_SSE));
            
           /* update derivative wrt j atom born radius */
            _mm_store_ps(dvda_align+j,
                         _mm_sub_ps( _mm_load_ps(dvda_align+j) , dvdatmp_SSE0 ));
            
            /* Lennard-Jones interaction */
            rinvsix_SSE0       = _mm_mul_ps(rinvsq_SSE0,_mm_mul_ps(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1       = _mm_mul_ps(rinvsq_SSE1,_mm_mul_ps(rinvsq_SSE1,rinvsq_SSE1));
            rinvsix_SSE2       = _mm_mul_ps(rinvsq_SSE2,_mm_mul_ps(rinvsq_SSE2,rinvsq_SSE2));
            rinvsix_SSE3       = _mm_mul_ps(rinvsq_SSE3,_mm_mul_ps(rinvsq_SSE3,rinvsq_SSE3));
            vvdw6_SSE0         = _mm_mul_ps(c6_SSE0,rinvsix_SSE0);
            vvdw6_SSE1         = _mm_mul_ps(c6_SSE1,rinvsix_SSE1);
            vvdw6_SSE2         = _mm_mul_ps(c6_SSE2,rinvsix_SSE2);
            vvdw6_SSE3         = _mm_mul_ps(c6_SSE3,rinvsix_SSE3);
            vvdw12_SSE0        = _mm_mul_ps(c12_SSE0,_mm_mul_ps(rinvsix_SSE0,rinvsix_SSE0));
            vvdw12_SSE1        = _mm_mul_ps(c12_SSE1,_mm_mul_ps(rinvsix_SSE1,rinvsix_SSE1));
            vvdw12_SSE2        = _mm_mul_ps(c12_SSE2,_mm_mul_ps(rinvsix_SSE2,rinvsix_SSE2));
            vvdw12_SSE3        = _mm_mul_ps(c12_SSE3,_mm_mul_ps(rinvsix_SSE3,rinvsix_SSE3));
            
            vvdwtot_SSE        = _mm_add_ps(vvdwtot_SSE, gmx_mm_sum4_ps(_mm_sub_ps(vvdw12_SSE0,vvdw6_SSE0),
                                                                     _mm_sub_ps(vvdw12_SSE1,vvdw6_SSE1),
                                                                     _mm_sub_ps(vvdw12_SSE2,vvdw6_SSE2),
                                                                     _mm_sub_ps(vvdw12_SSE3,vvdw6_SSE3)));
            
            fscal_SSE0         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE0, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE0),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE0))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE0,fscal_SSE0),rinv_SSE0 ));
            
            fscal_SSE1         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE1, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE1),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE1))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE1,fscal_SSE1),rinv_SSE1 ));
            
            fscal_SSE2         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE2, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE2),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE2))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE2,fscal_SSE2),rinv_SSE2 ));
            
            fscal_SSE3         = _mm_sub_ps(_mm_mul_ps(rinvsq_SSE3, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve_SSE,vvdw12_SSE3),
                                                                  _mm_mul_ps(six_SSE,vvdw6_SSE3))),
                                            _mm_mul_ps( _mm_sub_ps( fijGB_SSE3,fscal_SSE3),rinv_SSE3 ));

             /* Calculate temporary vectorial force */
            tx_SSE0            = _mm_mul_ps(fscal_SSE0,dx_SSE0);
            tx_SSE1            = _mm_mul_ps(fscal_SSE1,dx_SSE1);
            tx_SSE2            = _mm_mul_ps(fscal_SSE2,dx_SSE2);
            tx_SSE3            = _mm_mul_ps(fscal_SSE3,dx_SSE3);
            ty_SSE0            = _mm_mul_ps(fscal_SSE0,dy_SSE0);
            ty_SSE1            = _mm_mul_ps(fscal_SSE1,dy_SSE1);
            ty_SSE2            = _mm_mul_ps(fscal_SSE2,dy_SSE2);
            ty_SSE3            = _mm_mul_ps(fscal_SSE3,dy_SSE3);
            tz_SSE0            = _mm_mul_ps(fscal_SSE0,dz_SSE0);
            tz_SSE1            = _mm_mul_ps(fscal_SSE1,dz_SSE1);
            tz_SSE2            = _mm_mul_ps(fscal_SSE2,dz_SSE2);
            tz_SSE3            = _mm_mul_ps(fscal_SSE3,dz_SSE3);
            
            /* Increment i atom force */
            fix_SSE0          = _mm_add_ps(fix_SSE0,tx_SSE0);
            fix_SSE1          = _mm_add_ps(fix_SSE1,tx_SSE1);
            fix_SSE2          = _mm_add_ps(fix_SSE2,tx_SSE2);
            fix_SSE3          = _mm_add_ps(fix_SSE3,tx_SSE3);
            fiy_SSE0          = _mm_add_ps(fiy_SSE0,ty_SSE0);
            fiy_SSE1          = _mm_add_ps(fiy_SSE1,ty_SSE1);
            fiy_SSE2          = _mm_add_ps(fiy_SSE2,ty_SSE2);
            fiy_SSE3          = _mm_add_ps(fiy_SSE3,ty_SSE3);
            fiz_SSE0          = _mm_add_ps(fiz_SSE0,tz_SSE0);
            fiz_SSE1          = _mm_add_ps(fiz_SSE1,tz_SSE1);
            fiz_SSE2          = _mm_add_ps(fiz_SSE2,tz_SSE2);
            fiz_SSE3          = _mm_add_ps(fiz_SSE3,tz_SSE3);
            
            /* Decrement j atom force */
            _mm_store_ps(fx_align+j,
                         _mm_sub_ps( _mm_load_ps(fx_align+j) , gmx_mm_sum4_ps(tx_SSE0,tx_SSE1,tx_SSE2,tx_SSE3) ));
            _mm_store_ps(fy_align+j,
                         _mm_sub_ps( _mm_load_ps(fy_align+j) , gmx_mm_sum4_ps(ty_SSE0,ty_SSE1,ty_SSE2,ty_SSE3) ));
            _mm_store_ps(fz_align+j,
                         _mm_sub_ps( _mm_load_ps(fz_align+j) , gmx_mm_sum4_ps(tz_SSE0,tz_SSE1,tz_SSE2,tz_SSE3) ));
            
            /* Inner loop uses 38 flops/iteration */
        }
        
		/* Add i forces to mem and shifted force list */
        _MM_TRANSPOSE4_PS(fix_SSE0,fix_SSE1,fix_SSE2,fix_SSE3);
        fix_SSE0 = _mm_add_ps(fix_SSE0,fix_SSE1);
        fix_SSE2 = _mm_add_ps(fix_SSE2,fix_SSE3);
        fix_SSE0 = _mm_add_ps(fix_SSE0,fix_SSE2);
        _mm_store_ps(fx_align+i, _mm_add_ps(fix_SSE0, _mm_load_ps(fx_align+i)));
        
        _MM_TRANSPOSE4_PS(fiy_SSE0,fiy_SSE1,fiy_SSE2,fiy_SSE3);
        fiy_SSE0 = _mm_add_ps(fiy_SSE0,fiy_SSE1);
        fiy_SSE2 = _mm_add_ps(fiy_SSE2,fiy_SSE3);
        fiy_SSE0 = _mm_add_ps(fiy_SSE0,fiy_SSE2);
        _mm_store_ps(fy_align+i, _mm_add_ps(fiy_SSE0, _mm_load_ps(fy_align+i)));
        
        _MM_TRANSPOSE4_PS(fiz_SSE0,fiz_SSE1,fiz_SSE2,fiz_SSE3);
        fiz_SSE0 = _mm_add_ps(fiz_SSE0,fiz_SSE1);
        fiz_SSE2 = _mm_add_ps(fiz_SSE2,fiz_SSE3);
        fiz_SSE0 = _mm_add_ps(fiz_SSE0,fiz_SSE2);
        _mm_store_ps(fz_align+i, _mm_add_ps(fiz_SSE0, _mm_load_ps(fz_align+i)));
		
        _MM_TRANSPOSE4_PS(dvdasum_SSE0,dvdasum_SSE1,dvdasum_SSE2,dvdasum_SSE3);
        dvdasum_SSE0 = _mm_add_ps(dvdasum_SSE0,dvdasum_SSE1);
        dvdasum_SSE2 = _mm_add_ps(dvdasum_SSE2,dvdasum_SSE3);
        dvdasum_SSE0 = _mm_add_ps(dvdasum_SSE0,dvdasum_SSE2);
        dvdasum_SSE0 = _mm_mul_ps(dvdasum_SSE0, _mm_mul_ps(isai_SSE,isai_SSE));
        _mm_store_ps(dvda_align+i, _mm_add_ps(dvdasum_SSE0, _mm_load_ps(dvda_align+i)));

		/* Add potential energies to the group for this list */
		ggid             = 0;         
        
		_mm_storeu_ps(tmpsum,vctot_SSE);
		vc[ggid]         = vc[ggid] + tmpsum[0]+tmpsum[1]+tmpsum[2]+tmpsum[3];
		
		_mm_storeu_ps(tmpsum,vvdwtot_SSE);
		vvdw[ggid]       = vvdw[ggid] + tmpsum[0]+tmpsum[1]+tmpsum[2]+tmpsum[3];
	
        _mm_storeu_ps(tmpsum,vgbtot_SSE);
		vpol[ggid]       = vpol[ggid] + tmpsum[0]+tmpsum[1]+tmpsum[2]+tmpsum[3];
		
		/* Outer loop uses 6 flops/iteration */
	}    

	for(i=ni0;i<ni1+1+natoms/2;i++)
	{
        k = i%natoms;
		f[3*k]   += fx_align[i];
		f[3*k+1] += fy_align[i];
		f[3*k+2] += fz_align[i];
        fr->dvda[k] += dvda_align[i];
	}
    
    /* Write outer/inner iteration count to pointers */
    *outeriter       = ni1-ni0;         
    *inneriter       = (ni1-ni0)*natoms/2;         
}

#undef SIMD_WIDTH
#undef UNROLLI   
#undef UNROLLJ   

