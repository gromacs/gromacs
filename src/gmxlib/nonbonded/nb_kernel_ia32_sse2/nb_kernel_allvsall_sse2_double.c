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

#include "nb_kernel_allvsall_sse2_double.h"
#include "gmx_sse2_double.h"


#include <xmmintrin.h>
#include <emmintrin.h>

#define SIMD_WIDTH 2
#define UNROLLI    2
#define UNROLLJ    2

typedef struct 
{
	double *     x_align;
	double *     y_align;
	double *     z_align;
	double *     q_align;
	double *     fx_align;
	double *     fy_align;
	double *     fz_align;	
    double **    pvdwaram_align;
    double **    ppvdw;
    int *        jindex;
    int *        imask;
    int **       prologue_mask;
    int **       epilogue_mask;
    int          ntype;
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
setup_exclusions_and_indices_double(gmx_allvsall_data_t *   aadata,
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
    snew(pi,2*(natoms+UNROLLI+2*SIMD_WIDTH));
    aadata->imask = (int *) (((size_t) pi + 16) & (~((size_t) 15)));
    for(i=0;i<natoms+UNROLLI;i++)
    {
        aadata->imask[2*i]   = (i>=start && i<end) ? 0xFFFFFFFF : 0;
        aadata->imask[2*i+1] = (i>=start && i<end) ? 0xFFFFFFFF : 0;
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
            snew(pi,2*(nj+2*SIMD_WIDTH));
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
                    aadata->prologue_mask[i][2*k]   = 0xFFFFFFFF;
                    aadata->prologue_mask[i][2*k+1] = 0xFFFFFFFF;
                }
                else
                {
                    aadata->prologue_mask[i][2*k]   = 0;
                    aadata->prologue_mask[i][2*k+1] = 0;
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
                        aadata->prologue_mask[i][2*k]   = 0;
                        aadata->prologue_mask[i][2*k+1] = 0;
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
            snew(pi,2*(nj+2*SIMD_WIDTH));
            aadata->epilogue_mask[i] = (int *) (((size_t) pi + 16) & (~((size_t) 15)));
            
            max_offset = calc_maxoffset(i,natoms);
            
            for(k=0;k<nj;k++)
            {
                j = aadata->jindex[4*i+2] + k;
                aadata->epilogue_mask[i][2*k]   = (j <= i+max_offset) ? 0xFFFFFFFF : 0;
                aadata->epilogue_mask[i][2*k+1] = (j <= i+max_offset) ? 0xFFFFFFFF : 0;
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
             double *                  pvdwaram)
{
	int i,j,k,idx;
	gmx_allvsall_data_t *aadata;
	double *pr;
    double *p;
	double *c6tmp,*c12tmp;
        
	snew(aadata,1);
	*p_aadata = aadata;
    
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->x_align = (double *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->y_align = (double *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->z_align = (double *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->fx_align = (double *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->fy_align = (double *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->fz_align = (double *) (((size_t) pr + 16) & (~((size_t) 15)));
	snew(pr,2*natoms+2*SIMD_WIDTH);
	aadata->q_align = (double *) (((size_t) pr + 16) & (~((size_t) 15)));	
    
    for(i=0;i<2*natoms+SIMD_WIDTH;i++)
	{
		aadata->x_align[i] = 0.0;
		aadata->y_align[i] = 0.0;
		aadata->z_align[i] = 0.0;
		aadata->q_align[i] = 0.0;
		aadata->fx_align[i] = 0.0;
		aadata->fy_align[i] = 0.0;
		aadata->fz_align[i] = 0.0;
	}
    
    /* Generate vdw params */
    snew(aadata->pvdwaram_align,ntype);
    snew(c6tmp,2*natoms+SIMD_WIDTH);
    snew(c12tmp,2*natoms+SIMD_WIDTH);
        
    for(i=0;i<ntype;i++)
    {
        /* Note that this 4 does NOT refer to SIMD_WIDTH, but to c6 & c12 params for 2*natoms! */
    	snew(pr,4*natoms+4*SIMD_WIDTH);
        aadata->pvdwaram_align[i] = (double *) (((size_t) pr + 16) & (~((size_t) 15)));
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
    
    setup_exclusions_and_indices_double(aadata,excl,start,end,natoms);
}




void
nb_kernel_allvsall_sse2_double(t_forcerec *           fr,
                               t_mdatoms *            mdatoms,
                               t_blocka *             excl,    
                               double *               x,
                               double *               f,
                               double *               Vc,
                               double *               Vvdw,
                               int *                  outeriter,
                               int *                  inneriter,
                               void *                 work)
{
	int        natoms;
	int        ni0,ni1;
	int        nj0,nj1,nj2,nj3;
	int        i,j,k;
	double *   charge;
	int *      type;
    double     facel;
	double **  pvdwaram_align;
    double **  ppvdw;
    double *   pvdw0;
    double *   pvdw1;
	int        ggid;
	gmx_allvsall_data_t *aadata;
	double *   x_align;
	double *   y_align;
	double *   z_align;
	double *   fx_align;
	double *   fy_align;
	double *   fz_align;
	double *   q_align;
	int        maxoffset;
    int *      jindex;
    int **     prologue_mask;
    int **     epilogue_mask;
    int *      pmask0;
    int *      pmask1;
    int *      emask0;
    int *      emask1;
    int *      imask;
    double     tmpsum[2];
    
    __m128d    ix_SSE0,iy_SSE0,iz_SSE0;
    __m128d    ix_SSE1,iy_SSE1,iz_SSE1;
	__m128d    fix_SSE0,fiy_SSE0,fiz_SSE0;
	__m128d    fix_SSE1,fiy_SSE1,fiz_SSE1;
	__m128d    fjxSSE,fjySSE,fjzSSE;
	__m128d    jxSSE,jySSE,jzSSE,jqSSE;
	__m128d    dx_SSE0,dy_SSE0,dz_SSE0;
	__m128d    dx_SSE1,dy_SSE1,dz_SSE1;
	__m128d    tx_SSE0,ty_SSE0,tz_SSE0;
	__m128d    tx_SSE1,ty_SSE1,tz_SSE1;
	__m128d    rsq_SSE0,rinv_SSE0,rinvsq_SSE0,rinvsix_SSE0;
	__m128d    rsq_SSE1,rinv_SSE1,rinvsq_SSE1,rinvsix_SSE1;
	__m128d    qq_SSE0,iq_SSE0;
	__m128d    qq_SSE1,iq_SSE1;
	__m128d    vcoul_SSE0,Vvdw6_SSE0,Vvdw12_SSE0,fscal_SSE0;
	__m128d    vcoul_SSE1,Vvdw6_SSE1,Vvdw12_SSE1,fscal_SSE1;
    __m128d    c6_SSE0,c12_SSE0;
    __m128d    c6_SSE1,c12_SSE1;
    
	__m128d    vctotSSE,VvdwtotSSE;
	__m128d    sixSSE,twelveSSE;
	__m128d    imask_SSE0,imask_SSE1;
    __m128d    jmask_SSE0,jmask_SSE1;
    __m128d    tmpSSE;
    
	charge              = mdatoms->chargeA;
	type                = mdatoms->typeA;
	facel               = fr->epsfac;
    
	natoms              = mdatoms->nr;
	ni0                 = (mdatoms->start/SIMD_WIDTH)*SIMD_WIDTH;
	ni1                 = mdatoms->start+mdatoms->homenr;
    
    sixSSE    = _mm_set1_pd(6.0);
	twelveSSE = _mm_set1_pd(12.0);
    
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
    
    prologue_mask = aadata->prologue_mask;
    epilogue_mask = aadata->epilogue_mask;
    jindex        = aadata->jindex;
    imask         = aadata->imask;
    
	for(i=ni0;i<ni1+1+natoms/2;i++)
    {
        k = i%natoms;
		x_align[i]  = x[3*k];
		y_align[i]  = x[3*k+1];
		z_align[i]  = x[3*k+2];
		q_align[i]  = charge[k];
		fx_align[i] = 0;
		fy_align[i] = 0;
		fz_align[i] = 0;
	}
    
	for(i=ni0; i<ni1; i+=UNROLLI)
	{
		/* We assume shifts are NOT used for all-vs-all interactions */
		
		/* Load i atom data */
		ix_SSE0          = _mm_load1_pd(x_align+i);
		ix_SSE1          = _mm_load1_pd(x_align+i+1);
		iy_SSE0          = _mm_load1_pd(y_align+i);
		iy_SSE1          = _mm_load1_pd(y_align+i+1);
		iz_SSE0          = _mm_load1_pd(z_align+i);
		iz_SSE1          = _mm_load1_pd(z_align+i+1);
		iq_SSE0          = _mm_set1_pd(facel*q_align[i]);
		iq_SSE1          = _mm_set1_pd(facel*q_align[i+1]);

        pvdw0            = ppvdw[i];
        pvdw1            = ppvdw[i+1];

		/* Zero the potential energy for this list */
		VvdwtotSSE       = _mm_setzero_pd();
		vctotSSE         = _mm_setzero_pd();

		/* Clear i atom forces */
		fix_SSE0           = _mm_setzero_pd();
		fix_SSE1           = _mm_setzero_pd();
		fiy_SSE0           = _mm_setzero_pd();
		fiy_SSE1           = _mm_setzero_pd();
		fiz_SSE0           = _mm_setzero_pd();
		fiz_SSE1           = _mm_setzero_pd();
        
		/* Load limits for loop over neighbors */
		nj0              = jindex[4*i];
		nj1              = jindex[4*i+1];
        nj2              = jindex[4*i+2];
        nj3              = jindex[4*i+3];

        pmask0           = prologue_mask[i];
        pmask1           = prologue_mask[i+1];
        emask0           = epilogue_mask[i]; 
        emask1           = epilogue_mask[i+1];
        imask_SSE0       = _mm_load1_pd((double *)(imask+2*i));
        imask_SSE1       = _mm_load1_pd((double *)(imask+2*i+2));
            
        printf("OUTER i=%d %d     imask=%d/%d %d/%d\n",i,i+1,*(imask+2*i),*(imask+2*i+1),*(imask+2*i+2),*(imask+2*i+3));
        
        for(j=nj0; j<nj1; j+=UNROLLJ)
        {                        
            jmask_SSE0 = _mm_load_pd((double *)pmask0);
            jmask_SSE1 = _mm_load_pd((double *)pmask1);
            printf("PROLOGUE j=%d %d    mask=%d/%d %d/%d\n",j,j+1,*(pmask0),*(pmask0+1),*(pmask1),*(pmask1+1));

            pmask0 += 2*UNROLLJ;
            pmask1 += 2*UNROLLJ;
            
            
            /* load j atom coordinates */
            jxSSE            = _mm_load_pd(x_align+j);
            jySSE            = _mm_load_pd(y_align+j);
            jzSSE            = _mm_load_pd(z_align+j);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_pd(ix_SSE0,jxSSE);
            dy_SSE0            = _mm_sub_pd(iy_SSE0,jySSE);
            dz_SSE0            = _mm_sub_pd(iz_SSE0,jzSSE);
            dx_SSE1            = _mm_sub_pd(ix_SSE1,jxSSE);
            dy_SSE1            = _mm_sub_pd(iy_SSE1,jySSE);
            dz_SSE1            = _mm_sub_pd(iz_SSE1,jzSSE);
            
            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_pd(dx_SSE0,dy_SSE0,dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_pd(dx_SSE1,dy_SSE1,dz_SSE1);

            /* Combine masks */
            jmask_SSE0         = _mm_and_pd(jmask_SSE0,imask_SSE0);
            jmask_SSE1         = _mm_and_pd(jmask_SSE1,imask_SSE1);
            
            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_pd(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_pd(rsq_SSE1);
            
            /* Apply mask */
            rinv_SSE0          = _mm_and_pd(rinv_SSE0,jmask_SSE0);
            rinv_SSE1          = _mm_and_pd(rinv_SSE1,jmask_SSE1);
            
            /* Load parameters for j atom */
            jqSSE             = _mm_load_pd(q_align+j);
            qq_SSE0            = _mm_mul_pd(iq_SSE0,jqSSE);
            qq_SSE1            = _mm_mul_pd(iq_SSE1,jqSSE);

            c6_SSE0            = _mm_load_pd(pvdw0+2*j);
            c6_SSE1            = _mm_load_pd(pvdw1+2*j);
            c12_SSE0           = _mm_load_pd(pvdw0+2*j+UNROLLJ);
            c12_SSE1           = _mm_load_pd(pvdw1+2*j+UNROLLJ);
            
            rinvsq_SSE0        = _mm_mul_pd(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1        = _mm_mul_pd(rinv_SSE1,rinv_SSE1);
            
            /* Coulomb interaction */
            vcoul_SSE0         = _mm_mul_pd(qq_SSE0,rinv_SSE0);
            vcoul_SSE1         = _mm_mul_pd(qq_SSE1,rinv_SSE1);

            vctotSSE           = _mm_add_pd(vctotSSE, _mm_add_pd(vcoul_SSE0,vcoul_SSE1));
            
            /* Lennard-Jones interaction */
            rinvsix_SSE0       = _mm_mul_pd(rinvsq_SSE0,_mm_mul_pd(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1       = _mm_mul_pd(rinvsq_SSE1,_mm_mul_pd(rinvsq_SSE1,rinvsq_SSE1));
            Vvdw6_SSE0         = _mm_mul_pd(c6_SSE0,rinvsix_SSE0);
            Vvdw6_SSE1         = _mm_mul_pd(c6_SSE1,rinvsix_SSE1);
            Vvdw12_SSE0        = _mm_mul_pd(c12_SSE0,_mm_mul_pd(rinvsix_SSE0,rinvsix_SSE0));
            Vvdw12_SSE1        = _mm_mul_pd(c12_SSE1,_mm_mul_pd(rinvsix_SSE1,rinvsix_SSE1));

            VvdwtotSSE         = _mm_add_pd(VvdwtotSSE, _mm_add_pd(_mm_sub_pd(Vvdw12_SSE0,Vvdw6_SSE0),
                                                                   _mm_sub_pd(Vvdw12_SSE1,Vvdw6_SSE1)));
            
            fscal_SSE0         = _mm_mul_pd(rinvsq_SSE0, 
                                          _mm_add_pd(vcoul_SSE0,
                                                     _mm_sub_pd(_mm_mul_pd(twelveSSE,Vvdw12_SSE0),
                                                                _mm_mul_pd(sixSSE,Vvdw6_SSE0))));
            fscal_SSE1         = _mm_mul_pd(rinvsq_SSE1, 
                                          _mm_add_pd(vcoul_SSE1,
                                                     _mm_sub_pd(_mm_mul_pd(twelveSSE,Vvdw12_SSE1),
                                                                _mm_mul_pd(sixSSE,Vvdw6_SSE1))));
          
            /* Calculate temporary vectorial force */
            tx_SSE0            = _mm_mul_pd(fscal_SSE0,dx_SSE0);
            tx_SSE1            = _mm_mul_pd(fscal_SSE1,dx_SSE1);
            ty_SSE0            = _mm_mul_pd(fscal_SSE0,dy_SSE0);
            ty_SSE1            = _mm_mul_pd(fscal_SSE1,dy_SSE1);
            tz_SSE0            = _mm_mul_pd(fscal_SSE0,dz_SSE0);
            tz_SSE1            = _mm_mul_pd(fscal_SSE1,dz_SSE1);
            
            /* Increment i atom force */
            fix_SSE0          = _mm_add_pd(fix_SSE0,tx_SSE0);
            fix_SSE1          = _mm_add_pd(fix_SSE1,tx_SSE1);
            fiy_SSE0          = _mm_add_pd(fiy_SSE0,ty_SSE0);
            fiy_SSE1          = _mm_add_pd(fiy_SSE1,ty_SSE1);
            fiz_SSE0          = _mm_add_pd(fiz_SSE0,tz_SSE0);
            fiz_SSE1          = _mm_add_pd(fiz_SSE1,tz_SSE1);
            
            /* Decrement j atom force */
            _mm_store_pd(fx_align+j,
                         _mm_sub_pd( _mm_load_pd(fx_align+j) , _mm_add_pd(tx_SSE0,tx_SSE1) ));
            _mm_store_pd(fy_align+j,
                         _mm_sub_pd( _mm_load_pd(fy_align+j) , _mm_add_pd(ty_SSE0,ty_SSE1) ));
            _mm_store_pd(fz_align+j,
                         _mm_sub_pd( _mm_load_pd(fz_align+j) , _mm_add_pd(tz_SSE0,tz_SSE1) ));
            
            
            /* Inner loop uses 38 flops/iteration */
        }

        for(j=nj1; j<nj2; j+=UNROLLJ)
        {                      
            printf("MAIN j=%d %d\n",j,j+1);

            /* load j atom coordinates */
            jxSSE            = _mm_load_pd(x_align+j);
            jySSE            = _mm_load_pd(y_align+j);
            jzSSE            = _mm_load_pd(z_align+j);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_pd(ix_SSE0,jxSSE);
            dy_SSE0            = _mm_sub_pd(iy_SSE0,jySSE);
            dz_SSE0            = _mm_sub_pd(iz_SSE0,jzSSE);
            dx_SSE1            = _mm_sub_pd(ix_SSE1,jxSSE);
            dy_SSE1            = _mm_sub_pd(iy_SSE1,jySSE);
            dz_SSE1            = _mm_sub_pd(iz_SSE1,jzSSE);
           
            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_pd(dx_SSE0,dy_SSE0,dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_pd(dx_SSE1,dy_SSE1,dz_SSE1);
            
            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_pd(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_pd(rsq_SSE1);
            
            rinv_SSE0          = _mm_and_pd(rinv_SSE0,imask_SSE0);
            rinv_SSE1          = _mm_and_pd(rinv_SSE1,imask_SSE1);

            /* Load parameters for j atom */
            jqSSE             = _mm_load_pd(q_align+j);
            qq_SSE0            = _mm_mul_pd(iq_SSE0,jqSSE);
            qq_SSE1            = _mm_mul_pd(iq_SSE1,jqSSE);
            
            c6_SSE0            = _mm_load_pd(pvdw0+2*j);
            c6_SSE1            = _mm_load_pd(pvdw1+2*j);
            c12_SSE0           = _mm_load_pd(pvdw0+2*j+UNROLLJ);
            c12_SSE1           = _mm_load_pd(pvdw1+2*j+UNROLLJ);
            
            rinvsq_SSE0        = _mm_mul_pd(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1        = _mm_mul_pd(rinv_SSE1,rinv_SSE1);
            
            /* Coulomb interaction */
            vcoul_SSE0         = _mm_mul_pd(qq_SSE0,rinv_SSE0);
            vcoul_SSE1         = _mm_mul_pd(qq_SSE1,rinv_SSE1);
            
            /* Lennard-Jones interaction */
            rinvsix_SSE0       = _mm_mul_pd(rinvsq_SSE0,_mm_mul_pd(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1       = _mm_mul_pd(rinvsq_SSE1,_mm_mul_pd(rinvsq_SSE1,rinvsq_SSE1));
            Vvdw6_SSE0         = _mm_mul_pd(c6_SSE0,rinvsix_SSE0);
            Vvdw6_SSE1         = _mm_mul_pd(c6_SSE1,rinvsix_SSE1);
            Vvdw12_SSE0        = _mm_mul_pd(c12_SSE0,_mm_mul_pd(rinvsix_SSE0,rinvsix_SSE0));
            Vvdw12_SSE1        = _mm_mul_pd(c12_SSE1,_mm_mul_pd(rinvsix_SSE1,rinvsix_SSE1));
            
            vctotSSE           = _mm_add_pd(vctotSSE, _mm_add_pd(vcoul_SSE0,vcoul_SSE1));
            VvdwtotSSE         = _mm_add_pd(VvdwtotSSE, _mm_add_pd(_mm_sub_pd(Vvdw12_SSE0,Vvdw6_SSE0),
                                                                   _mm_sub_pd(Vvdw12_SSE1,Vvdw6_SSE1)));
                                                                            
            fscal_SSE0         = _mm_mul_pd(rinvsq_SSE0, 
                                           _mm_add_pd(vcoul_SSE0,
                                                      _mm_sub_pd(_mm_mul_pd(twelveSSE,Vvdw12_SSE0),
                                                                 _mm_mul_pd(sixSSE,Vvdw6_SSE0))));
            fscal_SSE1         = _mm_mul_pd(rinvsq_SSE1, 
                                           _mm_add_pd(vcoul_SSE1,
                                                      _mm_sub_pd(_mm_mul_pd(twelveSSE,Vvdw12_SSE1),
                                                                 _mm_mul_pd(sixSSE,Vvdw6_SSE1))));

            /* Calculate temporary vectorial force */
            tx_SSE0            = _mm_mul_pd(fscal_SSE0,dx_SSE0);
            tx_SSE1            = _mm_mul_pd(fscal_SSE1,dx_SSE1);
            ty_SSE0            = _mm_mul_pd(fscal_SSE0,dy_SSE0);
            ty_SSE1            = _mm_mul_pd(fscal_SSE1,dy_SSE1);
            tz_SSE0            = _mm_mul_pd(fscal_SSE0,dz_SSE0);
            tz_SSE1            = _mm_mul_pd(fscal_SSE1,dz_SSE1);
            
            /* Increment i atom force */
            fix_SSE0          = _mm_add_pd(fix_SSE0,tx_SSE0);
            fix_SSE1          = _mm_add_pd(fix_SSE1,tx_SSE1);
            fiy_SSE0          = _mm_add_pd(fiy_SSE0,ty_SSE0);
            fiy_SSE1          = _mm_add_pd(fiy_SSE1,ty_SSE1);
            fiz_SSE0          = _mm_add_pd(fiz_SSE0,tz_SSE0);
            fiz_SSE1          = _mm_add_pd(fiz_SSE1,tz_SSE1);
            
            /* Decrement j atom force */
            _mm_store_pd(fx_align+j,
                         _mm_sub_pd( _mm_load_pd(fx_align+j) , _mm_add_pd(tx_SSE0,tx_SSE1) ));
            _mm_store_pd(fy_align+j,
                         _mm_sub_pd( _mm_load_pd(fy_align+j) , _mm_add_pd(ty_SSE0,ty_SSE1) ));
            _mm_store_pd(fz_align+j,
                         _mm_sub_pd( _mm_load_pd(fz_align+j) , _mm_add_pd(tz_SSE0,tz_SSE1) ));
            
            /* Inner loop uses 38 flops/iteration */
        }

        for(j=nj2; j<nj3; j+=UNROLLJ)
        {
            jmask_SSE0 = _mm_load_pd((double *)emask0);
            jmask_SSE1 = _mm_load_pd((double *)emask1);
            printf("EPILOGUE j=%d %d    mask=%d/%d %d/%d\n",j,j+1,*(emask0),*(emask0+1),*(emask1),*(emask1+1));

            emask0 += 2*UNROLLJ;
            emask1 += 2*UNROLLJ;
            
            /* load j atom coordinates */
            jxSSE            = _mm_load_pd(x_align+j);
            jySSE            = _mm_load_pd(y_align+j);
            jzSSE            = _mm_load_pd(z_align+j);
            
            /* Calculate distance */
            dx_SSE0            = _mm_sub_pd(ix_SSE0,jxSSE);
            dy_SSE0            = _mm_sub_pd(iy_SSE0,jySSE);
            dz_SSE0            = _mm_sub_pd(iz_SSE0,jzSSE);
            dx_SSE1            = _mm_sub_pd(ix_SSE1,jxSSE);
            dy_SSE1            = _mm_sub_pd(iy_SSE1,jySSE);
            dz_SSE1            = _mm_sub_pd(iz_SSE1,jzSSE);
            
            /* rsq = dx*dx+dy*dy+dz*dz */
            rsq_SSE0           = gmx_mm_calc_rsq_pd(dx_SSE0,dy_SSE0,dz_SSE0);
            rsq_SSE1           = gmx_mm_calc_rsq_pd(dx_SSE1,dy_SSE1,dz_SSE1);
            
            jmask_SSE0         = _mm_and_pd(jmask_SSE0,imask_SSE0);
            jmask_SSE1         = _mm_and_pd(jmask_SSE1,imask_SSE1);
            
            /* Calculate 1/r and 1/r2 */
            rinv_SSE0          = gmx_mm_invsqrt_pd(rsq_SSE0);
            rinv_SSE1          = gmx_mm_invsqrt_pd(rsq_SSE1);
            rinv_SSE0          = _mm_and_pd(rinv_SSE0,jmask_SSE0);
            rinv_SSE1          = _mm_and_pd(rinv_SSE1,jmask_SSE1);
            
            /* Load parameters for j atom */
            jqSSE             = _mm_load_pd(q_align+j);
            qq_SSE0            = _mm_mul_pd(iq_SSE0,jqSSE);
            qq_SSE1            = _mm_mul_pd(iq_SSE1,jqSSE);
            
            c6_SSE0            = _mm_load_pd(pvdw0+2*j);
            c6_SSE1            = _mm_load_pd(pvdw1+2*j);
            c12_SSE0           = _mm_load_pd(pvdw0+2*j+UNROLLJ);
            c12_SSE1           = _mm_load_pd(pvdw1+2*j+UNROLLJ);
            
            rinvsq_SSE0        = _mm_mul_pd(rinv_SSE0,rinv_SSE0);
            rinvsq_SSE1        = _mm_mul_pd(rinv_SSE1,rinv_SSE1);
            
            /* Coulomb interaction */
            vcoul_SSE0         = _mm_mul_pd(qq_SSE0,rinv_SSE0);
            vcoul_SSE1         = _mm_mul_pd(qq_SSE1,rinv_SSE1);

            vctotSSE           = _mm_add_pd(vctotSSE, _mm_add_pd(vcoul_SSE0,vcoul_SSE1));
            
            /* Lennard-Jones interaction */
            rinvsix_SSE0       = _mm_mul_pd(rinvsq_SSE0,_mm_mul_pd(rinvsq_SSE0,rinvsq_SSE0));
            rinvsix_SSE1       = _mm_mul_pd(rinvsq_SSE1,_mm_mul_pd(rinvsq_SSE1,rinvsq_SSE1));
            Vvdw6_SSE0         = _mm_mul_pd(c6_SSE0,rinvsix_SSE0);
            Vvdw6_SSE1         = _mm_mul_pd(c6_SSE1,rinvsix_SSE1);
            Vvdw12_SSE0        = _mm_mul_pd(c12_SSE0,_mm_mul_pd(rinvsix_SSE0,rinvsix_SSE0));
            Vvdw12_SSE1        = _mm_mul_pd(c12_SSE1,_mm_mul_pd(rinvsix_SSE1,rinvsix_SSE1));
            
            VvdwtotSSE         = _mm_add_pd(VvdwtotSSE, _mm_add_pd(_mm_sub_pd(Vvdw12_SSE0,Vvdw6_SSE0),
                                                                   _mm_sub_pd(Vvdw12_SSE1,Vvdw6_SSE1)));
            
            fscal_SSE0         = _mm_mul_pd(rinvsq_SSE0, 
                                           _mm_add_pd(vcoul_SSE0,
                                                      _mm_sub_pd(_mm_mul_pd(twelveSSE,Vvdw12_SSE0),
                                                                 _mm_mul_pd(sixSSE,Vvdw6_SSE0))));
            fscal_SSE1         = _mm_mul_pd(rinvsq_SSE1, 
                                           _mm_add_pd(vcoul_SSE1,
                                                      _mm_sub_pd(_mm_mul_pd(twelveSSE,Vvdw12_SSE1),
                                                                 _mm_mul_pd(sixSSE,Vvdw6_SSE1))));
            /* Calculate temporary vectorial force */
            tx_SSE0            = _mm_mul_pd(fscal_SSE0,dx_SSE0);
            tx_SSE1            = _mm_mul_pd(fscal_SSE1,dx_SSE1);
            ty_SSE0            = _mm_mul_pd(fscal_SSE0,dy_SSE0);
            ty_SSE1            = _mm_mul_pd(fscal_SSE1,dy_SSE1);
            tz_SSE0            = _mm_mul_pd(fscal_SSE0,dz_SSE0);
            tz_SSE1            = _mm_mul_pd(fscal_SSE1,dz_SSE1);
            
            /* Increment i atom force */
            fix_SSE0          = _mm_add_pd(fix_SSE0,tx_SSE0);
            fix_SSE1          = _mm_add_pd(fix_SSE1,tx_SSE1);
            fiy_SSE0          = _mm_add_pd(fiy_SSE0,ty_SSE0);
            fiy_SSE1          = _mm_add_pd(fiy_SSE1,ty_SSE1);
            fiz_SSE0          = _mm_add_pd(fiz_SSE0,tz_SSE0);
            fiz_SSE1          = _mm_add_pd(fiz_SSE1,tz_SSE1);
            
            /* Decrement j atom force */
            _mm_store_pd(fx_align+j,
                         _mm_sub_pd( _mm_load_pd(fx_align+j) , _mm_add_pd(tx_SSE0,tx_SSE1) ));
            _mm_store_pd(fy_align+j,
                         _mm_sub_pd( _mm_load_pd(fy_align+j) , _mm_add_pd(ty_SSE0,ty_SSE1) ));
            _mm_store_pd(fz_align+j,
                         _mm_sub_pd( _mm_load_pd(fz_align+j) , _mm_add_pd(tz_SSE0,tz_SSE1) ));
            /* Inner loop uses 38 flops/iteration */
        }
        
		/* Add i forces to mem and shifted force list */
        tmpSSE   = fix_SSE0;
        fix_SSE0 = _mm_unpacklo_pd(fix_SSE0,fix_SSE1);
        fix_SSE1 = _mm_unpackhi_pd(tmpSSE,fix_SSE1);
        fix_SSE0 = _mm_add_pd(fix_SSE0,fix_SSE1);
        _mm_store_pd(fx_align+i, _mm_add_pd(fix_SSE0, _mm_load_pd(fx_align+i)));
        
        tmpSSE   = fiy_SSE0;
        fiy_SSE0 = _mm_unpacklo_pd(fiy_SSE0,fiy_SSE1);
        fiy_SSE1 = _mm_unpackhi_pd(tmpSSE,fiy_SSE1);
        fiy_SSE0 = _mm_add_pd(fiy_SSE0,fiy_SSE1);
        _mm_store_pd(fy_align+i, _mm_add_pd(fiy_SSE0, _mm_load_pd(fy_align+i)));
        
        tmpSSE   = fiy_SSE0;
        fiz_SSE0 = _mm_unpacklo_pd(fiz_SSE0,fiz_SSE1);
        fiz_SSE1 = _mm_unpackhi_pd(tmpSSE,fiz_SSE1);
        fiz_SSE0 = _mm_add_pd(fiz_SSE0,fiz_SSE1);
        _mm_store_pd(fz_align+i, _mm_add_pd(fiz_SSE0, _mm_load_pd(fz_align+i)));
		
		/* Add potential energies to the group for this list */
		ggid             = 0;         
        
		_mm_storeu_pd(tmpsum,vctotSSE);
		Vc[ggid]         = Vc[ggid] + tmpsum[0]+tmpsum[1];
		
		_mm_storeu_pd(tmpsum,VvdwtotSSE);
		Vvdw[ggid]       = Vvdw[ggid] + tmpsum[0]+tmpsum[1];
		
		/* Outer loop uses 6 flops/iteration */
	}    
    
	for(i=ni0;i<ni1+1+natoms/2;i++)
	{
        k = i%natoms;
		f[3*k]   += fx_align[i];
		f[3*k+1] += fy_align[i];
		f[3*k+2] += fz_align[i];
	}
    
    
    /* Write outer/inner iteration count to pointers */
    *outeriter       = ni1-ni0;         
    *inneriter       = (ni1-ni0)*natoms/2;         
}

#undef SIMD_WIDTH
#undef UNROLLI   
#undef UNROLLJ   

