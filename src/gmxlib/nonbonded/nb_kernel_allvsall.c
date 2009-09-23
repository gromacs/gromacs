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

#include "nb_kernel_allvsall.h"


#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>

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
    real **    pvdwaram_align;
    real **    ppvdw;
    int *      jindex;
    int *      imask;
    int **     prologue_mask;
    int **     epilogue_mask;
    int        ntype;
} 
gmx_allvsall_data_t;
		

void
setup_exclusions_and_indices_float(gmx_allvsall_data_t *   aadata,
                                   t_blocka *              excl,  
                                   int                     start,
                                   int                     end,
                                   int                     natoms)
{
    int i,j,k;
    int ni0,ni1,nj0,nj1,nj;
    int imin,imax,maxoffset;
    int firstinteraction[UNROLLI];
    int ibase;
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
	
    /* First zero everything to avoid uninitialized data */
    for(i=0;i<natoms+UNROLLI;i++)
    {
        aadata->prologue_mask[i] = NULL;
    }
    
    /* Calculate the largest exclusion range we need for each UNROLLI-tuplet of i atoms. */
    for(ibase=ni0;ibase<ni1;ibase+=UNROLLI)
	{
        maxoffset = -1;
        
        /* First find maxoffset for the next 4 atoms (or fewer if we are close to end) */
        imax = ((ibase+UNROLLI) < end) ? (ibase+UNROLLI) : end;
        
        for(i=ibase;i<imax;i++)
        {
            /* Before exclusions, which atom is the first we (might) interact with? */
            firstinteraction[i-ibase] = i+1;
            
            nj0   = excl->index[i];
            nj1   = excl->index[i+1];
            for(j=nj0; j<nj1; j++)
            {
                k = excl->a[j] - ibase;
                if(k>maxoffset)
                {
                    maxoffset=k;                    
                }
                /* Exclusions are sorted, so this can be done iteratively */
                if(excl->a[j] == firstinteraction[i-ibase])
                {
                    firstinteraction[i-ibase]++;
                }
            }
        }
        /* round up to j unrolling factor */
        maxoffset = (maxoffset/UNROLLJ+1)*UNROLLJ;

        imin = firstinteraction[0];
        for(i=ibase;i<imax;i++)
        {
            imin = (imin < firstinteraction[i-ibase]) ? imin : firstinteraction[i-ibase];
        }
        imin = (imin/UNROLLJ)*UNROLLJ;
        
        /* Set all the prologue masks length to this value (even for i>end) */
        for(i=ibase;i<ibase+UNROLLI;i++)
        {
            aadata->jindex[4*i]   = imin;
            aadata->jindex[4*i+1] = ibase+maxoffset;
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
         
            maxoffset = natoms/2;
            if(natoms%2==0 && i>=natoms/2)
            {
                maxoffset--;
            }

            /* Include interactions i+1 <= j < i+maxoffset */
            for(k=0;k<nj;k++)
            {
                j = imin + k;
                                
                if( (j>i) && (j<=i+maxoffset) )
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
                    k = excl->a[j] - imin;
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
        maxoffset = natoms/2;
        if(natoms%2==0 && imin>=natoms/2)
        {
            maxoffset--;
        }
        imin = imin + 1 + maxoffset;
        
        /* Find largest index for which we need to use the epilogue */
        imax = ibase + UNROLLI-1;
        imax = (imax < end) ? imax : end; 
        
        maxoffset = natoms/2;
        if(natoms%2==0 && imax>=natoms/2)
        {
            maxoffset--;
        }
        imax = imax + 1 + maxoffset + UNROLLJ - 1;
        
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
            
            maxoffset = natoms/2;
            if(natoms%2==0 && i>=natoms/2)
            {
                maxoffset--;
            }
            
            for(k=0;k<nj;k++)
            {
                j = aadata->jindex[4*i+2] + k;
                aadata->epilogue_mask[i][k] = (j <= i+maxoffset) ? 0xFFFFFFFF : 0;
            }
        }
    }
}

void
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


static inline void
accumulate_SSE_force(__m128 fix1, __m128 fiy1, __m128 fiz1,float *fptr)
{
    __m128 t1,t2,t3;
    
    /* SSE2 */
    /* transpose data */
    t1 = fix1;
    _MM_TRANSPOSE4_PS(fix1,t1,fiy1,fiz1);
    fix1 = _mm_add_ps(_mm_add_ps(fix1,t1), _mm_add_ps(fiy1,fiz1));

    t2 = _mm_load_ss(fptr);
    t2 = _mm_loadh_pi(t2,(__m64 *)(fptr+1));
    
    t2 = _mm_add_ps(t2,fix1);
    
    _mm_store_ss(fptr,t2);
    _mm_storeh_pi((__m64 *)(fptr+1),t2);
}

static inline __m128
gmx_mm_invsqrt_ps(__m128 x)
{
    const __m128 half  = {0.5,0.5,0.5,0.5};
    const __m128 three = {3.0,3.0,3.0,3.0};
    
    __m128 lu = _mm_rsqrt_ps(x);
    
    return _mm_mul_ps(half,_mm_mul_ps(_mm_sub_ps(three,_mm_mul_ps(_mm_mul_ps(lu,lu),x)),lu));
}

void
printxmm(const char *s,__m128 xmm)
{
    float f[4];
    
    _mm_storeu_ps(f,xmm);
    printf("%s: %g %g %g %g\n",s,f[0],f[1],f[2],f[3]);      
}

void
printxmmi(const char *s,__m128i xmmi)
{
    int i[4];
    
    _mm_storeu_si128((__m128i *)i,xmmi);
    printf("%s: %d %d %d %d\n",s,i[0],i[1],i[2],i[3]);      
}



void
nb_kernel_allvsall(t_forcerec *           fr,
				   t_mdatoms *            mdatoms,
				   t_blocka *             excl,    
				   real *                 x,
				   real *                 f,
				   real *                 Vc,
				   real *                 Vvdw,
				   int *                  outeriter,
				   int *                  inneriter,
				   void *                 work)
{
	int        natoms;
	int        ni0,ni1;
	int        nj0,nj1,nj2,nj3;
	int        i,j,k;
	real *     charge;
	int *      type;
    real       facel;
	real **    pvdwaram_align;
    real **    ppvdw;
    real *     pvdw0;
    real *     pvdw1;
    real *     pvdw2;
    real *     pvdw3;
	int        ggid;
	gmx_allvsall_data_t *aadata;
	real *     x_align;
	real *     y_align;
	real *     z_align;
	real *     fx_align;
	real *     fy_align;
	real *     fz_align;
	real *     q_align;
	int        prologue0,prologue1,main0,main1,epilogue0,epilogue1;
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
    real       tmpsum[4];
    
    __m128     ix0SSE,iy0SSE,iz0SSE;
    __m128     ix1SSE,iy1SSE,iz1SSE;
    __m128     ix2SSE,iy2SSE,iz2SSE;
    __m128     ix3SSE,iy3SSE,iz3SSE;
	__m128     fix0SSE,fiy0SSE,fiz0SSE;
	__m128     fix1SSE,fiy1SSE,fiz1SSE;
	__m128     fix2SSE,fiy2SSE,fiz2SSE;
	__m128     fix3SSE,fiy3SSE,fiz3SSE;
	__m128     fjxSSE,fjySSE,fjzSSE;
	__m128     jxSSE,jySSE,jzSSE,jqSSE;
	__m128     dx0SSE,dy0SSE,dz0SSE;
	__m128     dx1SSE,dy1SSE,dz1SSE;
	__m128     dx2SSE,dy2SSE,dz2SSE;
	__m128     dx3SSE,dy3SSE,dz3SSE;
	__m128     tx0SSE,ty0SSE,tz0SSE;
	__m128     tx1SSE,ty1SSE,tz1SSE;
	__m128     tx2SSE,ty2SSE,tz2SSE;
	__m128     tx3SSE,ty3SSE,tz3SSE;
	__m128     rsq0SSE,rinv0SSE,rinvsq0SSE,rinvsix0SSE;
	__m128     rsq1SSE,rinv1SSE,rinvsq1SSE,rinvsix1SSE;
	__m128     rsq2SSE,rinv2SSE,rinvsq2SSE,rinvsix2SSE;
	__m128     rsq3SSE,rinv3SSE,rinvsq3SSE,rinvsix3SSE;
	__m128     tmpA0SSE,tmpB0SSE,qq0SSE,iq0SSE;
	__m128     tmpA1SSE,tmpB1SSE,qq1SSE,iq1SSE;
	__m128     tmpA2SSE,tmpB2SSE,qq2SSE,iq2SSE;
	__m128     tmpA3SSE,tmpB3SSE,qq3SSE,iq3SSE;
	__m128     vcoul0SSE,Vvdw60SSE,Vvdw120SSE,fscal0SSE;
	__m128     vcoul1SSE,Vvdw61SSE,Vvdw121SSE,fscal1SSE;
	__m128     vcoul2SSE,Vvdw62SSE,Vvdw122SSE,fscal2SSE;
	__m128     vcoul3SSE,Vvdw63SSE,Vvdw123SSE,fscal3SSE;
    __m128     c60SSE,c120SSE;
    __m128     c61SSE,c121SSE;
    __m128     c62SSE,c122SSE;
    __m128     c63SSE,c123SSE;
    
	__m128     vctotSSE,VvdwtotSSE;
	__m128     sixSSE,twelveSSE;
	__m128i    ikSSE,imSSE,ifourSSE;
	__m128i    ioneSSE;
	__m128     imask0SSE,imask1SSE,imask2SSE,imask3SSE;
    __m128     jmask0SSE,jmask1SSE,jmask2SSE,jmask3SSE;
    
	charge              = mdatoms->chargeA;
	type                = mdatoms->typeA;
	facel               = fr->epsfac;
    
	natoms            = mdatoms->nr;
	ni0               = mdatoms->start;
	ni1               = mdatoms->homenr;
    
    sixSSE    = _mm_set1_ps(6.0);
	twelveSSE = _mm_set1_ps(12.0);
	ifourSSE  = _mm_set1_epi32(4);
	ioneSSE   = _mm_set1_epi32(1);	
    
    aadata = *((gmx_allvsall_data_t **)work);

	if(aadata==NULL)
	{
		setup_aadata(&aadata,excl,ni0,ni1,natoms,type,fr->ntype,fr->nbfp);
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
    
	for(i=0;i<natoms;i++)
	{
		x_align[i] = x[3*i];
		y_align[i] = x[3*i+1];
		z_align[i] = x[3*i+2];
		q_align[i] = charge[i];
		fx_align[i] = 0;
		fy_align[i] = 0;
		fz_align[i] = 0;
	}
    
    /* Copy again */
	for(i=0;i<natoms;i++)
	{
		x_align[natoms+i] = x_align[i];
		y_align[natoms+i] = y_align[i];
		z_align[natoms+i] = z_align[i];
		q_align[natoms+i] = q_align[i];
		fx_align[natoms+i]   = 0;
		fy_align[natoms+i]   = 0;
		fz_align[natoms+i]   = 0;
    }        
    
	for(i=ni0; i<ni1; i+=UNROLLI)
	{
		/* We assume shifts are NOT used for all-vs-all interactions */
		
		/* Load i atom data */
		ix0SSE            = _mm_load1_ps(x_align+i);
		ix1SSE            = _mm_load1_ps(x_align+i+1);
		ix2SSE            = _mm_load1_ps(x_align+i+2);
		ix3SSE            = _mm_load1_ps(x_align+i+3);
		iy0SSE            = _mm_load1_ps(y_align+i);
		iy1SSE            = _mm_load1_ps(y_align+i+1);
		iy2SSE            = _mm_load1_ps(y_align+i+2);
		iy3SSE            = _mm_load1_ps(y_align+i+3);
		iz0SSE            = _mm_load1_ps(z_align+i);
		iz1SSE            = _mm_load1_ps(z_align+i+1);
		iz2SSE            = _mm_load1_ps(z_align+i+2);
		iz3SSE            = _mm_load1_ps(z_align+i+3);
		iq0SSE            = _mm_set1_ps(facel*q_align[i]);
		iq1SSE            = _mm_set1_ps(facel*q_align[i+1]);
		iq2SSE            = _mm_set1_ps(facel*q_align[i+2]);
		iq3SSE            = _mm_set1_ps(facel*q_align[i+3]);

        pvdw0            = ppvdw[i];
        pvdw1            = ppvdw[i+1];
        pvdw2            = ppvdw[i+2];
        pvdw3            = ppvdw[i+3];

		/* Zero the potential energy for this list */
		VvdwtotSSE       = _mm_setzero_ps();
		vctotSSE         = _mm_setzero_ps();

		/* Clear i atom forces */
		fix0SSE           = _mm_setzero_ps();
		fix1SSE           = _mm_setzero_ps();
		fix2SSE           = _mm_setzero_ps();
		fix3SSE           = _mm_setzero_ps();
		fiy0SSE           = _mm_setzero_ps();
		fiy1SSE           = _mm_setzero_ps();
		fiy2SSE           = _mm_setzero_ps();
		fiy3SSE           = _mm_setzero_ps();
		fiz0SSE           = _mm_setzero_ps();
		fiz1SSE           = _mm_setzero_ps();
		fiz2SSE           = _mm_setzero_ps();
		fiz3SSE           = _mm_setzero_ps();
        
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
        imask0SSE        = _mm_load1_ps((real *)(imask+i));
        imask1SSE        = _mm_load1_ps((real *)(imask+i+1));
        imask2SSE        = _mm_load1_ps((real *)(imask+i+2));
        imask3SSE        = _mm_load1_ps((real *)(imask+i+3));
        
        for(j=nj0; j<nj1; j+=UNROLLJ)
        {            
            jmask0SSE = _mm_load_ps((real *)pmask0);
            jmask1SSE = _mm_load_ps((real *)pmask1);
            jmask2SSE = _mm_load_ps((real *)pmask2);
            jmask3SSE = _mm_load_ps((real *)pmask3);
            pmask0 += UNROLLJ;
            pmask1 += UNROLLJ;
            pmask2 += UNROLLJ;
            pmask3 += UNROLLJ;
            
            /* load j atom coordinates */
            jxSSE            = _mm_load_ps(x_align+j);
            jySSE            = _mm_load_ps(y_align+j);
            jzSSE            = _mm_load_ps(z_align+j);
            
            /* Calculate distance */
            dx0SSE            = _mm_sub_ps(ix0SSE,jxSSE);
            dy0SSE            = _mm_sub_ps(iy0SSE,jySSE);
            dz0SSE            = _mm_sub_ps(iz0SSE,jzSSE);
            dx1SSE            = _mm_sub_ps(ix1SSE,jxSSE);
            dy1SSE            = _mm_sub_ps(iy1SSE,jySSE);
            dz1SSE            = _mm_sub_ps(iz1SSE,jzSSE);
            dx2SSE            = _mm_sub_ps(ix2SSE,jxSSE);
            dy2SSE            = _mm_sub_ps(iy2SSE,jySSE);
            dz2SSE            = _mm_sub_ps(iz2SSE,jzSSE);
            dx3SSE            = _mm_sub_ps(ix3SSE,jxSSE);
            dy3SSE            = _mm_sub_ps(iy3SSE,jySSE);
            dz3SSE            = _mm_sub_ps(iz3SSE,jzSSE);
            
            rsq0SSE           = _mm_mul_ps(dx0SSE,dx0SSE);
            rsq1SSE           = _mm_mul_ps(dx1SSE,dx1SSE);
            rsq2SSE           = _mm_mul_ps(dx2SSE,dx2SSE);
            rsq3SSE           = _mm_mul_ps(dx3SSE,dx3SSE);
            tmpA0SSE          = _mm_mul_ps(dy0SSE,dy0SSE);
            tmpA1SSE          = _mm_mul_ps(dy1SSE,dy1SSE);
            tmpA2SSE          = _mm_mul_ps(dy2SSE,dy2SSE);
            tmpA3SSE          = _mm_mul_ps(dy3SSE,dy3SSE);
            tmpB0SSE          = _mm_mul_ps(dz0SSE,dz0SSE);
            tmpB1SSE          = _mm_mul_ps(dz1SSE,dz1SSE);
            tmpB2SSE          = _mm_mul_ps(dz2SSE,dz2SSE);
            tmpB3SSE          = _mm_mul_ps(dz3SSE,dz3SSE);
            rsq0SSE           = _mm_add_ps(rsq0SSE,tmpA0SSE);
            rsq1SSE           = _mm_add_ps(rsq1SSE,tmpA1SSE);
            rsq2SSE           = _mm_add_ps(rsq2SSE,tmpA2SSE);
            rsq3SSE           = _mm_add_ps(rsq3SSE,tmpA3SSE);
            rsq0SSE           = _mm_add_ps(rsq0SSE,tmpB0SSE);
            rsq1SSE           = _mm_add_ps(rsq1SSE,tmpB1SSE);
            rsq2SSE           = _mm_add_ps(rsq2SSE,tmpB2SSE);
            rsq3SSE           = _mm_add_ps(rsq3SSE,tmpB3SSE);
            
            jmask0SSE         = _mm_and_ps(jmask0SSE,imask0SSE);
            jmask1SSE         = _mm_and_ps(jmask1SSE,imask1SSE);
            jmask2SSE         = _mm_and_ps(jmask2SSE,imask2SSE);
            jmask3SSE         = _mm_and_ps(jmask3SSE,imask3SSE);
            
            /* Calculate 1/r and 1/r2 */
            rinv0SSE          = gmx_mm_invsqrt_ps(rsq0SSE);
            rinv1SSE          = gmx_mm_invsqrt_ps(rsq1SSE);
            rinv2SSE          = gmx_mm_invsqrt_ps(rsq2SSE);
            rinv3SSE          = gmx_mm_invsqrt_ps(rsq3SSE);
            rinv0SSE          = _mm_and_ps(rinv0SSE,jmask0SSE);
            rinv1SSE          = _mm_and_ps(rinv1SSE,jmask1SSE);
            rinv2SSE          = _mm_and_ps(rinv2SSE,jmask2SSE);
            rinv3SSE          = _mm_and_ps(rinv3SSE,jmask3SSE);
            
            /* Load parameters for j atom */
            jqSSE             = _mm_load_ps(q_align+j);
            qq0SSE            = _mm_mul_ps(iq0SSE,jqSSE);
            qq1SSE            = _mm_mul_ps(iq1SSE,jqSSE);
            qq2SSE            = _mm_mul_ps(iq2SSE,jqSSE);
            qq3SSE            = _mm_mul_ps(iq3SSE,jqSSE);

            c60SSE            = _mm_load_ps(pvdw0+2*j);
            c61SSE            = _mm_load_ps(pvdw1+2*j);
            c62SSE            = _mm_load_ps(pvdw2+2*j);
            c63SSE            = _mm_load_ps(pvdw3+2*j);
            c120SSE           = _mm_load_ps(pvdw0+2*j+UNROLLJ);
            c121SSE           = _mm_load_ps(pvdw1+2*j+UNROLLJ);
            c122SSE           = _mm_load_ps(pvdw2+2*j+UNROLLJ);
            c123SSE           = _mm_load_ps(pvdw3+2*j+UNROLLJ);
                        
            rinvsq0SSE        = _mm_mul_ps(rinv0SSE,rinv0SSE);
            rinvsq1SSE        = _mm_mul_ps(rinv1SSE,rinv1SSE);
            rinvsq2SSE        = _mm_mul_ps(rinv2SSE,rinv2SSE);
            rinvsq3SSE        = _mm_mul_ps(rinv3SSE,rinv3SSE);
            
            /* Coulomb interaction */
            vcoul0SSE         = _mm_mul_ps(qq0SSE,rinv0SSE);
            vcoul1SSE         = _mm_mul_ps(qq1SSE,rinv1SSE);
            vcoul2SSE         = _mm_mul_ps(qq2SSE,rinv2SSE);
            vcoul3SSE         = _mm_mul_ps(qq3SSE,rinv3SSE);

            tmpA0SSE          = _mm_add_ps(vcoul0SSE,vcoul1SSE);
            tmpA2SSE          = _mm_add_ps(vcoul2SSE,vcoul3SSE);
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA2SSE);
            vctotSSE          = _mm_add_ps(vctotSSE,tmpA0SSE);
            
            /* Lennard-Jones interaction */
            rinvsix0SSE       = _mm_mul_ps(rinvsq0SSE,_mm_mul_ps(rinvsq0SSE,rinvsq0SSE));
            rinvsix1SSE       = _mm_mul_ps(rinvsq1SSE,_mm_mul_ps(rinvsq1SSE,rinvsq1SSE));
            rinvsix2SSE       = _mm_mul_ps(rinvsq2SSE,_mm_mul_ps(rinvsq2SSE,rinvsq2SSE));
            rinvsix3SSE       = _mm_mul_ps(rinvsq3SSE,_mm_mul_ps(rinvsq3SSE,rinvsq3SSE));
            Vvdw60SSE         = _mm_mul_ps(c60SSE,rinvsix0SSE);
            Vvdw61SSE         = _mm_mul_ps(c61SSE,rinvsix1SSE);
            Vvdw62SSE         = _mm_mul_ps(c62SSE,rinvsix2SSE);
            Vvdw63SSE         = _mm_mul_ps(c63SSE,rinvsix3SSE);
            Vvdw120SSE        = _mm_mul_ps(c120SSE,_mm_mul_ps(rinvsix0SSE,rinvsix0SSE));
            Vvdw121SSE        = _mm_mul_ps(c121SSE,_mm_mul_ps(rinvsix1SSE,rinvsix1SSE));
            Vvdw122SSE        = _mm_mul_ps(c122SSE,_mm_mul_ps(rinvsix2SSE,rinvsix2SSE));
            Vvdw123SSE        = _mm_mul_ps(c123SSE,_mm_mul_ps(rinvsix3SSE,rinvsix3SSE));

            tmpA0SSE          = _mm_sub_ps(Vvdw120SSE,Vvdw60SSE);
            tmpA1SSE          = _mm_sub_ps(Vvdw121SSE,Vvdw61SSE);
            tmpA2SSE          = _mm_sub_ps(Vvdw122SSE,Vvdw62SSE);
            tmpA3SSE          = _mm_sub_ps(Vvdw123SSE,Vvdw63SSE);
            
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA1SSE);
            tmpA2SSE          = _mm_add_ps(tmpA2SSE,tmpA3SSE);
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA2SSE);
            VvdwtotSSE       = _mm_add_ps(VvdwtotSSE, tmpA0SSE);
            
            fscal0SSE         = _mm_mul_ps(rinvsq0SSE, 
                                          _mm_add_ps(vcoul0SSE,
                                                     _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw120SSE),
                                                                _mm_mul_ps(sixSSE,Vvdw60SSE))));
            fscal1SSE         = _mm_mul_ps(rinvsq1SSE, 
                                          _mm_add_ps(vcoul1SSE,
                                                     _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw121SSE),
                                                                _mm_mul_ps(sixSSE,Vvdw61SSE))));
            fscal2SSE         = _mm_mul_ps(rinvsq2SSE, 
                                          _mm_add_ps(vcoul2SSE,
                                                     _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw122SSE),
                                                                _mm_mul_ps(sixSSE,Vvdw62SSE))));
            fscal3SSE         = _mm_mul_ps(rinvsq3SSE, 
                                          _mm_add_ps(vcoul3SSE,
                                                     _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw123SSE),
                                                                _mm_mul_ps(sixSSE,Vvdw63SSE))));
            
            /* Calculate temporary vectorial force */
            tx0SSE            = _mm_mul_ps(fscal0SSE,dx0SSE);
            tx1SSE            = _mm_mul_ps(fscal1SSE,dx1SSE);
            tx2SSE            = _mm_mul_ps(fscal2SSE,dx2SSE);
            tx3SSE            = _mm_mul_ps(fscal3SSE,dx3SSE);
            ty0SSE            = _mm_mul_ps(fscal0SSE,dy0SSE);
            ty1SSE            = _mm_mul_ps(fscal1SSE,dy1SSE);
            ty2SSE            = _mm_mul_ps(fscal2SSE,dy2SSE);
            ty3SSE            = _mm_mul_ps(fscal3SSE,dy3SSE);
            tz0SSE            = _mm_mul_ps(fscal0SSE,dz0SSE);
            tz1SSE            = _mm_mul_ps(fscal1SSE,dz1SSE);
            tz2SSE            = _mm_mul_ps(fscal2SSE,dz2SSE);
            tz3SSE            = _mm_mul_ps(fscal3SSE,dz3SSE);
            
            /* Increment i atom force */
            fix0SSE          = _mm_add_ps(fix0SSE,tx0SSE);
            fix1SSE          = _mm_add_ps(fix1SSE,tx1SSE);
            fix2SSE          = _mm_add_ps(fix2SSE,tx2SSE);
            fix3SSE          = _mm_add_ps(fix3SSE,tx3SSE);
            fiy0SSE          = _mm_add_ps(fiy0SSE,ty0SSE);
            fiy1SSE          = _mm_add_ps(fiy1SSE,ty1SSE);
            fiy2SSE          = _mm_add_ps(fiy2SSE,ty2SSE);
            fiy3SSE          = _mm_add_ps(fiy3SSE,ty3SSE);
            fiz0SSE          = _mm_add_ps(fiz0SSE,tz0SSE);
            fiz1SSE          = _mm_add_ps(fiz1SSE,tz1SSE);
            fiz2SSE          = _mm_add_ps(fiz2SSE,tz2SSE);
            fiz3SSE          = _mm_add_ps(fiz3SSE,tz3SSE);
            
            /* Decrement j atom force */
            fjxSSE           = _mm_load_ps(fx_align+j);
            fjySSE           = _mm_load_ps(fy_align+j);
            fjzSSE           = _mm_load_ps(fz_align+j);

            tx0SSE           = _mm_add_ps(tx0SSE,tx1SSE);
            tx2SSE           = _mm_add_ps(tx2SSE,tx3SSE);
            tx0SSE           = _mm_add_ps(tx0SSE,tx2SSE);
            ty0SSE           = _mm_add_ps(ty0SSE,ty1SSE);
            ty2SSE           = _mm_add_ps(ty2SSE,ty3SSE);
            ty0SSE           = _mm_add_ps(ty0SSE,ty2SSE);
            tz0SSE           = _mm_add_ps(tz0SSE,tz1SSE);
            tz2SSE           = _mm_add_ps(tz2SSE,tz3SSE);
            tz0SSE           = _mm_add_ps(tz0SSE,tz2SSE);
            
            fjxSSE           = _mm_sub_ps(fjxSSE,tx0SSE);
            fjySSE           = _mm_sub_ps(fjySSE,ty0SSE);
            fjzSSE           = _mm_sub_ps(fjzSSE,tz0SSE);
            
            _mm_store_ps(fx_align+j,fjxSSE);
            _mm_store_ps(fy_align+j,fjySSE);
            _mm_store_ps(fz_align+j,fjzSSE);
            
            /* Inner loop uses 38 flops/iteration */
        }

        for(j=nj1; j<nj2; j+=UNROLLJ)
        {                      
            /* load j atom coordinates */
            jxSSE            = _mm_load_ps(x_align+j);
            jySSE            = _mm_load_ps(y_align+j);
            jzSSE            = _mm_load_ps(z_align+j);
            
            /* Calculate distance */
            dx0SSE            = _mm_sub_ps(ix0SSE,jxSSE);
            dy0SSE            = _mm_sub_ps(iy0SSE,jySSE);
            dz0SSE            = _mm_sub_ps(iz0SSE,jzSSE);
            dx1SSE            = _mm_sub_ps(ix1SSE,jxSSE);
            dy1SSE            = _mm_sub_ps(iy1SSE,jySSE);
            dz1SSE            = _mm_sub_ps(iz1SSE,jzSSE);
            dx2SSE            = _mm_sub_ps(ix2SSE,jxSSE);
            dy2SSE            = _mm_sub_ps(iy2SSE,jySSE);
            dz2SSE            = _mm_sub_ps(iz2SSE,jzSSE);
            dx3SSE            = _mm_sub_ps(ix3SSE,jxSSE);
            dy3SSE            = _mm_sub_ps(iy3SSE,jySSE);
            dz3SSE            = _mm_sub_ps(iz3SSE,jzSSE);
            
            rsq0SSE           = _mm_mul_ps(dx0SSE,dx0SSE);
            rsq1SSE           = _mm_mul_ps(dx1SSE,dx1SSE);
            rsq2SSE           = _mm_mul_ps(dx2SSE,dx2SSE);
            rsq3SSE           = _mm_mul_ps(dx3SSE,dx3SSE);
            tmpA0SSE          = _mm_mul_ps(dy0SSE,dy0SSE);
            tmpA1SSE          = _mm_mul_ps(dy1SSE,dy1SSE);
            tmpA2SSE          = _mm_mul_ps(dy2SSE,dy2SSE);
            tmpA3SSE          = _mm_mul_ps(dy3SSE,dy3SSE);
            tmpB0SSE          = _mm_mul_ps(dz0SSE,dz0SSE);
            tmpB1SSE          = _mm_mul_ps(dz1SSE,dz1SSE);
            tmpB2SSE          = _mm_mul_ps(dz2SSE,dz2SSE);
            tmpB3SSE          = _mm_mul_ps(dz3SSE,dz3SSE);
            rsq0SSE           = _mm_add_ps(rsq0SSE,tmpA0SSE);
            rsq1SSE           = _mm_add_ps(rsq1SSE,tmpA1SSE);
            rsq2SSE           = _mm_add_ps(rsq2SSE,tmpA2SSE);
            rsq3SSE           = _mm_add_ps(rsq3SSE,tmpA3SSE);
            rsq0SSE           = _mm_add_ps(rsq0SSE,tmpB0SSE);
            rsq1SSE           = _mm_add_ps(rsq1SSE,tmpB1SSE);
            rsq2SSE           = _mm_add_ps(rsq2SSE,tmpB2SSE);
            rsq3SSE           = _mm_add_ps(rsq3SSE,tmpB3SSE);
            
            /* Calculate 1/r and 1/r2 */
            rinv0SSE          = gmx_mm_invsqrt_ps(rsq0SSE);
            rinv1SSE          = gmx_mm_invsqrt_ps(rsq1SSE);
            rinv2SSE          = gmx_mm_invsqrt_ps(rsq2SSE);
            rinv3SSE          = gmx_mm_invsqrt_ps(rsq3SSE);
            
            rinv0SSE          = _mm_and_ps(rinv0SSE,imask0SSE);
            rinv1SSE          = _mm_and_ps(rinv1SSE,imask1SSE);
            rinv2SSE          = _mm_and_ps(rinv2SSE,imask2SSE);
            rinv3SSE          = _mm_and_ps(rinv3SSE,imask3SSE);

            /* Load parameters for j atom */
            jqSSE             = _mm_load_ps(q_align+j);
            qq0SSE            = _mm_mul_ps(iq0SSE,jqSSE);
            qq1SSE            = _mm_mul_ps(iq1SSE,jqSSE);
            qq2SSE            = _mm_mul_ps(iq2SSE,jqSSE);
            qq3SSE            = _mm_mul_ps(iq3SSE,jqSSE);
            
            c60SSE            = _mm_load_ps(pvdw0+2*j);
            c61SSE            = _mm_load_ps(pvdw1+2*j);
            c62SSE            = _mm_load_ps(pvdw2+2*j);
            c63SSE            = _mm_load_ps(pvdw3+2*j);
            c120SSE           = _mm_load_ps(pvdw0+2*j+UNROLLJ);
            c121SSE           = _mm_load_ps(pvdw1+2*j+UNROLLJ);
            c122SSE           = _mm_load_ps(pvdw2+2*j+UNROLLJ);
            c123SSE           = _mm_load_ps(pvdw3+2*j+UNROLLJ);
            
            rinvsq0SSE        = _mm_mul_ps(rinv0SSE,rinv0SSE);
            rinvsq1SSE        = _mm_mul_ps(rinv1SSE,rinv1SSE);
            rinvsq2SSE        = _mm_mul_ps(rinv2SSE,rinv2SSE);
            rinvsq3SSE        = _mm_mul_ps(rinv3SSE,rinv3SSE);
            
            /* Coulomb interaction */
            vcoul0SSE         = _mm_mul_ps(qq0SSE,rinv0SSE);
            vcoul1SSE         = _mm_mul_ps(qq1SSE,rinv1SSE);
            vcoul2SSE         = _mm_mul_ps(qq2SSE,rinv2SSE);
            vcoul3SSE         = _mm_mul_ps(qq3SSE,rinv3SSE);

            tmpA0SSE          = _mm_add_ps(vcoul0SSE,vcoul1SSE);
            tmpA2SSE          = _mm_add_ps(vcoul2SSE,vcoul3SSE);
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA2SSE);
            vctotSSE          = _mm_add_ps(vctotSSE,tmpA0SSE);
            
            /* Lennard-Jones interaction */
            rinvsix0SSE       = _mm_mul_ps(rinvsq0SSE,_mm_mul_ps(rinvsq0SSE,rinvsq0SSE));
            rinvsix1SSE       = _mm_mul_ps(rinvsq1SSE,_mm_mul_ps(rinvsq1SSE,rinvsq1SSE));
            rinvsix2SSE       = _mm_mul_ps(rinvsq2SSE,_mm_mul_ps(rinvsq2SSE,rinvsq2SSE));
            rinvsix3SSE       = _mm_mul_ps(rinvsq3SSE,_mm_mul_ps(rinvsq3SSE,rinvsq3SSE));
            Vvdw60SSE         = _mm_mul_ps(c60SSE,rinvsix0SSE);
            Vvdw61SSE         = _mm_mul_ps(c61SSE,rinvsix1SSE);
            Vvdw62SSE         = _mm_mul_ps(c62SSE,rinvsix2SSE);
            Vvdw63SSE         = _mm_mul_ps(c63SSE,rinvsix3SSE);
            Vvdw120SSE        = _mm_mul_ps(c120SSE,_mm_mul_ps(rinvsix0SSE,rinvsix0SSE));
            Vvdw121SSE        = _mm_mul_ps(c121SSE,_mm_mul_ps(rinvsix1SSE,rinvsix1SSE));
            Vvdw122SSE        = _mm_mul_ps(c122SSE,_mm_mul_ps(rinvsix2SSE,rinvsix2SSE));
            Vvdw123SSE        = _mm_mul_ps(c123SSE,_mm_mul_ps(rinvsix3SSE,rinvsix3SSE));
            
            tmpA0SSE          = _mm_sub_ps(Vvdw120SSE,Vvdw60SSE);
            tmpA1SSE          = _mm_sub_ps(Vvdw121SSE,Vvdw61SSE);
            tmpA2SSE          = _mm_sub_ps(Vvdw122SSE,Vvdw62SSE);
            tmpA3SSE          = _mm_sub_ps(Vvdw123SSE,Vvdw63SSE);
            
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA1SSE);
            tmpA2SSE          = _mm_add_ps(tmpA2SSE,tmpA3SSE);
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA2SSE);
            VvdwtotSSE       = _mm_add_ps(VvdwtotSSE, tmpA0SSE);
            
            fscal0SSE         = _mm_mul_ps(rinvsq0SSE, 
                                           _mm_add_ps(vcoul0SSE,
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw120SSE),
                                                                 _mm_mul_ps(sixSSE,Vvdw60SSE))));
            fscal1SSE         = _mm_mul_ps(rinvsq1SSE, 
                                           _mm_add_ps(vcoul1SSE,
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw121SSE),
                                                                 _mm_mul_ps(sixSSE,Vvdw61SSE))));
            fscal2SSE         = _mm_mul_ps(rinvsq2SSE, 
                                           _mm_add_ps(vcoul2SSE,
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw122SSE),
                                                                 _mm_mul_ps(sixSSE,Vvdw62SSE))));
            fscal3SSE         = _mm_mul_ps(rinvsq3SSE, 
                                           _mm_add_ps(vcoul3SSE,
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw123SSE),
                                                                 _mm_mul_ps(sixSSE,Vvdw63SSE))));
            
            /* Calculate temporary vectorial force */
            tx0SSE            = _mm_mul_ps(fscal0SSE,dx0SSE);
            tx1SSE            = _mm_mul_ps(fscal1SSE,dx1SSE);
            tx2SSE            = _mm_mul_ps(fscal2SSE,dx2SSE);
            tx3SSE            = _mm_mul_ps(fscal3SSE,dx3SSE);
            ty0SSE            = _mm_mul_ps(fscal0SSE,dy0SSE);
            ty1SSE            = _mm_mul_ps(fscal1SSE,dy1SSE);
            ty2SSE            = _mm_mul_ps(fscal2SSE,dy2SSE);
            ty3SSE            = _mm_mul_ps(fscal3SSE,dy3SSE);
            tz0SSE            = _mm_mul_ps(fscal0SSE,dz0SSE);
            tz1SSE            = _mm_mul_ps(fscal1SSE,dz1SSE);
            tz2SSE            = _mm_mul_ps(fscal2SSE,dz2SSE);
            tz3SSE            = _mm_mul_ps(fscal3SSE,dz3SSE);
            
            /* Increment i atom force */
            fix0SSE          = _mm_add_ps(fix0SSE,tx0SSE);
            fix1SSE          = _mm_add_ps(fix1SSE,tx1SSE);
            fix2SSE          = _mm_add_ps(fix2SSE,tx2SSE);
            fix3SSE          = _mm_add_ps(fix3SSE,tx3SSE);
            fiy0SSE          = _mm_add_ps(fiy0SSE,ty0SSE);
            fiy1SSE          = _mm_add_ps(fiy1SSE,ty1SSE);
            fiy2SSE          = _mm_add_ps(fiy2SSE,ty2SSE);
            fiy3SSE          = _mm_add_ps(fiy3SSE,ty3SSE);
            fiz0SSE          = _mm_add_ps(fiz0SSE,tz0SSE);
            fiz1SSE          = _mm_add_ps(fiz1SSE,tz1SSE);
            fiz2SSE          = _mm_add_ps(fiz2SSE,tz2SSE);
            fiz3SSE          = _mm_add_ps(fiz3SSE,tz3SSE);
            
            /* Decrement j atom force */
            fjxSSE           = _mm_load_ps(fx_align+j);
            fjySSE           = _mm_load_ps(fy_align+j);
            fjzSSE           = _mm_load_ps(fz_align+j);
            
            tx0SSE           = _mm_add_ps(tx0SSE,tx1SSE);
            tx2SSE           = _mm_add_ps(tx2SSE,tx3SSE);
            tx0SSE           = _mm_add_ps(tx0SSE,tx2SSE);
            ty0SSE           = _mm_add_ps(ty0SSE,ty1SSE);
            ty2SSE           = _mm_add_ps(ty2SSE,ty3SSE);
            ty0SSE           = _mm_add_ps(ty0SSE,ty2SSE);
            tz0SSE           = _mm_add_ps(tz0SSE,tz1SSE);
            tz2SSE           = _mm_add_ps(tz2SSE,tz3SSE);
            tz0SSE           = _mm_add_ps(tz0SSE,tz2SSE);
            
            fjxSSE           = _mm_sub_ps(fjxSSE,tx0SSE);
            fjySSE           = _mm_sub_ps(fjySSE,ty0SSE);
            fjzSSE           = _mm_sub_ps(fjzSSE,tz0SSE);
            
            _mm_store_ps(fx_align+j,fjxSSE);
            _mm_store_ps(fy_align+j,fjySSE);
            _mm_store_ps(fz_align+j,fjzSSE);
            
            /* Inner loop uses 38 flops/iteration */
        }

        for(j=nj2; j<nj3; j+=UNROLLJ)
        {
            jmask0SSE = _mm_load_ps((real *)emask0);
            jmask1SSE = _mm_load_ps((real *)emask1);
            jmask2SSE = _mm_load_ps((real *)emask2);
            jmask3SSE = _mm_load_ps((real *)emask3);
            emask0 += UNROLLJ;
            emask1 += UNROLLJ;
            emask2 += UNROLLJ;
            emask3 += UNROLLJ;
            
            /* load j atom coordinates */
            jxSSE            = _mm_load_ps(x_align+j);
            jySSE            = _mm_load_ps(y_align+j);
            jzSSE            = _mm_load_ps(z_align+j);
            
            /* Calculate distance */
            dx0SSE            = _mm_sub_ps(ix0SSE,jxSSE);
            dy0SSE            = _mm_sub_ps(iy0SSE,jySSE);
            dz0SSE            = _mm_sub_ps(iz0SSE,jzSSE);
            dx1SSE            = _mm_sub_ps(ix1SSE,jxSSE);
            dy1SSE            = _mm_sub_ps(iy1SSE,jySSE);
            dz1SSE            = _mm_sub_ps(iz1SSE,jzSSE);
            dx2SSE            = _mm_sub_ps(ix2SSE,jxSSE);
            dy2SSE            = _mm_sub_ps(iy2SSE,jySSE);
            dz2SSE            = _mm_sub_ps(iz2SSE,jzSSE);
            dx3SSE            = _mm_sub_ps(ix3SSE,jxSSE);
            dy3SSE            = _mm_sub_ps(iy3SSE,jySSE);
            dz3SSE            = _mm_sub_ps(iz3SSE,jzSSE);
            
            rsq0SSE           = _mm_mul_ps(dx0SSE,dx0SSE);
            rsq1SSE           = _mm_mul_ps(dx1SSE,dx1SSE);
            rsq2SSE           = _mm_mul_ps(dx2SSE,dx2SSE);
            rsq3SSE           = _mm_mul_ps(dx3SSE,dx3SSE);
            tmpA0SSE          = _mm_mul_ps(dy0SSE,dy0SSE);
            tmpA1SSE          = _mm_mul_ps(dy1SSE,dy1SSE);
            tmpA2SSE          = _mm_mul_ps(dy2SSE,dy2SSE);
            tmpA3SSE          = _mm_mul_ps(dy3SSE,dy3SSE);
            tmpB0SSE          = _mm_mul_ps(dz0SSE,dz0SSE);
            tmpB1SSE          = _mm_mul_ps(dz1SSE,dz1SSE);
            tmpB2SSE          = _mm_mul_ps(dz2SSE,dz2SSE);
            tmpB3SSE          = _mm_mul_ps(dz3SSE,dz3SSE);
            rsq0SSE           = _mm_add_ps(rsq0SSE,tmpA0SSE);
            rsq1SSE           = _mm_add_ps(rsq1SSE,tmpA1SSE);
            rsq2SSE           = _mm_add_ps(rsq2SSE,tmpA2SSE);
            rsq3SSE           = _mm_add_ps(rsq3SSE,tmpA3SSE);
            rsq0SSE           = _mm_add_ps(rsq0SSE,tmpB0SSE);
            rsq1SSE           = _mm_add_ps(rsq1SSE,tmpB1SSE);
            rsq2SSE           = _mm_add_ps(rsq2SSE,tmpB2SSE);
            rsq3SSE           = _mm_add_ps(rsq3SSE,tmpB3SSE);
            
            jmask0SSE         = _mm_and_ps(jmask0SSE,imask0SSE);
            jmask1SSE         = _mm_and_ps(jmask1SSE,imask1SSE);
            jmask2SSE         = _mm_and_ps(jmask2SSE,imask2SSE);
            jmask3SSE         = _mm_and_ps(jmask3SSE,imask3SSE);
            
            /* Calculate 1/r and 1/r2 */
            rinv0SSE          = gmx_mm_invsqrt_ps(rsq0SSE);
            rinv1SSE          = gmx_mm_invsqrt_ps(rsq1SSE);
            rinv2SSE          = gmx_mm_invsqrt_ps(rsq2SSE);
            rinv3SSE          = gmx_mm_invsqrt_ps(rsq3SSE);
            rinv0SSE          = _mm_and_ps(rinv0SSE,jmask0SSE);
            rinv1SSE          = _mm_and_ps(rinv1SSE,jmask1SSE);
            rinv2SSE          = _mm_and_ps(rinv2SSE,jmask2SSE);
            rinv3SSE          = _mm_and_ps(rinv3SSE,jmask3SSE);
            
            /* Load parameters for j atom */
            jqSSE             = _mm_load_ps(q_align+j);
            qq0SSE            = _mm_mul_ps(iq0SSE,jqSSE);
            qq1SSE            = _mm_mul_ps(iq1SSE,jqSSE);
            qq2SSE            = _mm_mul_ps(iq2SSE,jqSSE);
            qq3SSE            = _mm_mul_ps(iq3SSE,jqSSE);
            
            c60SSE            = _mm_load_ps(pvdw0+2*j);
            c61SSE            = _mm_load_ps(pvdw1+2*j);
            c62SSE            = _mm_load_ps(pvdw2+2*j);
            c63SSE            = _mm_load_ps(pvdw3+2*j);
            c120SSE           = _mm_load_ps(pvdw0+2*j+UNROLLJ);
            c121SSE           = _mm_load_ps(pvdw1+2*j+UNROLLJ);
            c122SSE           = _mm_load_ps(pvdw2+2*j+UNROLLJ);
            c123SSE           = _mm_load_ps(pvdw3+2*j+UNROLLJ);
            
            rinvsq0SSE        = _mm_mul_ps(rinv0SSE,rinv0SSE);
            rinvsq1SSE        = _mm_mul_ps(rinv1SSE,rinv1SSE);
            rinvsq2SSE        = _mm_mul_ps(rinv2SSE,rinv2SSE);
            rinvsq3SSE        = _mm_mul_ps(rinv3SSE,rinv3SSE);
            
            /* Coulomb interaction */
            vcoul0SSE         = _mm_mul_ps(qq0SSE,rinv0SSE);
            vcoul1SSE         = _mm_mul_ps(qq1SSE,rinv1SSE);
            vcoul2SSE         = _mm_mul_ps(qq2SSE,rinv2SSE);
            vcoul3SSE         = _mm_mul_ps(qq3SSE,rinv3SSE);

            tmpA0SSE          = _mm_add_ps(vcoul0SSE,vcoul1SSE);
            tmpA2SSE          = _mm_add_ps(vcoul2SSE,vcoul3SSE);
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA2SSE);
            vctotSSE          = _mm_add_ps(vctotSSE,tmpA0SSE);
            
            /* Lennard-Jones interaction */
            rinvsix0SSE       = _mm_mul_ps(rinvsq0SSE,_mm_mul_ps(rinvsq0SSE,rinvsq0SSE));
            rinvsix1SSE       = _mm_mul_ps(rinvsq1SSE,_mm_mul_ps(rinvsq1SSE,rinvsq1SSE));
            rinvsix2SSE       = _mm_mul_ps(rinvsq2SSE,_mm_mul_ps(rinvsq2SSE,rinvsq2SSE));
            rinvsix3SSE       = _mm_mul_ps(rinvsq3SSE,_mm_mul_ps(rinvsq3SSE,rinvsq3SSE));
            Vvdw60SSE         = _mm_mul_ps(c60SSE,rinvsix0SSE);
            Vvdw61SSE         = _mm_mul_ps(c61SSE,rinvsix1SSE);
            Vvdw62SSE         = _mm_mul_ps(c62SSE,rinvsix2SSE);
            Vvdw63SSE         = _mm_mul_ps(c63SSE,rinvsix3SSE);
            Vvdw120SSE        = _mm_mul_ps(c120SSE,_mm_mul_ps(rinvsix0SSE,rinvsix0SSE));
            Vvdw121SSE        = _mm_mul_ps(c121SSE,_mm_mul_ps(rinvsix1SSE,rinvsix1SSE));
            Vvdw122SSE        = _mm_mul_ps(c122SSE,_mm_mul_ps(rinvsix2SSE,rinvsix2SSE));
            Vvdw123SSE        = _mm_mul_ps(c123SSE,_mm_mul_ps(rinvsix3SSE,rinvsix3SSE));
            
            tmpA0SSE          = _mm_sub_ps(Vvdw120SSE,Vvdw60SSE);
            tmpA1SSE          = _mm_sub_ps(Vvdw121SSE,Vvdw61SSE);
            tmpA2SSE          = _mm_sub_ps(Vvdw122SSE,Vvdw62SSE);
            tmpA3SSE          = _mm_sub_ps(Vvdw123SSE,Vvdw63SSE);
            
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA1SSE);
            tmpA2SSE          = _mm_add_ps(tmpA2SSE,tmpA3SSE);
            tmpA0SSE          = _mm_add_ps(tmpA0SSE,tmpA2SSE);
            VvdwtotSSE       = _mm_add_ps(VvdwtotSSE, tmpA0SSE);
            
            fscal0SSE         = _mm_mul_ps(rinvsq0SSE, 
                                           _mm_add_ps(vcoul0SSE,
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw120SSE),
                                                                 _mm_mul_ps(sixSSE,Vvdw60SSE))));
            fscal1SSE         = _mm_mul_ps(rinvsq1SSE, 
                                           _mm_add_ps(vcoul1SSE,
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw121SSE),
                                                                 _mm_mul_ps(sixSSE,Vvdw61SSE))));
            fscal2SSE         = _mm_mul_ps(rinvsq2SSE, 
                                           _mm_add_ps(vcoul2SSE,
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw122SSE),
                                                                 _mm_mul_ps(sixSSE,Vvdw62SSE))));
            fscal3SSE         = _mm_mul_ps(rinvsq3SSE, 
                                           _mm_add_ps(vcoul3SSE,
                                                      _mm_sub_ps(_mm_mul_ps(twelveSSE,Vvdw123SSE),
                                                                 _mm_mul_ps(sixSSE,Vvdw63SSE))));
            
            /* Calculate temporary vectorial force */
            tx0SSE            = _mm_mul_ps(fscal0SSE,dx0SSE);
            tx1SSE            = _mm_mul_ps(fscal1SSE,dx1SSE);
            tx2SSE            = _mm_mul_ps(fscal2SSE,dx2SSE);
            tx3SSE            = _mm_mul_ps(fscal3SSE,dx3SSE);
            ty0SSE            = _mm_mul_ps(fscal0SSE,dy0SSE);
            ty1SSE            = _mm_mul_ps(fscal1SSE,dy1SSE);
            ty2SSE            = _mm_mul_ps(fscal2SSE,dy2SSE);
            ty3SSE            = _mm_mul_ps(fscal3SSE,dy3SSE);
            tz0SSE            = _mm_mul_ps(fscal0SSE,dz0SSE);
            tz1SSE            = _mm_mul_ps(fscal1SSE,dz1SSE);
            tz2SSE            = _mm_mul_ps(fscal2SSE,dz2SSE);
            tz3SSE            = _mm_mul_ps(fscal3SSE,dz3SSE);
            
            /* Increment i atom force */
            fix0SSE          = _mm_add_ps(fix0SSE,tx0SSE);
            fix1SSE          = _mm_add_ps(fix1SSE,tx1SSE);
            fix2SSE          = _mm_add_ps(fix2SSE,tx2SSE);
            fix3SSE          = _mm_add_ps(fix3SSE,tx3SSE);
            fiy0SSE          = _mm_add_ps(fiy0SSE,ty0SSE);
            fiy1SSE          = _mm_add_ps(fiy1SSE,ty1SSE);
            fiy2SSE          = _mm_add_ps(fiy2SSE,ty2SSE);
            fiy3SSE          = _mm_add_ps(fiy3SSE,ty3SSE);
            fiz0SSE          = _mm_add_ps(fiz0SSE,tz0SSE);
            fiz1SSE          = _mm_add_ps(fiz1SSE,tz1SSE);
            fiz2SSE          = _mm_add_ps(fiz2SSE,tz2SSE);
            fiz3SSE          = _mm_add_ps(fiz3SSE,tz3SSE);
            
            /* Decrement j atom force */
            fjxSSE           = _mm_load_ps(fx_align+j);
            fjySSE           = _mm_load_ps(fy_align+j);
            fjzSSE           = _mm_load_ps(fz_align+j);
            
            tx0SSE           = _mm_add_ps(tx0SSE,tx1SSE);
            tx2SSE           = _mm_add_ps(tx2SSE,tx3SSE);
            tx0SSE           = _mm_add_ps(tx0SSE,tx2SSE);
            ty0SSE           = _mm_add_ps(ty0SSE,ty1SSE);
            ty2SSE           = _mm_add_ps(ty2SSE,ty3SSE);
            ty0SSE           = _mm_add_ps(ty0SSE,ty2SSE);
            tz0SSE           = _mm_add_ps(tz0SSE,tz1SSE);
            tz2SSE           = _mm_add_ps(tz2SSE,tz3SSE);
            tz0SSE           = _mm_add_ps(tz0SSE,tz2SSE);
            
            fjxSSE           = _mm_sub_ps(fjxSSE,tx0SSE);
            fjySSE           = _mm_sub_ps(fjySSE,ty0SSE);
            fjzSSE           = _mm_sub_ps(fjzSSE,tz0SSE);
            
            _mm_store_ps(fx_align+j,fjxSSE);
            _mm_store_ps(fy_align+j,fjySSE);
            _mm_store_ps(fz_align+j,fjzSSE);
            
            /* Inner loop uses 38 flops/iteration */
        }
        
		/* Add i forces to mem and shifted force list */
        accumulate_SSE_force(fix0SSE,fiy0SSE,fiz0SSE,f+3*i);
        accumulate_SSE_force(fix1SSE,fiy1SSE,fiz1SSE,f+3*(i+1));
        accumulate_SSE_force(fix2SSE,fiy2SSE,fiz2SSE,f+3*(i+2));
        accumulate_SSE_force(fix3SSE,fiy3SSE,fiz3SSE,f+3*(i+3));
		
		/* Add potential energies to the group for this list */
		ggid             = 0;         
        
		_mm_storeu_ps(tmpsum,vctotSSE);
		Vc[ggid]         = Vc[ggid] + tmpsum[0]+tmpsum[1]+tmpsum[2]+tmpsum[3];
		
		_mm_storeu_ps(tmpsum,VvdwtotSSE);
		Vvdw[ggid]       = Vvdw[ggid] + tmpsum[0]+tmpsum[1]+tmpsum[2]+tmpsum[3];
		
		/* Outer loop uses 6 flops/iteration */
	}    

	for(i=0;i<natoms;i++)
	{
		f[3*i]   += fx_align[i] + fx_align[natoms+i];
		f[3*i+1] += fy_align[i] + fy_align[natoms+i];
		f[3*i+2] += fz_align[i] + fz_align[natoms+i];
	}
    
    /* Write outer/inner iteration count to pointers */
    *outeriter       = ni1-ni0;         
    *inneriter       = (ni1-ni0)*natoms/2;         
}

#undef SIMD_WIDTH
#undef UNROLLI   
#undef UNROLLJ   

