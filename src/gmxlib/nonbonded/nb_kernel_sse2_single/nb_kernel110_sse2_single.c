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
#if 0

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "sse_common_single.h"


/*
 * Gromacs nonbonded kernel nb_kernel112
 * Coulomb interaction:     Normal Coulomb
 * VdW interaction:         Lennard-Jones
 * water optimization:      pairs of SPC/TIP3P interactions
 * Calculate forces:        yes
 */
void nb_kernel110_sse2_single(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    float *         shiftvec,
                    float *         fshift,
                    int *           gid,
                    float *         pos,
                    float *         faction,
                    float *         charge,
                    float *         p_facel,
                    float *         p_krf,
                    float *         p_crf,
                    float *         Vc,
                    int *           type,
                    int *           p_ntype,
                    float *         vdwparam,
                    float *         Vvdw,
                    float *         p_tabscale,
                    float *         VFtab,
                    float *         invsqrta,
                    float *         dvda,
                    float *         p_gbtabscale,
                    float *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    float *         work)
{
    int           nri,ntype,nthreads,offset;
    float         facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    float         shX,shY,shZ;
    int           tj;
    float         qO,qH;
	int           j3A,j3B,j3C,j3D,jnrA,jnrB,jnrC,jnrD;
	float         tfloat1[4],tfloat2[4],tfloat3[4];

	__m128     mask1 = _mm_castsi128_ps( _mm_set_epi32(0, 0, 0, 0xffffffff) );
	__m128     mask2 = _mm_castsi128_ps( _mm_set_epi32(0, 0, 0xffffffff, 0xffffffff) );
	__m128     mask3 = _mm_castsi128_ps( _mm_set_epi32(0, 0xffffffff, 0xffffffff, 0xffffffff) );
	__m128     mask;
	
	__m128     qqOO,qqOH,qqHH,c6,c12,Vvdwtot,vctot;
	__m128     dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
	__m128     dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
	__m128     dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;
	__m128     fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
	__m128     fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
	__m128     rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
	__m128     tx,ty,tz;
	__m128     ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
	__m128     jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;  
	__m128     t1,t2,t3;
	__m128     rinv11,rinv12,rinv13,rinv21,rinv22,rinv23,rinv31,rinv32,rinv33;
	__m128     fscal11,fscal12,fscal13,fscal21,fscal22,fscal23,fscal31,fscal32,fscal33;
	
	__m128     half   = _mm_set1_ps(0.5);
	__m128     three  = _mm_set1_ps(3.0);
	__m128     six    = _mm_set1_ps(6.0);
	__m128     twelve = _mm_set1_ps(12.0);
	
    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    ii               = iinr[0];        
    qO               = charge[ii];     
    qH               = charge[ii+1];   
    qqOO             = _mm_set1_ps(facel*qO*qO);    
    qqOH             = _mm_set1_ps(facel*qO*qH);    
    qqHH             = _mm_set1_ps(facel*qH*qH);    
    tj               = 2*(ntype+1)*type[ii];
    c6               = _mm_set1_ps(vdwparam[tj]);   
    c12              = _mm_set1_ps(vdwparam[tj+1]); 

    nj1              = 0;              
    
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
		
		ix1      = _mm_set1_ps(shX+pos[ii3+0]);
		iy1      = _mm_set1_ps(shY+pos[ii3+1]);
		iz1      = _mm_set1_ps(shZ+pos[ii3+2]); 

        vctot            = _mm_setzero_ps();              
        Vvdwtot          = _mm_setzero_ps();               
        fix1             = _mm_setzero_ps(); 
        fiy1             = _mm_setzero_ps();               
        fiz1             = _mm_setzero_ps();               
        
        for(k=nj0; k<nj1-3; k+=4)
        {
			jnrA             = jjnr[k];
			jnrB             = jjnr[k+1];
			jnrC             = jjnr[k+2];
			jnrD             = jjnr[k+3];
			
			j3A              = jnrA * 3;
			j3B              = jnrB * 3;
			j3C              = jnrC * 3;
			j3D              = jnrD * 3;
			
			GMX_MM_LOAD_4JCOORD_1ATOM_PS(pos+j3A,pos+j3B,pos+j3C,pos+j3D,jx1,jy1,jz1);
			
			dx11   = _mm_sub_ps(ix1,jx1);
			dy11   = _mm_sub_ps(iy1,jy1);
			dz11   = _mm_sub_ps(iz1,jz1);
			
			rsq11  = gmx_mm_scalarprod_ps(dx11,dy11,dz11);
			
			/* invsqrt */								
			rinv11 = gmx_mm_invsqrt_ps(rsq11);

			/***********************************/
			/* INTERACTION SECTION STARTS HERE */
			/***********************************/
			
			fscal11 = _mm_add_ps(gmx_mm_int_coulomb_ps(rinv11,qqOO,&vctot),
								 gmx_mm_int_lj_ps(rinv11,c6,c12,&Vvdwtot));
			fscal11 = _mm_mul_ps(fscal11,_mm_mul_ps(rinv11,rinv11)); /*  *=rinvsq */
			
			/***********************************/
			/*  INTERACTION SECTION ENDS HERE  */
			/***********************************/
			
			/* Calculate vectorial forces from fscal */
			
			tx      = _mm_mul_ps(fscal11,dx11); 
			ty      = _mm_mul_ps(fscal11,dy11); 
			tz      = _mm_mul_ps(fscal11,dz11); 
			fix1    = _mm_add_ps(fix1,tx);      
			fiy1    = _mm_add_ps(fiy1,ty);      
			fiz1    = _mm_add_ps(fiz1,tz);     
			fjx1    = tx;
			fjy1    = ty;
			fjz1    = tz;
			
			/* Store j forces back */
			GMX_MM_UPDATE_4JCOORD_1ATOM_PS(faction+j3A,faction+j3B,faction+j3C,faction+j3D,fjx1,fjy1,fjz1);
        }
		
        /* Treat offset particles here */
		if(k<nj1)
        {
	
			if(k<nj1-2)
			{
				jnrA             = jjnr[k];
				jnrB             = jjnr[k+1];
				jnrC             = jjnr[k+2];
				j3A              = jnrA * 3;
				j3B              = jnrB * 3;
				j3C              = jnrC * 3;
				GMX_MM_LOAD_3JCOORD_1ATOM_PS(pos+j3A,pos+j3B,pos+j3C,jx1,jy1,jz1);                         
			}
			else if(k<nj1-1)
			{
				jnrA             = jjnr[k];
				jnrB             = jjnr[k+1];
				j3A              = jnrA * 3;
				j3B              = jnrB * 3;
				GMX_MM_LOAD_2JCOORD_1ATOM_PS(pos+j3A,pos+j3B,jx1,jy1,jz1);
			}
			else
			{
				jnrA             = jjnr[k];
				j3A              = jnrA * 3;
				GMX_MM_LOAD_1JCOORD_1ATOM_PS(pos+j3A,jx1,jy1,jz1);
			}
			
			dx11   = _mm_sub_ps(ix1,jx1);
			dy11   = _mm_sub_ps(iy1,jy1);
			dz11   = _mm_sub_ps(iz1,jz1);
			
			rsq11  = gmx_mm_scalarprod_ps(dx11,dy11,dz11);
			
			/* invsqrt */					
			rinv11 = _mm_and_ps(mask,gmx_mm_invsqrt_ps(rsq11));
			
			/***********************************/
			/* INTERACTION SECTION STARTS HERE */
			/***********************************/
			
			fscal11 = _mm_add_ps(gmx_mm_int_coulomb_ps(rinv11,qqOO,&vctot),
								 gmx_mm_int_lj_ps(rinv11,c6,c12,&Vvdwtot));
			fscal11 = _mm_mul_ps(fscal11,_mm_mul_ps(rinv11,rinv11)); /*  *=rinvsq */
			
			/***********************************/
			/*  INTERACTION SECTION ENDS HERE  */
			/***********************************/
			
			/* Calculate vectorial forces */
			tx      = _mm_mul_ps(fscal11,dx11);
			ty      = _mm_mul_ps(fscal11,dy11);
			tz      = _mm_mul_ps(fscal11,dz11);
			fix1    = _mm_add_ps(fix1,tx);
			fiy1    = _mm_add_ps(fiy1,ty);
			fiz1    = _mm_add_ps(fiz1,tz);
			fjx1    = tx;
			fjy1    = ty;
			fjz1    = tz;

			if(k<nj1-2)
			{
				GMX_MM_UPDATE_3JCOORD_1ATOM_PS(faction+j3A,faction+j3B,faction+j3C,fjx1,fjy1,fjz1);
			}
			else if(k<nj1-1)
			{
				GMX_MM_UPDATE_2JCOORD_1ATOM_PS(faction+j3A,faction+j3B,fjx1,fjy1,fjz1);
			}
			else
			{
				GMX_MM_UPDATE_1JCOORD_1ATOM_PS(faction+j3A,fjx1,fjy1,fjz1);
			}
		}
		
		gmx_mm_update_iforce_1atom_ps(&fix1,&fiy1,&fiz1,faction+ii3,fshift+is3);
		
        ggid             = gid[n];         
	GMX_MM_UPDATE_2POT_PS(vctot,Vc+ggid,Vvdwtot,Vvdw+ggid);
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}

#else
int nb_kernel110_sse2_single_dummy;
#endif
