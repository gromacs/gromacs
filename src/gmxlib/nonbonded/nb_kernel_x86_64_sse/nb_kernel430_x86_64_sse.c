/*
 * Copyright (c) Erik Lindahl, David van der Spoel 2003
 * 
 * This file is generated automatically at compile time
 * by the program mknb in the Gromacs distribution.
 *
 * Options used when generation this file:
 * Language:         c
 * Precision:        single
 * Threads:          no
 * Software invsqrt: yes
 * Prefetch forces:  no
 * Comments:         no
 */
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#include<math.h>
#include<vec.h>


#include <xmmintrin.h>
#include <emmintrin.h>

#include <gmx_sse2_single.h>

/* get gmx_gbdata_t */
#include "../nb_kerneltype.h"

#include "nb_kernel430_x86_64_sse.h"



void nb_kernel430_x86_64_sse(int *           p_nri,
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
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale;
    int           n,ii,is3,ii3,k,nj0,nj1,ggid;
    float         shX,shY,shZ;
	int			  offset,nti;
	int           jnr,jnr2,jnr3,jnr4,j3,j23,j33,j43;
	int           tj,tj2,tj3,tj4;
	gmx_gbdata_t *gbdata;
	float *        gpol;
	
	__m128   iq,qq,q,isai;
	__m128   ix,iy,iz;
	__m128   jx,jy,jz;
	__m128   xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8;
	__m128   dx1,dy1,dz1;
	__m128   vctot,Vvdwtot,dvdasum;
	__m128   fix,fiy,fiz,rsq;
	__m128   t1,t2,t3;
	__m128   rinv,isaj,isaprod;
	__m128   vcoul,fscal,gbscale,c6,c12;
	__m128   rinvsq,dvdaj,r,rt;
	__m128   eps,eps2,Y,F,G,H,Geps,Heps2;
	__m128   Fp,VV,FF,vgb,fijC,fijD,fijR,dvdatmp;
	__m128   rinvsix,Vvdw6,Vvdw12,Vvdwtmp,vgbtot,n0f;
	__m128   fac_sse,tabscale_sse,gbtabscale_sse;
	
	__m128i  n0, nnn;
	const __m128 neg    = _mm_set1_ps(-1.0f);
	const __m128 zero   = _mm_set1_ps(0.0f);
    const __m128 half   = _mm_set1_ps(0.5f);
	const __m128 two    = _mm_set1_ps(2.0f);
	const __m128 three  = _mm_set1_ps(3.0f);
	const __m128 six    = _mm_set1_ps(6.0f);
    const __m128 twelwe = _mm_set1_ps(12.0f);
	
	__m128i four        = _mm_set1_epi32(4);
	__m128i maski       = _mm_set_epi32(0, 0xffffffff, 0xffffffff, 0xffffffff);     
	__m128i mask        = _mm_set_epi32(0, 0xffffffff, 0xffffffff, 0xffffffff);   
	
	float vct,vdwt,vgbt,dva,isai_f;
	
	gbdata          = (gmx_gbdata_t *)work;
	gpol            = gbdata->gpol;
		
	nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = (*p_facel) * ((1.0/gbdata->epsilon_r) - (1.0/gbdata->gb_epsilon_solvent));       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    gbtabscale       = *p_gbtabscale;  
    nj1              = 0;

	/* Splat variables */
	fac_sse        = _mm_load1_ps(&facel);
	tabscale_sse   = _mm_load1_ps(&tabscale);
	gbtabscale_sse = _mm_load1_ps(&gbtabscale);
	
	
	/* Keep the compiler happy */
	Vvdwtmp = _mm_setzero_ps();
	dvdatmp = _mm_setzero_ps();
	vcoul   = _mm_setzero_ps();
	vgb     = _mm_setzero_ps();
	t1      = _mm_setzero_ps();
	t2      = _mm_setzero_ps();
	t3      = _mm_setzero_ps();
	xmm1    = _mm_setzero_ps();
	xmm2    = _mm_setzero_ps();
	xmm3    = _mm_setzero_ps();
	xmm4    = _mm_setzero_ps();
	j23     = j33  = 0;
	jnr2    = jnr3 = 0;
	
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
        
		ix               = _mm_set1_ps(shX+pos[ii3+0]);
		iy               = _mm_set1_ps(shY+pos[ii3+1]);
		iz               = _mm_set1_ps(shZ+pos[ii3+2]);
		
		iq               = _mm_load1_ps(charge+ii);
		iq               = _mm_mul_ps(iq,fac_sse);
		
        isai_f           = invsqrta[ii];   
		isai             = _mm_load1_ps(&isai_f);
		
		nti              = 2*ntype*type[ii];
		
		vctot            = _mm_setzero_ps();
		Vvdwtot          = _mm_setzero_ps();
		dvdasum          = _mm_setzero_ps();
		vgbtot           = _mm_setzero_ps();
		fix              = _mm_setzero_ps();
		fiy              = _mm_setzero_ps();
		fiz              = _mm_setzero_ps();
		
		offset           = (nj1-nj0)%4;
        
        for(k=nj0;k<nj1-offset;k+=4)
		{
			jnr     = jjnr[k];   
			jnr2    = jjnr[k+1];
			jnr3    = jjnr[k+2];
			jnr4    = jjnr[k+3];
			
            j3      = 3*jnr;  
			j23		= 3*jnr2;
			j33     = 3*jnr3;
			j43     = 3*jnr4;
			
			xmm1    = _mm_loadh_pi(xmm1, (__m64 *) (pos+j3));  /* x1 y1 - - */
			xmm2    = _mm_loadh_pi(xmm2, (__m64 *) (pos+j23)); /* x2 y2 - - */
			xmm3    = _mm_loadh_pi(xmm3, (__m64 *) (pos+j33)); /* x3 y3 - - */
			xmm4    = _mm_loadh_pi(xmm4, (__m64 *) (pos+j43)); /* x4 y4 - - */			
			
			xmm5    = _mm_load1_ps(pos+j3+2);  
			xmm6    = _mm_load1_ps(pos+j23+2); 
			xmm7    = _mm_load1_ps(pos+j33+2); 
			xmm8    = _mm_load1_ps(pos+j43+2);
			
			/* transpose */
			xmm5    = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0));
			xmm6    = _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(0,0,0,0));
			jz      = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(2,0,2,0));
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2));
			xmm2    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(3,2,3,2));
			jx      = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(2,0,2,0));
			jy      = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,1,3,1));
			
			dx1     = _mm_sub_ps(ix,jx);
			dy1     = _mm_sub_ps(iy,jy);
			dz1     = _mm_sub_ps(iz,jz);
			
			t1      = _mm_mul_ps(dx1,dx1);
			t2      = _mm_mul_ps(dy1,dy1);
			t3      = _mm_mul_ps(dz1,dz1);
			
			rsq     = _mm_add_ps(t1,t2);
			rsq     = _mm_add_ps(rsq,t3);
			rinv    = gmx_mm_invsqrt_ps(rsq);
			
			xmm1    = _mm_load_ss(invsqrta+jnr); 
			xmm2    = _mm_load_ss(invsqrta+jnr2);
			xmm3    = _mm_load_ss(invsqrta+jnr3);
			xmm4    = _mm_load_ss(invsqrta+jnr4);
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); /* j1 j1 j2 j2 */
			xmm3    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); /* j3 j3 j4 j4 */
			isaj    = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			
			isaprod = _mm_mul_ps(isai,isaj);
			
			xmm1    = _mm_load_ss(charge+jnr); 
			xmm2    = _mm_load_ss(charge+jnr2); 
			xmm3    = _mm_load_ss(charge+jnr3); 
			xmm4    = _mm_load_ss(charge+jnr4);
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); /* j1 j1 j2 j2 */
			xmm3    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); /* j3 j3 j4 j4 */
			q       = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			qq      = _mm_mul_ps(iq,q);
			
			vcoul   = _mm_mul_ps(qq,rinv);
			fscal   = _mm_mul_ps(vcoul,rinv);
			qq      = _mm_mul_ps(qq,neg);
			qq      = _mm_mul_ps(isaprod,qq);
			gbscale = _mm_mul_ps(isaprod,gbtabscale_sse);
			
            tj      = nti+2*type[jnr];
			tj2     = nti+2*type[jnr2];
			tj3     = nti+2*type[jnr3];
			tj4     = nti+2*type[jnr4];
			
			xmm1    = _mm_loadl_pi(xmm1,(__m64 *) (vdwparam+tj));
			xmm1    = _mm_loadh_pi(xmm1,(__m64 *) (vdwparam+tj2));
			xmm2    = _mm_loadl_pi(xmm2,(__m64 *) (vdwparam+tj3));
			xmm2    = _mm_loadh_pi(xmm2,(__m64 *) (vdwparam+tj4));
			
			xmm3    = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(0,2,0,2)); /* c61 c62 c61 c62 */
			xmm4    = _mm_shuffle_ps(xmm2,xmm2,_MM_SHUFFLE(0,2,0,2)); /* c63 c64 c63 c64 */
			c6      = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,1,0,1)); /* c61 c62 c63 c64 */
			
			xmm3    = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,3,1)); /* c121 c122 c121 c122 */
			xmm4    = _mm_shuffle_ps(xmm2,xmm2,_MM_SHUFFLE(3,1,3,1)); /* c123 c124 c123 c124 */
			c12     = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(1,0,1,0)); /* c121 c122 c123 c124 */		
				
			xmm1    = _mm_load_ss(dvda+jnr); 
			xmm2    = _mm_load_ss(dvda+jnr2); 
			xmm3    = _mm_load_ss(dvda+jnr3);
			xmm4    = _mm_load_ss(dvda+jnr4);
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); /* j1 j1 j2 j2 */
			xmm3    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); /* j3 j3 j4 j4 */
			dvdaj   = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			
			/* Calculate GB table index */
			r       = _mm_mul_ps(rsq,rinv);
			rt      = _mm_mul_ps(r,gbscale);
			
			n0      = _mm_cvttps_epi32(rt);
			n0f     = _mm_cvtepi32_ps(n0);
			eps     = _mm_sub_ps(rt,n0f);
			
			eps2    = _mm_mul_ps(eps,eps);
			nnn     = _mm_slli_epi32(n0,2);
		
			/* the tables are 16-byte aligned, so we can use _mm_load_ps */			
			xmm1    = _mm_load_ps(GBtab+(gmx_mm_extract_epi32(nnn,0)));  /* Y1,F1,G1,H1 */
			xmm2    = _mm_load_ps(GBtab+(gmx_mm_extract_epi32(nnn,1)));  /* Y2,F2,G2,H2 */
			xmm3    = _mm_load_ps(GBtab+(gmx_mm_extract_epi32(nnn,2)));  /* Y3,F3,G3,H3 */
			xmm4    = _mm_load_ps(GBtab+(gmx_mm_extract_epi32(nnn,3)));  /* Y4,F4,G4,H4 */
			
			/* transpose 4*4 */
			xmm5    = _mm_unpacklo_ps(xmm1,xmm2); /* Y1,Y2,F1,F2 */
			xmm6    = _mm_unpacklo_ps(xmm3,xmm4); /* Y3,Y4,F3,F4 */
			xmm7    = _mm_unpackhi_ps(xmm1,xmm2); /* G1,G2,H1,H2 */
			xmm8    = _mm_unpackhi_ps(xmm3,xmm4); /* G3,G4,H3,H4 */
			
			Y       = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); /* Y1 Y2 Y3 Y4 */
			F		= _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); /* F1 F2 F3 F4 */
			G		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(1,0,1,0)); /* G1 G2 G3 G4 */
			H		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(3,2,3,2)); /* H1 H2 H3 H4 */
			
			Geps    = _mm_mul_ps(G,eps);
			Heps2   = _mm_mul_ps(H,eps2);
			Fp      = _mm_add_ps(F,Geps);
			Fp      = _mm_add_ps(Fp,Heps2);
			
			VV      = _mm_mul_ps(Fp,eps);
			VV      = _mm_add_ps(Y,VV);
			
			xmm1    = _mm_mul_ps(two,Heps2);
			FF      = _mm_add_ps(Fp,Geps);
			FF      = _mm_add_ps(FF,xmm1);
			
			vgb     = _mm_mul_ps(qq,VV);
			fijC    = _mm_mul_ps(qq,FF);
			fijC    = _mm_mul_ps(fijC,gbscale);
			
			dvdatmp = _mm_mul_ps(fijC,r);
			dvdatmp = _mm_add_ps(vgb,dvdatmp);
			dvdatmp = _mm_mul_ps(half,dvdatmp);
			dvdatmp = _mm_mul_ps(neg,dvdatmp);
			
			dvdasum = _mm_add_ps(dvdasum,dvdatmp);
			
			xmm1    = _mm_mul_ps(dvdatmp,isaj);
			xmm1    = _mm_mul_ps(xmm1,isaj);
			dvdaj   = _mm_add_ps(dvdaj,xmm1);
			
			_mm_store_ss(dvda+jnr,dvdaj);
			xmm1    = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(0,3,2,1)); 
			_mm_store_ss(dvda+jnr2,xmm1);
			xmm1    = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(1,0,3,2)); 
			_mm_store_ss(dvda+jnr3,xmm1);
			xmm1    = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(2,1,0,3));
			_mm_store_ss(dvda+jnr4,xmm1);
			
			vctot   = _mm_add_ps(vctot,vcoul);
			vgbtot  = _mm_add_ps(vgbtot,vgb);
			
			/* Calculate VDW table index */
			rt      = _mm_mul_ps(r,tabscale_sse);
			n0      = _mm_cvttps_epi32(rt);
			n0f     = _mm_cvtepi32_ps(n0);
			eps     = _mm_sub_ps(rt,n0f);
			eps2    = _mm_mul_ps(eps,eps);
			nnn     = _mm_slli_epi32(n0,3);

			/* Tabulated VdW interaction - disperion */			
			xmm1    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,0)));  /* Y1,F1,G1,H1 */
			xmm2    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,1)));  /* Y2,F2,G2,H2 */
			xmm3    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,2)));  /* Y3,F3,G3,H3 */
			xmm4    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,3)));  /* Y4,F4,G4,H4 */
			
			/* transpose 4*4 */
			xmm5    = _mm_unpacklo_ps(xmm1,xmm2); /* Y1,Y2,F1,F2 */
			xmm6    = _mm_unpacklo_ps(xmm3,xmm4); /* Y3,Y4,F3,F4 */
			xmm7    = _mm_unpackhi_ps(xmm1,xmm2); /* G1,G2,H1,H2 */
			xmm8    = _mm_unpackhi_ps(xmm3,xmm4); /* G3,G4,H3,H4 */
			
			Y       = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); /* Y1 Y2 Y3 Y4 */
			F		= _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); /* F1 F2 F3 F4 */
			G		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(1,0,1,0)); /* G1 G2 G3 G4 */
			H		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(3,2,3,2)); /* H1 H2 H3 H4 */
			
			Geps    = _mm_mul_ps(G,eps);
			Heps2   = _mm_mul_ps(H,eps2);
			Fp      = _mm_add_ps(F,Geps);
			Fp      = _mm_add_ps(Fp,Heps2);
			VV      = _mm_mul_ps(Fp,eps);
			VV      = _mm_add_ps(Y,VV);
			xmm1    = _mm_mul_ps(two,Heps2);
			FF      = _mm_add_ps(Fp,Geps);
			FF      = _mm_add_ps(FF,xmm1);

			Vvdw6   = _mm_mul_ps(c6,VV);
			fijD    = _mm_mul_ps(c6,FF);
			
			/* Tabulated VdW interaction - repulsion */
			nnn     = _mm_add_epi32(nnn,four);
			
			xmm1    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,0)));  /* Y1,F1,G1,H1 */
			xmm2    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,1)));  /* Y2,F2,G2,H2 */
			xmm3    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,2)));  /* Y3,F3,G3,H3 */
			xmm4    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,3)));  /* Y4,F4,G4,H4 */
			
			/* transpose 4*4 */
			xmm5    = _mm_unpacklo_ps(xmm1,xmm2); /* Y1,Y2,F1,F2 */
			xmm6    = _mm_unpacklo_ps(xmm3,xmm4); /* Y3,Y4,F3,F4 */
			xmm7    = _mm_unpackhi_ps(xmm1,xmm2); /* G1,G2,H1,H2 */
			xmm8    = _mm_unpackhi_ps(xmm3,xmm4); /* G3,G4,H3,H4 */
			
			Y       = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); /* Y1 Y2 Y3 Y4 */
			F		= _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); /* F1 F2 F3 F4 */
			G		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(1,0,1,0)); /* G1 G2 G3 G4 */
			H		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(3,2,3,2)); /* H1 H2 H3 H4 */
			
			Geps    = _mm_mul_ps(G,eps);
			Heps2   = _mm_mul_ps(H,eps2);
			Fp      = _mm_add_ps(F,Geps);
			Fp      = _mm_add_ps(Fp,Heps2);
			VV      = _mm_mul_ps(Fp,eps);
			VV      = _mm_add_ps(Y,VV);
			xmm1    = _mm_mul_ps(two,Heps2);
			FF      = _mm_add_ps(Fp,Geps);
			FF      = _mm_add_ps(FF,xmm1);
			
			Vvdw12  = _mm_mul_ps(c12,VV);
			fijR    = _mm_mul_ps(c12,FF);
			
			Vvdwtmp = _mm_add_ps(Vvdw12,Vvdw6);
			Vvdwtot = _mm_add_ps(Vvdwtot,Vvdwtmp);
			
			xmm1    = _mm_add_ps(fijD,fijR);
			xmm1    = _mm_mul_ps(xmm1,tabscale_sse);
			xmm1    = _mm_add_ps(xmm1,fijC);
			xmm1    = _mm_sub_ps(xmm1,fscal);
			fscal   = _mm_mul_ps(xmm1,neg);
			fscal   = _mm_mul_ps(fscal,rinv);
	
			t1      = _mm_mul_ps(fscal,dx1);
			t2      = _mm_mul_ps(fscal,dy1);
			t3      = _mm_mul_ps(fscal,dz1);
			
			fix     = _mm_add_ps(fix,t1);
			fiy     = _mm_add_ps(fiy,t2);
			fiz     = _mm_add_ps(fiz,t3);
			
			xmm1    = _mm_loadh_pi(xmm1, (__m64 *) (faction+j3));  /* fx1 fy1 - - */
			xmm2    = _mm_loadh_pi(xmm2, (__m64 *) (faction+j23)); /* fx1 fy1 - - */
			xmm3    = _mm_loadh_pi(xmm3, (__m64 *) (faction+j33)); /* fx3 fy3 - - */
			xmm4    = _mm_loadh_pi(xmm4, (__m64 *) (faction+j43)); /* fx4 fy4 - - */
			
			xmm5    = _mm_load1_ps(faction+j3+2); 
			xmm6    = _mm_load1_ps(faction+j23+2);
			xmm7    = _mm_load1_ps(faction+j33+2);
			xmm8    = _mm_load1_ps(faction+j43+2);
			
			/* transpose the forces */
			xmm5    = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0));
			xmm6    = _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(0,0,0,0));
			xmm7    = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(2,0,2,0));
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2));
			xmm2    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(3,2,3,2));
			
			xmm5    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(2,0,2,0));
			xmm6    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,1,3,1));
			
			/* subtract partial terms */
			xmm5    = _mm_sub_ps(xmm5,t1);
			xmm6    = _mm_sub_ps(xmm6,t2);
			xmm7    = _mm_sub_ps(xmm7,t3);
			
			xmm1    = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0));
			xmm2    = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2));
			xmm1    = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,2,0));
			xmm2    = _mm_shuffle_ps(xmm2,xmm2,_MM_SHUFFLE(3,1,2,0));
			
			/* store fx's and fy's */
			_mm_storel_pi( (__m64 *) (faction+j3),xmm1);
			_mm_storeh_pi( (__m64 *) (faction+j23),xmm1);
			_mm_storel_pi( (__m64 *) (faction+j33),xmm2);
			_mm_storeh_pi( (__m64 *) (faction+j43),xmm2);
			
			/* ..then z */
			_mm_store_ss(faction+j3+2,xmm7);
			xmm7    = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
			_mm_store_ss(faction+j23+2,xmm7);
			xmm7    = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
			_mm_store_ss(faction+j33+2,xmm7);
			xmm7    = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
			_mm_store_ss(faction+j43+2,xmm7);	
			
		}
        
		if(offset!=0)
		{
			if(offset==1)
			{
				jnr   = jjnr[k];
				j3    = jnr*3;
				tj	  = nti+2*type[jnr];
				
				xmm1  = _mm_loadl_pi(xmm1, (__m64 *) (pos+j3));
				xmm5  = _mm_load1_ps(pos+j3+2);
				
				xmm6  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(0,0,0,0));
				xmm4  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(1,1,1,1));
				
				isaj  = _mm_load_ss(invsqrta+jnr);
				dvdaj = _mm_load_ss(dvda+jnr);
				q     = _mm_load_ss(charge+jnr);
				c6    = _mm_load_ss(vdwparam+tj);
				c12   = _mm_load_ss(vdwparam+tj+1);
				
				mask  =  _mm_set_epi32(0,0,0,0xffffffff);
				
			}
			else if(offset==2)
			{
				jnr   = jjnr[k];
				jnr2  = jjnr[k+1];
				
				j3    = jnr*3;
				j23   = jnr2*3;
				
				tj	  = nti+2*type[jnr];
				tj2	  = nti+2*type[jnr2];
				
				xmm1  = _mm_loadh_pi(xmm1, (__m64 *) (pos+j3));
				xmm2  = _mm_loadh_pi(xmm2, (__m64 *) (pos+j23));
				
				xmm5  = _mm_load1_ps(pos+j3+2);
				xmm6  = _mm_load1_ps(pos+j23+2);
				
				xmm5  = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0));
				xmm5  = _mm_shuffle_ps(xmm5,xmm5,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2));
				xmm6  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				xmm4  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,3,1));
				
				xmm1  = _mm_load_ss(invsqrta+jnr);
				xmm2  = _mm_load_ss(invsqrta+jnr2);
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0));
				isaj  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_load_ss(dvda+jnr);
				xmm2  = _mm_load_ss(dvda+jnr2);
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0));
				dvdaj = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_load_ss(charge+jnr);
				xmm2  = _mm_load_ss(charge+jnr2);
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0));
				q     = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_loadl_pi(xmm1,(__m64 *) (vdwparam+tj));
				xmm1  = _mm_loadh_pi(xmm1,(__m64 *) (vdwparam+tj2));
				
				c6    = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				c12   = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,3,1));
				
				mask  = _mm_set_epi32(0,0,0xffffffff,0xffffffff);
			}
			else
			{
				jnr   = jjnr[k];
				jnr2  = jjnr[k+1];
				jnr3  = jjnr[k+2];
			 	
				j3    = jnr*3;
				j23   = jnr2*3;
				j33   = jnr3*3;
				
				tj	  = nti+2*type[jnr];
				tj2	  = nti+2*type[jnr2];
				tj3	  = nti+2*type[jnr3];
				
				xmm1  = _mm_loadh_pi(xmm1,(__m64 *) (pos+j3)); 
				xmm2  = _mm_loadh_pi(xmm2,(__m64 *) (pos+j23)); 
				xmm3  = _mm_loadh_pi(xmm3,(__m64 *) (pos+j33)); 
				
				xmm5  = _mm_load1_ps(pos+j3+2); 
				xmm6  = _mm_load1_ps(pos+j23+2); 
				xmm7  = _mm_load1_ps(pos+j33+2); 
				
				xmm5  = _mm_shuffle_ps(xmm5,xmm6, _MM_SHUFFLE(0,0,0,0)); 
				xmm5  = _mm_shuffle_ps(xmm5,xmm7, _MM_SHUFFLE(3,1,3,1));	
				
				xmm1  = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(3,2,3,2)); 
				xmm2  = _mm_shuffle_ps(xmm3,xmm3, _MM_SHUFFLE(3,2,3,2)); 
				
				xmm6  = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(2,0,2,0)); 
				xmm4  = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(3,1,3,1)); 
				
				xmm1  = _mm_load_ss(invsqrta+jnr);
				xmm2  = _mm_load_ss(invsqrta+jnr2);
				xmm3  = _mm_load_ss(invsqrta+jnr3);
				
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
				xmm3  = _mm_shuffle_ps(xmm3,xmm3,_MM_SHUFFLE(0,0,0,0)); 
				isaj  = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_load_ss(dvda+jnr);
				xmm2  = _mm_load_ss(dvda+jnr2);
				xmm3  = _mm_load_ss(dvda+jnr3);
				
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
				xmm3  = _mm_shuffle_ps(xmm3,xmm3,_MM_SHUFFLE(0,0,0,0)); 
				dvdaj = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_load_ss(charge+jnr);
				xmm2  = _mm_load_ss(charge+jnr2);
				xmm3  = _mm_load_ss(charge+jnr3);
				
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
				xmm3  = _mm_shuffle_ps(xmm3,xmm3,_MM_SHUFFLE(0,0,0,0)); 
				q     = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_loadl_pi(xmm1, (__m64 *) (vdwparam+tj)); 
				xmm1  = _mm_loadh_pi(xmm1, (__m64 *) (vdwparam+tj2));
				xmm2  = _mm_loadl_pi(xmm2, (__m64 *) (vdwparam+tj3));
				
				xmm3  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(0,2,0,2));
				c6    = _mm_shuffle_ps(xmm3,xmm2,_MM_SHUFFLE(1,0,0,1));
				
				xmm3  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,3,1));
				c12   = _mm_shuffle_ps(xmm3,xmm2,_MM_SHUFFLE(1,1,1,0));
				
				mask  = _mm_set_epi32(0,0xffffffff,0xffffffff,0xffffffff);
			}
			
			jx      = _mm_and_ps( gmx_mm_castsi128_ps(mask), xmm6);
			jy      = _mm_and_ps( gmx_mm_castsi128_ps(mask), xmm4);
			jz      = _mm_and_ps( gmx_mm_castsi128_ps(mask), xmm5);
			
			c6      = _mm_and_ps( gmx_mm_castsi128_ps(mask), c6);
			c12     = _mm_and_ps( gmx_mm_castsi128_ps(mask), c12);
			dvdaj   = _mm_and_ps( gmx_mm_castsi128_ps(mask), dvdaj);
			isaj    = _mm_and_ps( gmx_mm_castsi128_ps(mask), isaj);			
			q       = _mm_and_ps( gmx_mm_castsi128_ps(mask), q);
			
			dx1     = _mm_sub_ps(ix,jx);
			dy1     = _mm_sub_ps(iy,jy);
			dz1     = _mm_sub_ps(iz,jz);
			
			t1      = _mm_mul_ps(dx1,dx1);
			t2      = _mm_mul_ps(dy1,dy1);
			t3      = _mm_mul_ps(dz1,dz1);
			
			rsq     = _mm_add_ps(t1,t2);
			rsq     = _mm_add_ps(rsq,t3);
			
			rinv    = gmx_mm_invsqrt_ps(rsq);
			
			isaprod = _mm_mul_ps(isai,isaj);
			qq      = _mm_mul_ps(iq,q);
			vcoul   = _mm_mul_ps(qq,rinv);
			fscal   = _mm_mul_ps(vcoul,rinv);
			
			qq      = _mm_mul_ps(qq,neg);
			qq      = _mm_mul_ps(isaprod,qq);
			
			gbscale = _mm_mul_ps(isaprod,gbtabscale_sse);
			rinvsq  = _mm_mul_ps(rinv,rinv);
			r       = _mm_mul_ps(rsq,rinv);
			rt      = _mm_mul_ps(r,gbscale);
			
			n0      = _mm_cvttps_epi32(rt);
			n0f     = _mm_cvtepi32_ps(n0);
			eps     = _mm_sub_ps(rt,n0f);
			
			eps2    = _mm_mul_ps(eps,eps);
			nnn     = _mm_slli_epi32(n0,2);
			
			/* the tables are 16-byte aligned, so we can use _mm_load_ps */			
			xmm1    = _mm_load_ps(GBtab+(gmx_mm_extract_epi32(nnn,0)));  /* Y1,F1,G1,H1 */
			xmm2    = _mm_load_ps(GBtab+(gmx_mm_extract_epi32(nnn,1)));  /* Y2,F2,G2,H2 */
			xmm3    = _mm_load_ps(GBtab+(gmx_mm_extract_epi32(nnn,2)));  /* Y3,F3,G3,H3 */
			xmm4    = _mm_load_ps(GBtab+(gmx_mm_extract_epi32(nnn,3)));  /* Y4,F4,G4,H4 */
			
			/* transpose 4*4 */
			xmm5    = _mm_unpacklo_ps(xmm1,xmm2); /* Y1,Y2,F1,F2 */
			xmm6    = _mm_unpacklo_ps(xmm3,xmm4); /* Y3,Y4,F3,F4 */
			xmm7    = _mm_unpackhi_ps(xmm1,xmm2); /* G1,G2,H1,H2 */
			xmm8    = _mm_unpackhi_ps(xmm3,xmm4); /* G3,G4,H3,H4 */
			
			Y       = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); 
			F		= _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); 
			G		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(1,0,1,0)); 
			H		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(3,2,3,2)); 
			
			Geps = _mm_mul_ps(G,eps);
			Heps2   = _mm_mul_ps(H,eps2);
			Fp      = _mm_add_ps(F,Geps);
			Fp      = _mm_add_ps(Fp,Heps2);
			VV      = _mm_mul_ps(Fp,eps);
			VV      = _mm_add_ps(Y,VV);
			
			xmm1    = _mm_mul_ps(two,Heps2);
			FF      = _mm_add_ps(Fp,Geps);
			FF      = _mm_add_ps(FF,xmm1);
			vgb     = _mm_mul_ps(qq,VV);
			fijC    = _mm_mul_ps(qq,FF);
			fijC    = _mm_mul_ps(fijC,gbscale);
			
			dvdatmp = _mm_mul_ps(fijC,r);
			dvdatmp = _mm_add_ps(vgb,dvdatmp);
			dvdatmp = _mm_mul_ps(half,dvdatmp);
			dvdatmp = _mm_mul_ps(neg,dvdatmp);
			
			dvdasum = _mm_add_ps(dvdasum,dvdatmp);
			
			xmm1    = _mm_mul_ps(dvdatmp,isaj);
			xmm1    = _mm_mul_ps(xmm1,isaj);
			dvdaj   = _mm_add_ps(dvdaj,xmm1);
			
			vcoul   = _mm_and_ps( gmx_mm_castsi128_ps(mask), vcoul);
			vgb     = _mm_and_ps( gmx_mm_castsi128_ps(mask), vgb);
			
			vctot   = _mm_add_ps(vctot,vcoul);
			vgbtot  = _mm_add_ps(vgbtot,vgb);
						
			/* Calculate VDW table index */
			rt      = _mm_mul_ps(r,tabscale_sse);
			n0      = _mm_cvttps_epi32(rt);
			n0f     = _mm_cvtepi32_ps(n0);
			eps     = _mm_sub_ps(rt,n0f);
			eps2    = _mm_mul_ps(eps,eps);
			nnn     = _mm_slli_epi32(n0,3);
			
			/* Tabulated VdW interaction - disperion */	
			xmm1    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,0)));  /* Y1,F1,G1,H1 */
			xmm2    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,1)));  /* Y2,F2,G2,H2 */
			xmm3    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,2)));  /* Y3,F3,G3,H3 */
			xmm4    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,3)));  /* Y4,F4,G4,H4 */
		
			/* transpose 4*4 */
			xmm5    = _mm_unpacklo_ps(xmm1,xmm2); /* Y1,Y2,F1,F2 */
			xmm6    = _mm_unpacklo_ps(xmm3,xmm4); /* Y3,Y4,F3,F4 */
			xmm7    = _mm_unpackhi_ps(xmm1,xmm2); /* G1,G2,H1,H2 */
			xmm8    = _mm_unpackhi_ps(xmm3,xmm4); /* G3,G4,H3,H4 */
			
			Y       = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); /* Y1 Y2 Y3 Y4 */
			F		= _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); /* F1 F2 F3 F4 */
			G		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(1,0,1,0)); /* G1 G2 G3 G4 */
			H		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(3,2,3,2)); /* H1 H2 H3 H4 */
			
			Geps    = _mm_mul_ps(G,eps);
			Heps2   = _mm_mul_ps(H,eps2);
			Fp      = _mm_add_ps(F,Geps);
			Fp      = _mm_add_ps(Fp,Heps2);
			VV      = _mm_mul_ps(Fp,eps);
			VV      = _mm_add_ps(Y,VV);
			xmm1    = _mm_mul_ps(two,Heps2);
			FF      = _mm_add_ps(Fp,Geps);
			FF      = _mm_add_ps(FF,xmm1);
			
			Vvdw6   = _mm_mul_ps(c6,VV);
			fijD    = _mm_mul_ps(c6,FF);
			
			/* Tabulated VdW interaction - repulsion */
			nnn     = _mm_add_epi32(nnn,four);
					
			xmm1    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,0)));  /* Y1,F1,G1,H1 */
			xmm2    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,1)));  /* Y2,F2,G2,H2 */
			xmm3    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,2)));  /* Y3,F3,G3,H3 */
			xmm4    = _mm_load_ps(VFtab+(gmx_mm_extract_epi32(nnn,3)));  /* Y4,F4,G4,H4 */
			
			/* transpose 4*4 */
			xmm5    = _mm_unpacklo_ps(xmm1,xmm2); /* Y1,Y2,F1,F2 */
			xmm6    = _mm_unpacklo_ps(xmm3,xmm4); /* Y3,Y4,F3,F4 */
			xmm7    = _mm_unpackhi_ps(xmm1,xmm2); /* G1,G2,H1,H2 */
			xmm8    = _mm_unpackhi_ps(xmm3,xmm4); /* G3,G4,H3,H4 */
			
			Y       = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); /* Y1 Y2 Y3 Y4 */
			F		= _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); /* F1 F2 F3 F4 */
			G		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(1,0,1,0)); /* G1 G2 G3 G4 */
			H		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(3,2,3,2)); /* H1 H2 H3 H4 */
			
			Geps    = _mm_mul_ps(G,eps);
			Heps2   = _mm_mul_ps(H,eps2);
			Fp      = _mm_add_ps(F,Geps);
			Fp      = _mm_add_ps(Fp,Heps2);
			VV      = _mm_mul_ps(Fp,eps);
			VV      = _mm_add_ps(Y,VV);
			xmm1    = _mm_mul_ps(two,Heps2);
			FF      = _mm_add_ps(Fp,Geps);
			FF      = _mm_add_ps(FF,xmm1);
			
			Vvdw12  = _mm_mul_ps(c12,VV);
			fijR    = _mm_mul_ps(c12,FF);
			
			Vvdwtmp = _mm_add_ps(Vvdw12,Vvdw6);
			Vvdwtot = _mm_add_ps(Vvdwtot,Vvdwtmp);
			
			xmm1    = _mm_add_ps(fijD,fijR);
			xmm1    = _mm_mul_ps(xmm1,tabscale_sse);
			xmm1    = _mm_add_ps(xmm1,fijC);
			xmm1    = _mm_sub_ps(xmm1,fscal);
			fscal   = _mm_mul_ps(xmm1,neg);
			fscal   = _mm_mul_ps(fscal,rinv);			
			
			t1      = _mm_mul_ps(fscal,dx1);
			t2      = _mm_mul_ps(fscal,dy1);
			t3      = _mm_mul_ps(fscal,dz1);
			
			if(offset==1)
			{
				_mm_store_ss(dvda+jnr,dvdaj);
				
				xmm1 = _mm_loadl_pi(xmm1, (__m64 *) (faction+j3));
				xmm7 = _mm_load1_ps(faction+j3+2);
				
				xmm5 = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(0,0,0,0));
				xmm6 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(1,1,1,1));
				
				xmm5 = _mm_sub_ps(xmm5,t1);
				xmm6 = _mm_sub_ps(xmm6,t2);
				xmm7 = _mm_sub_ps(xmm7,t3);
				
				_mm_store_ss(faction+j3 , xmm5);
				_mm_store_ss(faction+j3+1,xmm6);
				_mm_store_ss(faction+j3+2,xmm7);
			}
			else if(offset==2)
			{
				_mm_store_ss(dvda+jnr,dvdaj);
				xmm1 = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(0,3,2,1)); 
				_mm_store_ss(dvda+jnr2,xmm1);
				
				xmm1 = _mm_loadh_pi(xmm1, (__m64 *) (faction+j3));
				xmm2 = _mm_loadh_pi(xmm2, (__m64 *) (faction+j23)); 
				
				xmm5 = _mm_load1_ps(faction+j3+2); 
				xmm6 = _mm_load1_ps(faction+j23+2);
				
				xmm5 = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0)); 
				xmm7 = _mm_shuffle_ps(xmm5,xmm5,_MM_SHUFFLE(2,0,2,0)); 
				
				xmm1 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2)); 
				xmm5 = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0)); 
				xmm6 = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,3,1)); 
				
				xmm5 = _mm_sub_ps(xmm5, t1);
				xmm6 = _mm_sub_ps(xmm6, t2);
				xmm7 = _mm_sub_ps(xmm7, t3);
				
				xmm1 = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); 
				xmm5 = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,2,0)); 
				
				_mm_storel_pi( (__m64 *) (faction+j3), xmm5);
				_mm_storeh_pi( (__m64 *) (faction+j23), xmm5);
				
				_mm_store_ss(faction+j3+2,xmm7);
				xmm7 = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
				_mm_store_ss(faction+j23+2,xmm7);
			}
			else
			{
				_mm_store_ss(dvda+jnr,dvdaj);
				xmm1 = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(0,3,2,1)); 
				_mm_store_ss(dvda+jnr2,xmm1);
				xmm1 = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(1,0,3,2)); 
				_mm_store_ss(dvda+jnr3,xmm1);
				
				xmm1 = _mm_loadh_pi(xmm1, (__m64 *) (faction+j3)); 
				xmm2 = _mm_loadh_pi(xmm2, (__m64 *) (faction+j23));  
				xmm3 = _mm_loadh_pi(xmm3, (__m64 *) (faction+j33));
				
				xmm5 = _mm_load1_ps(faction+j3+2);
				xmm6 = _mm_load1_ps(faction+j23+2); 
				xmm7 = _mm_load1_ps(faction+j33+2); 
				
				xmm5 = _mm_shuffle_ps(xmm5,xmm6, _MM_SHUFFLE(0,0,0,0)); 
				xmm6 = _mm_shuffle_ps(xmm7,xmm7, _MM_SHUFFLE(0,0,0,0)); 
				xmm7 = _mm_shuffle_ps(xmm5,xmm6, _MM_SHUFFLE(2,0,2,0)); 
				
				xmm1 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2)); 
				xmm2 = _mm_shuffle_ps(xmm3,xmm3,_MM_SHUFFLE(3,2,3,2)); 
				
				xmm5 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(2,0,2,0)); 
				xmm6 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,1,3,1)); 
				
				xmm5 = _mm_sub_ps(xmm5, t1);
				xmm6 = _mm_sub_ps(xmm6, t2);
				xmm7 = _mm_sub_ps(xmm7, t3);
				
				xmm1 = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); 
				xmm2 = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); 
				xmm1 = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,2,0)); 
				xmm2 = _mm_shuffle_ps(xmm2,xmm2,_MM_SHUFFLE(3,1,2,0)); 
				
				_mm_storel_pi( (__m64 *) (faction+j3), xmm1);
				_mm_storeh_pi( (__m64 *) (faction+j23), xmm1);
				_mm_storel_pi( (__m64 *) (faction+j33), xmm2);
				
				_mm_store_ss(faction+j3+2,xmm7); 
				xmm7 = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
				_mm_store_ss(faction+j23+2,xmm7); 
				xmm7 = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
				_mm_store_ss(faction+j33+2,xmm7); 
			}
			
			t1 = _mm_and_ps( gmx_mm_castsi128_ps(mask), t1);
			t2 = _mm_and_ps( gmx_mm_castsi128_ps(mask), t2);
			t3 = _mm_and_ps( gmx_mm_castsi128_ps(mask), t3);
			
			fix = _mm_add_ps(fix,t1);
			fiy = _mm_add_ps(fiy,t2);
			fiz = _mm_add_ps(fiz,t3);
		}
		
		t1      = _mm_movehl_ps(t1,fix);
		t2      = _mm_movehl_ps(t2,fiy);
		t3      = _mm_movehl_ps(t3,fiz);
		
		fix     = _mm_add_ps(fix,t1);
		fiy     = _mm_add_ps(fiy,t2);
		fiz     = _mm_add_ps(fiz,t3);
		
		t1      = _mm_shuffle_ps(fix,fix,_MM_SHUFFLE(1,1,1,1));
		t2      = _mm_shuffle_ps(fiy,fiy,_MM_SHUFFLE(1,1,1,1));
		t3      = _mm_shuffle_ps(fiz,fiz,_MM_SHUFFLE(1,1,1,1));
		
		fix     = _mm_add_ps(fix,t1);
		fiy     = _mm_add_ps(fiy,t2);
		fiz     = _mm_add_ps(fiz,t3);
		
		xmm2    = _mm_unpacklo_ps(fix,fiy); /* fx, fy, - - */
		xmm2    = _mm_movelh_ps(xmm2,fiz); 
		xmm2    = _mm_and_ps( gmx_mm_castsi128_ps(maski), xmm2);
		
		/* load i force from memory */
		xmm4    = _mm_loadl_pi(xmm4, (__m64 *) (faction+ii3));
		xmm5    = _mm_load1_ps(faction+ii3+2);
		xmm4    = _mm_shuffle_ps(xmm4,xmm5,_MM_SHUFFLE(3,2,1,0));
		
		/* add to i force */
		xmm4    = _mm_add_ps(xmm4,xmm2);
		
		/* store i force to memory */
		_mm_storel_pi( (__m64 *) (faction+ii3),xmm4);
		xmm4    = _mm_shuffle_ps(xmm4,xmm4,_MM_SHUFFLE(2,2,2,2));
		_mm_store_ss(faction+ii3+2,xmm4);
		
    /* Load, add and store i shift forces */
    xmm4 = _mm_loadl_pi(xmm4, (__m64 *) (fshift+is3));
    xmm5 = _mm_load1_ps(fshift+is3+2);
    xmm4 = _mm_shuffle_ps(xmm4,xmm5,_MM_SHUFFLE(3,2,1,0));
      
    xmm4 = _mm_add_ps(xmm4,xmm2);
      
    _mm_storel_pi( (__m64 *) (fshift+is3),xmm4);
    xmm4 = _mm_shuffle_ps(xmm4,xmm4,_MM_SHUFFLE(2,2,2,2));
    _mm_store_ss(fshift+is3+2,xmm4);
      
    /* Coulomb potential */
    ggid             = gid[n];         
		
		vcoul   = _mm_movehl_ps(vcoul,vctot);
		vctot   = _mm_add_ps(vctot,vcoul);
		vcoul   = _mm_shuffle_ps(vctot,vctot,_MM_SHUFFLE(1,1,1,1));
		vctot   = _mm_add_ss(vctot,vcoul);
		
		_mm_store_ss(&vct,vctot);
		Vc[ggid] = Vc[ggid] + vct;
		
		Vvdwtmp  = _mm_movehl_ps(Vvdwtmp,Vvdwtot);
		Vvdwtot  = _mm_add_ps(Vvdwtot,Vvdwtmp);
		Vvdwtmp  = _mm_shuffle_ps(Vvdwtot,Vvdwtot,_MM_SHUFFLE(1,1,1,1));
		Vvdwtot  = _mm_add_ss(Vvdwtot,Vvdwtmp);
		
		_mm_store_ss(&vdwt,Vvdwtot);
		Vvdw[ggid] = Vvdw[ggid] + vdwt;
		
		/* dvda */
		dvdatmp = _mm_movehl_ps(dvdatmp,dvdasum);
		dvdasum = _mm_add_ps(dvdasum,dvdatmp);
		dvdatmp = _mm_shuffle_ps(dvdasum,dvdasum,_MM_SHUFFLE(1,1,1,1));
		dvdasum = _mm_add_ss(dvdasum,dvdatmp);
		
		_mm_store_ss(&dva,dvdasum);
		dvda[ii] = dvda[ii] + dva*isai_f*isai_f;
		
		/* Store the GB potential to the work array */
		vgb     = _mm_movehl_ps(vgb,vgbtot);
		vgbtot  = _mm_add_ps(vgbtot,vgb);
		vgb     = _mm_shuffle_ps(vgbtot,vgbtot,_MM_SHUFFLE(1,1,1,1));
		vgbtot  = _mm_add_ss(vgbtot,vgb);
		
		_mm_store_ss(&vgbt,vgbtot);
		gpol[ggid] = gpol[ggid] + vgbt;
    }
	
	*outeriter       = nri;            
    *inneriter       = nj1;            
}


/*
 * Gromacs nonbonded kernel nb_kernel430nf
 * Coulomb interaction:     Generalized-Born
 * VdW interaction:         Tabulated
 * water optimization:      No
 * Calculate forces:        no
 */
void nb_kernel430nf_x86_64_sse(
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
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale,vgb,fgb;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    float         shX,shY,shZ;
    float         iq;
    float         qq,vcoul,vctot;
    int           nti;
    int           tj;
    float         Vvdw6,Vvdwtot;
    float         Vvdw12;
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         isai,isaj,isaprod,gbscale;
    float         ix1,iy1,iz1;
    float         jx1,jy1,jz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         c6,c12;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    gbtabscale       = *p_gbtabscale;  
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
        ix1              = shX + pos[ii3+0];
        iy1              = shY + pos[ii3+1];
        iz1              = shZ + pos[ii3+2];
        iq               = facel*charge[ii];
        isai             = invsqrta[ii];   
        nti              = 2*ntype*type[ii];
        vctot            = 0;              
        Vvdwtot          = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            jx1              = pos[j3+0];      
            jy1              = pos[j3+1];      
            jz1              = pos[j3+2];      
            dx11             = ix1 - jx1;      
            dy11             = iy1 - jy1;      
            dz11             = iz1 - jz1;      
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
            rinv11           = gmx_invsqrt(rsq11);
            isaj             = invsqrta[jnr];  
            isaprod          = isai*isaj;      
            qq               = iq*charge[jnr]; 
            vcoul            = qq*rinv11;      
            qq               = isaprod*(-qq);  
            gbscale          = isaprod*gbtabscale;
            tj               = nti+2*type[jnr];
            c6               = vdwparam[tj];   
            c12              = vdwparam[tj+1]; 
            r                = rsq11*rinv11;   
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
            vgb              = qq*VV;          
            vctot            = vctot + vcoul;  
            r                = rsq11*rinv11;   
            rt               = r*tabscale;     
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 8*n0;           
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            Vvdw6            = c6*VV;          
            nnn              = nnn+4;          
            Y                = VFtab[nnn];     
            F                = VFtab[nnn+1];   
            Geps             = eps*VFtab[nnn+2];
            Heps2            = eps2*VFtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            Vvdw12           = c12*VV;         
            Vvdwtot          = Vvdwtot+ Vvdw6 + Vvdw12;
        }
        
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
        Vvdw[ggid]       = Vvdw[ggid] + Vvdwtot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}


