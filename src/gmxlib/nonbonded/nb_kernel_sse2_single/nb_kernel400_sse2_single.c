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

/* to extract single integers from a __m128i datatype */
#define _mm_extract_epi32(x, imm) \
    _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

static inline __m128
my_invrsq_ps(__m128 x)
{
	const __m128 three = (const __m128) {3.0f, 3.0f, 3.0f, 3.0f};
	const __m128 half  = (const __m128) {0.5f, 0.5f, 0.5f, 0.5f};
	
	__m128 t1 = _mm_rsqrt_ps(x);
	
	return (__m128) _mm_mul_ps(half,_mm_mul_ps(t1,_mm_sub_ps(three,_mm_mul_ps(x,_mm_mul_ps(t1,t1)))));
}

void nb_kernel400_sse2_single(int *           p_nri,
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
	float         krf,facel,crf,tabscl,gbtabscl,vct,vgbt;
	int           n,ii,is3,ii3,k,nj0,nj1,jnr1,jnr2,jnr3,jnr4,j13,j23,j33,j43,ggid;
	float         shX,shY,shZ,isai_f,dva;
	
	__m128        ix,iy,iz,jx,jy,jz;
	__m128		  dx,dy,dz,t1,t2,t3;
	__m128		  fix,fiy,fiz,rsq11,rinv,r,fscal,rt,eps,eps2;
	__m128		  q,iq,qq,isai,isaj,isaprod,vcoul,gbscale,dvdai,dvdaj;
	__m128        Y,F,G,H,Fp,VV,FF,vgb,fijC,dvdatmp,dvdasum,vctot,n0f,vgbtot;
	__m128        fac,tabscale,gbtabscale;
	__m128		  xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8;
	__m128i       n0,nnn;
	
	const __m128 neg    = {-1.0f,-1.0f,-1.0f,-1.0f};
	const __m128 zero   = {0.0f,0.0f,0.0f,0.0f};
	const __m128 half   = {0.5f,0.5f,0.5f,0.5f};
	const __m128 two    = {2.0f,2.0f,2.0f,2.0f};
	const __m128 three  = {3.0f,3.0f,3.0f,3.0f};
	
	__m128i mask        = _mm_set_epi32(0, 0xffffffff,0xffffffff,0xffffffff);
	__m128i maski       = _mm_set_epi32(0, 0xffffffff, 0xffffffff, 0xffffffff);
	
	nri        = *p_nri;
	ntype      = *p_ntype;
	nthreads   = *p_nthreads; 
	facel      = *p_facel;
	krf        = *p_krf;
	crf        = *p_crf;
	tabscl     = *p_tabscale;
	gbtabscl   = *p_gbtabscale;
	nj1        = 0;
	
	/* Splat variables */
	fac        = _mm_load1_ps(&facel);
	tabscale   = _mm_load1_ps(&tabscl);
	gbtabscale = _mm_load1_ps(&gbtabscl);
	
	/* Keep the compiler happy */
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
	
	for(n=0;n<nri;n++)
	{
		is3     = 3*shift[n];
		shX     = shiftvec[is3];
		shY     = shiftvec[is3+1];
		shZ     = shiftvec[is3+2];
		nj0     = jindex[n];      
        nj1     = jindex[n+1];  
		offset  = (nj1-nj0)%4;
		ii      = iinr[n];
		ii3     = ii*3;
		ix      = _mm_set1_ps(shX+pos[ii3+0]);
		iy      = _mm_set1_ps(shX+pos[ii3+1]);
		iz      = _mm_set1_ps(shX+pos[ii3+2]); 
		q       = _mm_set1_ps(charge[ii]);
		iq      = _mm_mul_ps(fac,q); 
		isai_f  = invsqrta[ii];
		isai    = _mm_load1_ps(&isai_f);
		
		fix     = _mm_setzero_ps();
		fiy     = _mm_setzero_ps();
		fiz     = _mm_setzero_ps();
		dvdasum = _mm_setzero_ps();
		vctot   = _mm_setzero_ps();
		vgbtot  = _mm_setzero_ps();
		
		for(k=nj0;k<nj1-offset; k+=4)
		{
			jnr1    = jjnr[k];
			jnr2    = jjnr[k+1];
			jnr3    = jjnr[k+2];
			jnr4    = jjnr[k+3];
			
			j13     = jnr1 * 3;
			j23     = jnr2 * 3;
			j33     = jnr3 * 3;
			j43     = jnr4 * 3;
			
			/* Load coordinates */
			xmm1    = _mm_loadh_pi(xmm1, (__m64 *) (pos+j13)); /* x1 y1 - - */ 
			xmm2    = _mm_loadh_pi(xmm2, (__m64 *) (pos+j23)); /* x2 y2 - - */
			xmm3    = _mm_loadh_pi(xmm3, (__m64 *) (pos+j33)); /* x3 y3 - - */
			xmm4    = _mm_loadh_pi(xmm4, (__m64 *) (pos+j43)); /* x4 y4 - - */			
			
			xmm5    = _mm_load1_ps(pos+j13+2);  
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
			
			/* distances */																
			dx      = _mm_sub_ps(ix,jx);
			dy      = _mm_sub_ps(iy,jy);
			dz      = _mm_sub_ps(iz,jz);
			
			t1      = _mm_mul_ps(dx,dx);
			t2      = _mm_mul_ps(dy,dy);
			t3      = _mm_mul_ps(dz,dz);
			
			rsq11   = _mm_add_ps(t1,t2);
			rsq11   = _mm_add_ps(rsq11,t3);
			
			rinv    = my_invrsq_ps(rsq11);
			
			xmm1    = _mm_load_ss(invsqrta+jnr1); 
			xmm2    = _mm_load_ss(invsqrta+jnr2); 
			xmm3    = _mm_load_ss(invsqrta+jnr3); 
			xmm4    = _mm_load_ss(invsqrta+jnr4);
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
			xmm3    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); 
			isaj    = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			
			isaprod = _mm_mul_ps(isai,isaj);
			
			xmm1    = _mm_load_ss(charge+jnr1); 
			xmm2    = _mm_load_ss(charge+jnr2);  
			xmm3    = _mm_load_ss(charge+jnr3); 
			xmm4    = _mm_load_ss(charge+jnr4);
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
			xmm3    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); 
			q       = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			
			qq      = _mm_mul_ps(iq,q);
			vcoul   = _mm_mul_ps(qq,rinv);
			fscal   = _mm_mul_ps(vcoul,rinv);
			qq      = _mm_mul_ps(isaprod, qq);
			qq      = _mm_mul_ps(qq,neg);
			gbscale = _mm_mul_ps(isaprod,gbtabscale);
			
			/* load dvdaj */
			xmm1    = _mm_load_ss(dvda+jnr1); 
			xmm2    = _mm_load_ss(dvda+jnr2);  
			xmm3    = _mm_load_ss(dvda+jnr3); 
			xmm4    = _mm_load_ss(dvda+jnr4);
			
			xmm1    = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
			xmm3    = _mm_shuffle_ps(xmm3,xmm4,_MM_SHUFFLE(0,0,0,0)); 
			dvdaj   = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
			
			r       = _mm_mul_ps(rsq11,rinv);
			rt      = _mm_mul_ps(r,gbscale);
			n0      = _mm_cvttps_epi32(rt); 
			n0f     = _mm_cvtepi32_ps(n0); 
			eps     = _mm_sub_ps(rt,n0f);
			eps2    = _mm_mul_ps(eps,eps);
			
			nnn     = _mm_slli_epi32(n0,2); 
			
			/* the tables are 16-byte aligned, so we can use _mm_load_ps */			
			xmm1    = _mm_load_ps(GBtab+(_mm_extract_epi32(nnn,0)));  /* Y1,F1,G1,H1 */
			xmm2    = _mm_load_ps(GBtab+(_mm_extract_epi32(nnn,1)));  /* Y2,F2,G2,H2 */ 
			xmm3    = _mm_load_ps(GBtab+(_mm_extract_epi32(nnn,2)));  /* Y3,F3,G3,H3 */
			xmm4    = _mm_load_ps(GBtab+(_mm_extract_epi32(nnn,3)));  /* Y4,F4,G4,H4 */
			
			/* transpose 4*4 */
			xmm5    = _mm_unpacklo_ps(xmm1,xmm2); /* Y1,Y2,F1,F2 */
			xmm6    = _mm_unpacklo_ps(xmm3,xmm4); /* Y3,Y4,F3,F4 */
			xmm7    = _mm_unpackhi_ps(xmm1,xmm2); /* G1,G2,H1,H2 */
			xmm8    = _mm_unpackhi_ps(xmm3,xmm4); /* G3,G4,H3,H4 */
			
			Y       = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); /* Y1 Y2 Y3 Y4 */
			F		= _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); /* F1 F2 F3 F4 */
			G		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(1,0,1,0)); /* G1 G2 G3 G4 */
			H		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(3,2,3,2)); /* H1 H2 H3 H4 */
			
			G       = _mm_mul_ps(G,eps);
			H       = _mm_mul_ps(H,eps2); 
			Fp      = _mm_add_ps(F,G);
			Fp      = _mm_add_ps(Fp,H); 
			VV      = _mm_mul_ps(Fp,eps);
			VV      = _mm_add_ps(Y,VV); 
			xmm1    = _mm_mul_ps(two,H);
			FF      = _mm_add_ps(Fp,G);
			FF      = _mm_add_ps(FF,xmm1); 
			vgb     = _mm_mul_ps(qq,VV);
			fijC    = _mm_mul_ps(qq,FF);
			fijC    = _mm_mul_ps(fijC,gbscale); 
			
			dvdatmp = _mm_mul_ps(fijC,r);
			dvdatmp = _mm_add_ps(vgb,dvdatmp);
			dvdatmp = _mm_mul_ps(neg,dvdatmp);
			dvdatmp = _mm_mul_ps(dvdatmp,half);
			dvdasum = _mm_add_ps(dvdasum,dvdatmp);
			
			xmm1    = _mm_mul_ps(dvdatmp,isaj);
			xmm1    = _mm_mul_ps(xmm1,isaj);
			dvdaj   = _mm_add_ps(dvdaj,xmm1);
			
			/* store dvdaj */
			_mm_store_ss(dvda+jnr1,dvdaj);
			xmm1    = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(0,3,2,1)); 
			_mm_store_ss(dvda+jnr2,xmm1);
			xmm1    = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(1,0,3,2)); 
			_mm_store_ss(dvda+jnr3,xmm1);
			xmm1    = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(2,1,0,3));
			_mm_store_ss(dvda+jnr4,xmm1);
			
			vctot   = _mm_add_ps(vctot,vcoul);
			vgbtot  = _mm_add_ps(vgbtot,vgb);
			
			fscal   = _mm_sub_ps(fijC,fscal);
			fscal   = _mm_mul_ps(neg,fscal);
			fscal   = _mm_mul_ps(fscal,rinv);
			
			/* calculate partial force terms */			
			t1      = _mm_mul_ps(fscal,dx); /* fx1,fx2,fx3,fx4 */
			t2      = _mm_mul_ps(fscal,dy); /* fy1,fy2,fy3,fy4 */
			t3      = _mm_mul_ps(fscal,dz); /* fz1,fz2,fz3,fz4 */
			
			/* update the i force */
			fix     = _mm_add_ps(fix,t1);
			fiy     = _mm_add_ps(fiy,t2);
			fiz     = _mm_add_ps(fiz,t3);
			
			/* accumulate fx's and fy's from memory */
			xmm1    = _mm_loadh_pi(xmm1, (__m64 *) (faction+j13)); /* fx1 fy1 - - */
			xmm2    = _mm_loadh_pi(xmm2, (__m64 *) (faction+j23)); /* fx2 fy2 - - */
			xmm3    = _mm_loadh_pi(xmm3, (__m64 *) (faction+j33)); /* fx3 fy3 - - */
			xmm4    = _mm_loadh_pi(xmm4, (__m64 *) (faction+j43)); /* fx4 fy4 - - */
			
			xmm5    = _mm_load1_ps(faction+j13+2); 
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
			_mm_storel_pi( (__m64 *) (faction+j13),xmm1);
			_mm_storeh_pi( (__m64 *) (faction+j23),xmm1);
			_mm_storel_pi( (__m64 *) (faction+j33),xmm2);
			_mm_storeh_pi( (__m64 *) (faction+j43),xmm2);
			
			/* ..then z */
			_mm_store_ss(faction+j13+2,xmm7);
			xmm7    = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
			_mm_store_ss(faction+j23+2,xmm7);
			xmm7    = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
			_mm_store_ss(faction+j33+2,xmm7);
			xmm7    = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
			_mm_store_ss(faction+j43+2,xmm7);
		}
		
		if(offset!=0)
		{
			jnr1=jnr2=jnr3=jnr4=0;
			j13=j23=j33=j43=0;
			
			if(offset==1)
			{
				jnr1  = jjnr[k];
				j13   = jnr1*3;
				
				xmm1  = _mm_loadl_pi(xmm1, (__m64 *) (pos+j13));
				xmm5  = _mm_load1_ps(pos+j13+2);
				
				xmm6  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(0,0,0,0));
				xmm4  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(1,1,1,1));
				
				isaj  = _mm_load_ss(invsqrta+jnr1);
				dvdaj = _mm_load_ss(dvda+jnr1);
				q     = _mm_load_ss(charge+jnr1);
				
				mask  =  _mm_set_epi32(0,0,0,0xffffffff);
			}
			else if(offset==2)
			{
				jnr1  = jjnr[k];
				jnr2  = jjnr[k+1];
				
				j13   = jnr1*3;
				j23   = jnr2*3;
				
				xmm1  = _mm_loadh_pi(xmm1, (__m64 *) (pos+j13));
				xmm2  = _mm_loadh_pi(xmm2, (__m64 *) (pos+j23));
				
				xmm5  = _mm_load1_ps(pos+j13+2);
				xmm6  = _mm_load1_ps(pos+j23+2);
				
				xmm5  = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(0,0,0,0));
				xmm5  = _mm_shuffle_ps(xmm5,xmm5,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(3,2,3,2));
				xmm6  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				xmm4  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(3,1,3,1));
				
				xmm1  = _mm_load_ss(invsqrta+jnr1);
				xmm2  = _mm_load_ss(invsqrta+jnr2);
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0));
				isaj  = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_load_ss(dvda+jnr1);
				xmm2  = _mm_load_ss(dvda+jnr2);
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0));
				dvdaj = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				
				xmm1 = _mm_load_ss(charge+jnr1);
				xmm2 = _mm_load_ss(charge+jnr2);
				xmm1 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0));
				q    = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(2,0,2,0));
				
				mask  = _mm_set_epi32(0,0,0xffffffff,0xffffffff);
			}
			else
			{
				jnr1  = jjnr[k];
				jnr2  = jjnr[k+1];
				jnr3  = jjnr[k+2];
			 	
				j13   = jnr1*3;
				j23   = jnr2*3;
				j33   = jnr3*3;
				
				xmm1 = _mm_loadh_pi(xmm1,(__m64 *) (pos+j13)); 
				xmm2 = _mm_loadh_pi(xmm2,(__m64 *) (pos+j23)); 
				xmm3 = _mm_loadh_pi(xmm3,(__m64 *) (pos+j33)); 
				
				xmm5 = _mm_load1_ps(pos+j13+2); 
				xmm6 = _mm_load1_ps(pos+j23+2); 
				xmm7 = _mm_load1_ps(pos+j33+2); 
				
				xmm5 = _mm_shuffle_ps(xmm5,xmm6, _MM_SHUFFLE(0,0,0,0));
				xmm5 = _mm_shuffle_ps(xmm5,xmm7, _MM_SHUFFLE(3,1,3,1));						
				
				xmm1 = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(3,2,3,2));
				xmm2 = _mm_shuffle_ps(xmm3,xmm3, _MM_SHUFFLE(3,2,3,2));
				
				xmm6 = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(2,0,2,0)); 
				xmm4 = _mm_shuffle_ps(xmm1,xmm2, _MM_SHUFFLE(3,1,3,1));
				
				xmm1  = _mm_load_ss(invsqrta+jnr1);
				xmm2  = _mm_load_ss(invsqrta+jnr2);
				xmm3  = _mm_load_ss(invsqrta+jnr3);
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
				xmm3  = _mm_shuffle_ps(xmm3,xmm3,_MM_SHUFFLE(0,0,0,0)); 
				isaj  = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_load_ss(dvda+jnr1);
				xmm2  = _mm_load_ss(dvda+jnr2);
				xmm3  = _mm_load_ss(dvda+jnr3);
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
				xmm3  = _mm_shuffle_ps(xmm3,xmm3,_MM_SHUFFLE(0,0,0,0)); 
				dvdaj = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
				
				xmm1  = _mm_load_ss(charge+jnr1);
				xmm2  = _mm_load_ss(charge+jnr2);
				xmm3  = _mm_load_ss(charge+jnr3);
				xmm1  = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(0,0,0,0)); 
				xmm3  = _mm_shuffle_ps(xmm3,xmm3,_MM_SHUFFLE(0,0,0,0)); 
				q     = _mm_shuffle_ps(xmm1,xmm3,_MM_SHUFFLE(2,0,2,0));
				
				mask  = _mm_set_epi32(0,0xffffffff,0xffffffff,0xffffffff);
			}	
			
			jx      = _mm_and_ps( (__m128) mask, xmm6);
			jy      = _mm_and_ps( (__m128) mask, xmm4);
			jz      = _mm_and_ps( (__m128) mask, xmm5);
			
			dvdaj   = _mm_and_ps( (__m128) mask, dvdaj);
			isaj    = _mm_and_ps( (__m128) mask, isaj);			
			q       = _mm_and_ps( (__m128) mask, q);
			
			dx      = _mm_sub_ps(ix,jx);
			dy      = _mm_sub_ps(iy,jy);
			dz      = _mm_sub_ps(iz,jz);
			
			t1      = _mm_mul_ps(dx,dx);
			t2      = _mm_mul_ps(dy,dy);
			t3      = _mm_mul_ps(dz,dz);
			
			rsq11   = _mm_add_ps(t1,t2);
			rsq11   = _mm_add_ps(rsq11,t3);
			
			rinv    = my_invrsq_ps(rsq11);
			
			isaprod = _mm_mul_ps(isai,isaj);
			qq      = _mm_mul_ps(iq,q);
			vcoul   = _mm_mul_ps(qq,rinv);
			fscal   = _mm_mul_ps(vcoul,rinv);
			
			qq      = _mm_mul_ps(isaprod, qq);
			qq      = _mm_mul_ps(qq,neg);
			gbscale = _mm_mul_ps(isaprod,gbtabscale);
			
			r       = _mm_mul_ps(rsq11,rinv);
			rt      = _mm_mul_ps(r,gbscale);
			n0      = _mm_cvttps_epi32(rt); 
			n0f     = _mm_cvtepi32_ps(n0); 
			eps     = _mm_sub_ps(rt,n0f);
			eps2    = _mm_mul_ps(eps,eps);
			
			nnn     = _mm_slli_epi32(n0,2); 
			
			/* the tables are 16-byte aligned, we can use _mm_load_ps */			
			xmm1    = _mm_load_ps(GBtab+(_mm_extract_epi32(nnn,0)));  /* Y1,F1,G1,H1 */
			xmm2    = _mm_load_ps(GBtab+(_mm_extract_epi32(nnn,1)));  /* Y2,F2,G2,H2 */
			xmm3    = _mm_load_ps(GBtab+(_mm_extract_epi32(nnn,2)));  /* Y3,F3,G3,H3 */
			xmm4    = _mm_load_ps(GBtab+(_mm_extract_epi32(nnn,3)));  /* Y4,F4,G4,H4 */
			
			/* transpose 4*4 */
			xmm5    = _mm_unpacklo_ps(xmm1,xmm2); /* Y1,Y2,F1,F2 */
			xmm6    = _mm_unpacklo_ps(xmm3,xmm4); /* Y3,Y4,F3,F4 */
			xmm7    = _mm_unpackhi_ps(xmm1,xmm2); /* G1,G2,H1,H2 */
			xmm8    = _mm_unpackhi_ps(xmm3,xmm4); /* G3,G4,H3,H4 */
			
			Y       = _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(1,0,1,0)); /* Y1 Y2 Y3 Y4 */
			F		= _mm_shuffle_ps(xmm5,xmm6,_MM_SHUFFLE(3,2,3,2)); /* F1 F2 F3 F4 */
			G		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(1,0,1,0)); /* G1 G2 G3 G4 */
			H		= _mm_shuffle_ps(xmm7,xmm8,_MM_SHUFFLE(3,2,3,2)); /* H1 H2 H3 H4 */
			
			G       = _mm_mul_ps(G,eps); 
			H       = _mm_mul_ps(H,eps2); 
			Fp      = _mm_add_ps(F,G);
			Fp      = _mm_add_ps(Fp,H); 
			VV      = _mm_mul_ps(Fp,eps);
			VV      = _mm_add_ps(Y,VV); 
			xmm1    = _mm_mul_ps(two,H);
			FF      = _mm_add_ps(Fp,G);
			FF      = _mm_add_ps(FF,xmm1); 
			vgb     = _mm_mul_ps(qq,VV);
			fijC    = _mm_mul_ps(qq,FF);
			fijC    = _mm_mul_ps(fijC,gbscale); 
			
			dvdatmp = _mm_mul_ps(fijC,r);
			dvdatmp = _mm_add_ps(vgb,dvdatmp);
			dvdatmp = _mm_mul_ps(neg,dvdatmp);
			dvdatmp = _mm_mul_ps(dvdatmp,half); 
			dvdasum = _mm_add_ps(dvdasum,dvdatmp);
			
			xmm1    = _mm_mul_ps(dvdatmp,isaj);
			xmm1    = _mm_mul_ps(xmm1,isaj);
			dvdaj   = _mm_add_ps(dvdaj,xmm1);
			
			vcoul   = _mm_and_ps( (__m128) mask, vcoul);
			vgb     = _mm_and_ps( (__m128) mask, vgb);		
			
			vctot   = _mm_add_ps(vctot,vcoul);
			vgbtot  = _mm_add_ps(vgbtot,vgb);
			
			fscal   = _mm_sub_ps(fijC,fscal);
			fscal   = _mm_mul_ps(neg,fscal);
			fscal   = _mm_mul_ps(fscal,rinv);
			
			/* partial forces */
			t1      = _mm_mul_ps(fscal,dx); /* fx1,fx2,fx3,fx4 */
			t2      = _mm_mul_ps(fscal,dy); /* fy1,fy2,fy3,fy4 */
			t3      = _mm_mul_ps(fscal,dz); /* fz1,fz2,fz3,fz4 */
			
			/* store the partial force and dvda */
			if(offset==1) {
				xmm1 = _mm_loadl_pi(xmm1, (__m64 *) (faction+j13));
				xmm7 = _mm_load1_ps(faction+j13+2);
				
				xmm5 = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(0,0,0,0));
				xmm6 = _mm_shuffle_ps(xmm1,xmm2,_MM_SHUFFLE(1,1,1,1));
				
				xmm5 = _mm_sub_ps(xmm5,t1);
				xmm6 = _mm_sub_ps(xmm6,t2);
				xmm7 = _mm_sub_ps(xmm7,t3);
				
				_mm_store_ss(faction+j13 , xmm5);
				_mm_store_ss(faction+j13+1,xmm6);
				_mm_store_ss(faction+j13+2,xmm7);
				
				_mm_store_ss(dvda+jnr1,dvdaj);
			}
			else if(offset==2) {
				xmm1 = _mm_loadh_pi(xmm1, (__m64 *) (faction+j13)); 
				xmm2 = _mm_loadh_pi(xmm2, (__m64 *) (faction+j23));  
				
				xmm5 = _mm_load1_ps(faction+j13+2); 
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
				
				_mm_storel_pi( (__m64 *) (faction+j13), xmm5);
				_mm_storeh_pi( (__m64 *) (faction+j23), xmm5);
				
				_mm_store_ss(faction+j13+2,xmm7);
				xmm7 = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
				_mm_store_ss(faction+j23+2,xmm7);
				
				_mm_store_ss(dvda+jnr1,dvdaj);
				xmm1 = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(0,3,2,1)); 
				_mm_store_ss(dvda+jnr2,xmm1);
			}
			else {
				xmm1 = _mm_loadh_pi(xmm1, (__m64 *) (faction+j13)); 
				xmm2 = _mm_loadh_pi(xmm2, (__m64 *) (faction+j23));  
				xmm3 = _mm_loadh_pi(xmm3, (__m64 *) (faction+j33)); 
				
				xmm5 = _mm_load1_ps(faction+j13+2); 
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
				
				_mm_storel_pi( (__m64 *) (faction+j13), xmm1);
				_mm_storeh_pi( (__m64 *) (faction+j23), xmm1);
				_mm_storel_pi( (__m64 *) (faction+j33), xmm2);
				
				_mm_store_ss(faction+j13+2,xmm7); 
				xmm7 = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
				_mm_store_ss(faction+j23+2,xmm7); 
				xmm7 = _mm_shuffle_ps(xmm7,xmm7,_MM_SHUFFLE(0,3,2,1));
				_mm_store_ss(faction+j33+2,xmm7); 
				
				_mm_store_ss(dvda+jnr1,dvdaj);
				xmm1 = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(0,3,2,1)); 
				_mm_store_ss(dvda+jnr2,xmm1);
				xmm1 = _mm_shuffle_ps(dvdaj,dvdaj,_MM_SHUFFLE(1,0,3,2)); 
				_mm_store_ss(dvda+jnr3,xmm1);
			}
			
			t1      = _mm_and_ps( (__m128) mask, t1);
			t2      = _mm_and_ps( (__m128) mask, t2);
			t3      = _mm_and_ps( (__m128) mask, t3);		
			
			/* add the i force */
			fix     = _mm_add_ps(fix,t1);
			fiy     = _mm_add_ps(fiy,t2);
			fiz     = _mm_add_ps(fiz,t3);
			
		}
		
		/* fix/fiy/fiz now contain four partial terms, that all should be
		 * added to the i particle forces
		 */
		
		t1   = _mm_movehl_ps(t1,fix);
		t2   = _mm_movehl_ps(t2,fiy);
		t3   = _mm_movehl_ps(t3,fiz);
		
		fix  = _mm_add_ps(fix,t1);
		fiy  = _mm_add_ps(fiy,t2);
		fiz  = _mm_add_ps(fiz,t3);
		
		t1   = _mm_shuffle_ps(fix,fix,_MM_SHUFFLE(1,1,1,1));
		t2   = _mm_shuffle_ps(fiy,fiy,_MM_SHUFFLE(1,1,1,1));
		t3   = _mm_shuffle_ps(fiz,fiz,_MM_SHUFFLE(1,1,1,1));
		
		fix  = _mm_add_ps(fix,t1);
		fiy  = _mm_add_ps(fiy,t2);
		fiz  = _mm_add_ps(fiz,t3);
		
		xmm2 = _mm_unpacklo_ps(fix,fiy); /* fx, fy, - - */
		xmm2 = _mm_movelh_ps(xmm2,fiz);  
		xmm2 = _mm_and_ps( (__m128) maski, xmm2);
		
		/* load i force from memory */
		xmm4 = _mm_loadl_pi(xmm4, (__m64 *) (faction+ii3));
		xmm5 = _mm_load1_ps(faction+ii3+2);
		xmm4 = _mm_shuffle_ps(xmm4,xmm5,_MM_SHUFFLE(3,2,1,0));
		
		/* add to i force */
		xmm4 = _mm_add_ps(xmm4,xmm2);
		
		/* store i force to memory */
		_mm_storel_pi( (__m64 *) (faction+ii3),xmm4);
		xmm4 = _mm_shuffle_ps(xmm4,xmm4,_MM_SHUFFLE(2,2,2,2));
		_mm_store_ss(faction+ii3+2,xmm4);
		
		/* now do dvda */
		dvdatmp = _mm_movehl_ps(dvdatmp,dvdasum);
		dvdasum = _mm_add_ps(dvdasum,dvdatmp);
		dvdatmp = _mm_shuffle_ps(dvdasum,dvdasum,_MM_SHUFFLE(1,1,1,1));
		dvdasum = _mm_add_ss(dvdasum,dvdatmp);
		
		_mm_store_ss(&dva,dvdasum);
		dvda[ii] = dvda[ii] + dva*isai_f*isai_f;
		
		ggid    = gid[n];
		
		/* Coulomb potential */
		vcoul   = _mm_movehl_ps(vcoul,vctot);
		vctot   = _mm_add_ps(vctot,vcoul);
		vcoul   = _mm_shuffle_ps(vctot,vctot,_MM_SHUFFLE(1,1,1,1));
		vctot   = _mm_add_ss(vctot,vcoul);
		
		_mm_store_ss(&vct,vctot);
		Vc[ggid] = Vc[ggid] + vct;
		
		/* Store the GB potential to the work array */
		vgb     = _mm_movehl_ps(vgb,vgbtot);
		vgbtot  = _mm_add_ps(vgbtot,vgb);
		vgb     = _mm_shuffle_ps(vgbtot,vgbtot,_MM_SHUFFLE(1,1,1,1));
		vgbtot  = _mm_add_ss(vgbtot,vgb);
		
		_mm_store_ss(&vgbt,vgbtot);
		work[ggid] = work[ggid] + vgbt;
	}
	
	*outeriter       = nri;            
    *inneriter       = nj1; 
}


/*
 * Gromacs nonbonded kernel nb_kernel400nf
 * Coulomb interaction:     Generalized-Born
 * VdW interaction:         Not calculated
 * water optimization:      No
 * Calculate forces:        no
 */
void nb_kernel400nf_sse2_single(
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
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         isai,isaj,isaprod,gbscale;
    float         ix1,iy1,iz1;
    float         jx1,jy1,jz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    const int     fractshift = 12;
    const int     fractmask = 8388607;
    const int     expshift = 23;
    const int     expmask = 2139095040;
    const int     explsb = 8388608;
    float         lu;
    int           iexp,addr;
    union { unsigned int bval; float fval; } bitpattern,result;
	
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
        vctot            = 0;              
        
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
            bitpattern.fval  = rsq11;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv11           = (0.5*lu*(3.0-((rsq11*lu)*lu)));
            isaj             = invsqrta[jnr];  
            isaprod          = isai*isaj;      
            qq               = iq*charge[jnr]; 
            vcoul            = qq*rinv11;      
            qq               = isaprod*(-qq);  
            gbscale          = isaprod*gbtabscale;
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
        }
        
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}


