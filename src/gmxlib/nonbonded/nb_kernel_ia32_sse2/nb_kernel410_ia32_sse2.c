/*
 * Copyright (c) Erik Lindahl, David van der Spoel 2003
 * 
 * This file is generated automatically at compile time
 * by the program mknb in the Gromacs distribution.
 *
 * Options used when generation this file:
 * Language:         c
 * Precision:        double
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

/* get gmx_gbdata_t */
#include "../nb_kerneltype.h"


static inline __m128d
my_invrsq_pd(__m128d x)
{
	const __m128d three = (const __m128d) {3.0f, 3.0f};
	const __m128d half  = (const __m128d) {0.5f, 0.5f};
	
	__m128  t  = _mm_rsqrt_ps(_mm_cvtpd_ps(x)); /* Convert to single precision and do _mm_rsqrt_ps() */
	__m128d t1 = _mm_cvtps_pd(t); /* Convert back to double precision */
	
	/* First Newton-Rapson step, accuracy is now 24 bits */
	__m128d t2 = _mm_mul_pd(half,_mm_mul_pd(t1,_mm_sub_pd(three,_mm_mul_pd(x,_mm_mul_pd(t1,t1)))));
	
	/* Return second Newton-Rapson step, accuracy 48 bits */
	return (__m128d) _mm_mul_pd(half,_mm_mul_pd(t2,_mm_sub_pd(three,_mm_mul_pd(x,_mm_mul_pd(t2,t2)))));
}

/* to extract single integers from a __m128i datatype */
#define _mm_extract_epi64(x, imm) \
_mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

void nb_kernel410_ia32_sse2(int *           p_nri,
							int *           iinr,
							int *           jindex,
							int *           jjnr,
							int *           shift,
							double *         shiftvec,
							double *         fshift,
							int *           gid,
							double *         pos,
							double *         faction,
							double *         charge,
							double *         p_facel,
							double *         p_krf,
							double *         p_crf,
							double *         Vc,
							int *           type,
							int *           p_ntype,
							double *         vdwparam,
							double *         Vvdw,
							double *         p_tabscale,
							double *         VFtab,
							double *         invsqrta,
							double *         dvda,
							double *         p_gbtabscale,
							double *         GBtab,
							int *           p_nthreads,
							int *           count,
							void *          mtx,
							int *           outeriter,
							int *           inneriter,
							double *         work)
{
	int           nri,ntype,nthreads,offset,tj,tj2,nti;
	int           n,ii,is3,ii3,k,nj0,nj1,jnr1,jnr2,j13,j23,ggid;
	double        facel,krf,crf,tabscl,gbtabscl,vct,vdwt,nt1,nt2;
	double        shX,shY,shZ,isai_d,dva,vgbt;
	gmx_gbdata_t *gbdata;
	double *        gpol;

	__m128d       ix,iy,iz,jx,jy,jz;
	__m128d		  dx,dy,dz,t1,t2,t3;
	__m128d		  fix,fiy,fiz,rsq11,rinv,r,fscal,rt,eps,eps2;
	__m128d		  q,iq,qq,isai,isaj,isaprod,vcoul,gbscale,dvdai,dvdaj;
	__m128d       Y,F,G,H,Fp,VV,FF,vgb,fijC,dvdatmp,dvdasum,vctot,vgbtot,n0d;
	__m128d		  xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8;
	__m128d       c6,c12,Vvdw6,Vvdw12,Vvdwtmp,Vvdwtot,rinvsq,rinvsix;
	__m128d       fac,tabscale,gbtabscale;
	__m128i       n0,nnn;
	
	const __m128d neg    = {-1.0f,-1.0f};
	const __m128d zero   = {0.0f,0.0f};
	const __m128d half   = {0.5f,0.5f};
	const __m128d two    = {2.0f,2.0f};
	const __m128d three  = {3.0f,3.0f};
	const __m128d six    = {6.0f,6.0f};
	const __m128d twelwe = {12.0f,12.0f};
	
	gbdata     = (gmx_gbdata_t *)work;
	gpol       = gbdata->gpol;

	nri        = *p_nri;
	ntype      = *p_ntype;
	nthreads   = *p_nthreads; 
    facel      = (*p_facel) * (1.0 - (1.0/gbdata->gb_epsilon_solvent));       
	krf        = *p_krf;
	crf        = *p_crf;
	tabscl     = *p_tabscale;
	gbtabscl   = *p_gbtabscale;
	nj1        = 0;
	
	/* Splat variables */
	fac        = _mm_load1_pd(&facel);
	tabscale   = _mm_load1_pd(&tabscl);
	gbtabscale = _mm_load1_pd(&gbtabscl);
		
	/* Keep compiler happy */
	Vvdwtmp = _mm_setzero_pd();
	Vvdwtot = _mm_setzero_pd();
	dvdatmp = _mm_setzero_pd();
	dvdaj   = _mm_setzero_pd();
	isaj    = _mm_setzero_pd();
	vcoul   = _mm_setzero_pd();
	vgb     = _mm_setzero_pd();
	t1      = _mm_setzero_pd();
	t2      = _mm_setzero_pd();
	t3      = _mm_setzero_pd();
	xmm1    = _mm_setzero_pd();
	xmm2    = _mm_setzero_pd();
	xmm3    = _mm_setzero_pd();
	xmm4    = _mm_setzero_pd();
	jnr1    = jnr2 = 0;
	j13     = j23  = 0;
	
	for(n=0;n<nri;n++)
	{
		is3     = 3*shift[n];
		shX     = shiftvec[is3];
		shY     = shiftvec[is3+1];
		shZ     = shiftvec[is3+2];
		
		nj0     = jindex[n];      
        nj1     = jindex[n+1];  
		offset  = (nj1-nj0)%2;
		
		ii      = iinr[n];
		ii3     = ii*3;
		
		ix      = _mm_set1_pd(shX+pos[ii3+0]);
		iy      = _mm_set1_pd(shX+pos[ii3+1]);
		iz      = _mm_set1_pd(shX+pos[ii3+2]); 
		q       = _mm_set1_pd(charge[ii]);
		
		iq      = _mm_mul_pd(fac,q); 
		isai_d  = invsqrta[ii];
		isai    = _mm_load1_pd(&isai_d);
		
		nti      = 2*ntype*type[ii];
		
		fix     = _mm_setzero_pd();
		fiy     = _mm_setzero_pd();
		fiz     = _mm_setzero_pd();
		dvdasum = _mm_setzero_pd();
		vctot   = _mm_setzero_pd();
		vgbtot  = _mm_setzero_pd();
		Vvdwtot = _mm_setzero_pd();
		
		for(k=nj0;k<nj1-offset; k+=2)
		{
			jnr1    = jjnr[k];
			jnr2    = jjnr[k+1];
			
			j13     = jnr1 * 3;
			j23     = jnr2 * 3;
			
			/* Load coordinates */
			xmm1    = _mm_loadu_pd(pos+j13); /* x1 y1 */
			xmm2    = _mm_loadu_pd(pos+j23); /* x2 y2 */
			
			xmm5    = _mm_load_sd(pos+j13+2); /* z1 - */
			xmm6    = _mm_load_sd(pos+j23+2); /* z2 - */
			
			/* transpose */
			jx      = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0)); 
			jy      = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1)); 
			jz      = _mm_shuffle_pd(xmm5,xmm6,_MM_SHUFFLE2(0,0)); 
			
			/* distances */
			dx      = _mm_sub_pd(ix,jx);
			dy		= _mm_sub_pd(iy,jy);
			dz		= _mm_sub_pd(iz,jz);
			
			rsq11   = _mm_add_pd( _mm_add_pd( _mm_mul_pd(dx,dx) , _mm_mul_pd(dy,dy) ) , _mm_mul_pd(dz,dz) );
			rinv    = my_invrsq_pd(rsq11);
			
			/* Load invsqrta */
			isaj	= _mm_loadl_pd(isaj,invsqrta+jnr1);
			isaj	= _mm_loadh_pd(isaj,invsqrta+jnr2);
			isaprod = _mm_mul_pd(isai,isaj);
			
			/* Load charges */
			q		= _mm_loadl_pd(q,charge+jnr1);
			q		= _mm_loadh_pd(q,charge+jnr2);
			qq		= _mm_mul_pd(iq,q);
			
			vcoul	= _mm_mul_pd(qq,rinv);
			fscal	= _mm_mul_pd(vcoul,rinv);
			qq		= _mm_mul_pd(isaprod,qq);
			qq		= _mm_mul_pd(qq,neg);
			gbscale	= _mm_mul_pd(isaprod,gbtabscale);
			
			/* Load VdW parameters */
			tj      = nti+2*type[jnr1];
			tj2     = nti+2*type[jnr2];
			
			xmm1      = _mm_loadu_pd(vdwparam+tj);
			xmm2     = _mm_loadu_pd(vdwparam+tj2);
			c6      = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0));
			c12     = _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1));
						
			rinvsq  = _mm_mul_pd(rinv,rinv);
			
			/* Load dvdaj */
			dvdaj	= _mm_loadl_pd(dvdaj, dvda+jnr1);
			dvdaj	= _mm_loadh_pd(dvdaj, dvda+jnr2);
			
			r		= _mm_mul_pd(rsq11,rinv);
			rt		= _mm_mul_pd(r,gbscale);
			n0		= _mm_cvttpd_epi32(rt);
			n0d		= _mm_cvtepi32_pd(n0);
			eps		= _mm_sub_pd(rt,n0d);
			eps2	= _mm_mul_pd(eps,eps);
			
			nnn		= _mm_slli_epi64(n0,2);
			
			xmm1	= _mm_load_pd(GBtab+(_mm_extract_epi64(nnn,0)));   /* Y1 F1 */
			xmm2	= _mm_load_pd(GBtab+(_mm_extract_epi64(nnn,1)));   /* Y2 F2 */
			xmm3	= _mm_load_pd(GBtab+(_mm_extract_epi64(nnn,0))+2); /* G1 H1 */
			xmm4	= _mm_load_pd(GBtab+(_mm_extract_epi64(nnn,1))+2); /* G2 H2 */
			
			Y		= _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0)); /* Y1 Y2 */
			F		= _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1)); /* F1 F2 */
			G		= _mm_shuffle_pd(xmm3,xmm4,_MM_SHUFFLE2(0,0)); /* G1 G2 */
			H		= _mm_shuffle_pd(xmm3,xmm4,_MM_SHUFFLE2(1,1)); /* H1 H2 */
			
			G		= _mm_mul_pd(G,eps);
			H		= _mm_mul_pd(H,eps2);
			Fp		= _mm_add_pd(F,G);
			Fp		= _mm_add_pd(Fp,H);
			VV		= _mm_mul_pd(Fp,eps);
			VV		= _mm_add_pd(Y,VV);
			H		= _mm_mul_pd(two,H);
			FF		= _mm_add_pd(Fp,G);
			FF		= _mm_add_pd(FF,H);
			vgb		= _mm_mul_pd(qq,VV);
			fijC	= _mm_mul_pd(qq,FF);
			fijC	= _mm_mul_pd(fijC,gbscale);
			
			dvdatmp = _mm_mul_pd(fijC,r);
			dvdatmp	= _mm_add_pd(vgb,dvdatmp);
			dvdatmp = _mm_mul_pd(dvdatmp,neg);
			dvdatmp = _mm_mul_pd(dvdatmp,half);
			dvdasum	= _mm_add_pd(dvdasum,dvdatmp);
			
			xmm1	= _mm_mul_pd(dvdatmp,isaj);
			xmm1	= _mm_mul_pd(xmm1,isaj);
			dvdaj	= _mm_add_pd(dvdaj,xmm1);
			
			/* store dvda */
			_mm_storel_pd(dvda+jnr1,dvdaj);
			_mm_storeh_pd(dvda+jnr2,dvdaj);
			
			vctot	= _mm_add_pd(vctot,vcoul);
			vgbtot  = _mm_add_pd(vgbtot,vgb);
			
			/* VdW interaction */
			rinvsix = _mm_mul_pd(rinvsq,rinvsq);
			rinvsix = _mm_mul_pd(rinvsix,rinvsq);
			
			Vvdw6   = _mm_mul_pd(c6,rinvsix);
			Vvdw12  = _mm_mul_pd(c12,rinvsix);
			Vvdw12  = _mm_mul_pd(Vvdw12,rinvsix);
			Vvdwtmp = _mm_sub_pd(Vvdw12,Vvdw6);
			Vvdwtot = _mm_add_pd(Vvdwtot,Vvdwtmp);
			
			xmm1    = _mm_mul_pd(twelwe,Vvdw12);
			xmm2    = _mm_mul_pd(six,Vvdw6);
			xmm1    = _mm_sub_pd(xmm1,xmm2);
			xmm1    = _mm_mul_pd(xmm1,rinvsq);
			
			/* Scalar force */
			fscal   = _mm_sub_pd(fijC,fscal);
			fscal   = _mm_mul_pd(fscal,rinv);
			fscal   = _mm_sub_pd(xmm1,fscal);
			
			/* calculate partial force terms */
			t1		= _mm_mul_pd(fscal,dx);
			t2		= _mm_mul_pd(fscal,dy);
			t3		= _mm_mul_pd(fscal,dz);
			
			/* update the i force */
			fix		= _mm_add_pd(fix,t1);
			fiy		= _mm_add_pd(fiy,t2);
			fiz		= _mm_add_pd(fiz,t3);
			
			/* accumulate forces from memory */
			xmm1	= _mm_loadu_pd(faction+j13); /* fx1 fy1 */
			xmm2	= _mm_loadu_pd(faction+j23); /* fx2 fy2 */
			
			xmm5	= _mm_load1_pd(faction+j13+2); /* fz1 fz1 */
			xmm6	= _mm_load1_pd(faction+j23+2); /* fz2 fz2 */
			
			/* transpose */
			xmm7	= _mm_shuffle_pd(xmm5,xmm6,_MM_SHUFFLE2(0,0)); /* fz1 fz2 */
			xmm5	= _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0)); /* fx1 fx2 */
			xmm6	= _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1)); /* fy1 fy2 */
			
			/* subtract partial forces */
			xmm5	= _mm_sub_pd(xmm5,t1);
			xmm6	= _mm_sub_pd(xmm6,t2);
			xmm7	= _mm_sub_pd(xmm7,t3);
			
			xmm1	= _mm_shuffle_pd(xmm5,xmm6,_MM_SHUFFLE2(0,0)); /* fx1 fy1 */
			xmm2	= _mm_shuffle_pd(xmm5,xmm6,_MM_SHUFFLE2(1,1)); /* fy1 fy2 */
			
			/* store fx and fy */
			_mm_storeu_pd(faction+j13,xmm1);
			_mm_storeu_pd(faction+j23,xmm2);
			
			/* .. then fz */
			_mm_storel_pd(faction+j13+2,xmm7);
			_mm_storeh_pd(faction+j23+2,xmm7);
		}
		
		/* In double precision, offset can only be either 0 or 1 */
		if(offset!=0)
		{
			jnr1	= jjnr[k];
			j13		= jnr1*3;
			
			jx      = _mm_load_sd(pos+j13);
			jy      = _mm_load_sd(pos+j13+1);
			jz      = _mm_load_sd(pos+j13+2);
			
			isaj	= _mm_load_sd(invsqrta+jnr1);
			isaprod = _mm_mul_sd(isai,isaj);
			dvdaj	= _mm_load_sd(dvda+jnr1);
			q		= _mm_load_sd(charge+jnr1);
			qq      = _mm_mul_sd(iq,q);
			
			dx      = _mm_sub_sd(ix,jx);
			dy		= _mm_sub_sd(iy,jy);
			dz		= _mm_sub_sd(iz,jz);
			
			rsq11   = _mm_add_pd( _mm_add_pd( _mm_mul_pd(dx,dx) , _mm_mul_pd(dy,dy) ) , _mm_mul_pd(dz,dz) );
			rinv    = my_invrsq_pd(rsq11);
			
			vcoul	= _mm_mul_sd(qq,rinv);
			fscal	= _mm_mul_sd(vcoul,rinv);
			qq		= _mm_mul_sd(isaprod,qq);
			qq		= _mm_mul_sd(qq,neg);
			gbscale	= _mm_mul_sd(isaprod,gbtabscale);
			
			/* Load VdW parameters */
			tj      = nti+2*type[jnr1];
			
			c6      = _mm_load_sd(vdwparam+tj);
			c12     = _mm_load_sd(vdwparam+tj+1);
						
			rinvsq  = _mm_mul_sd(rinv,rinv);
			
			r		= _mm_mul_sd(rsq11,rinv);
			rt		= _mm_mul_sd(r,gbscale);
			n0		= _mm_cvttpd_epi32(rt);
			n0d		= _mm_cvtepi32_pd(n0);
			eps		= _mm_sub_sd(rt,n0d);
			eps2	= _mm_mul_sd(eps,eps);
			
			nnn		= _mm_slli_epi64(n0,2);
			
			xmm1	= _mm_load_pd(GBtab+(_mm_extract_epi64(nnn,0))); 
			xmm2	= _mm_load_pd(GBtab+(_mm_extract_epi64(nnn,1))); 
			xmm3	= _mm_load_pd(GBtab+(_mm_extract_epi64(nnn,0))+2); 
			xmm4	= _mm_load_pd(GBtab+(_mm_extract_epi64(nnn,1))+2); 
			
			Y		= _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(0,0)); 
			F		= _mm_shuffle_pd(xmm1,xmm2,_MM_SHUFFLE2(1,1)); 
			G		= _mm_shuffle_pd(xmm3,xmm4,_MM_SHUFFLE2(0,0)); 
			H		= _mm_shuffle_pd(xmm3,xmm4,_MM_SHUFFLE2(1,1)); 
			
			G		= _mm_mul_sd(G,eps);
			H		= _mm_mul_sd(H,eps2);
			Fp		= _mm_add_sd(F,G);
			Fp		= _mm_add_sd(Fp,H);
			VV		= _mm_mul_sd(Fp,eps);
			VV		= _mm_add_sd(Y,VV);
			H		= _mm_mul_sd(two,H);
			FF		= _mm_add_sd(Fp,G);
			FF		= _mm_add_sd(FF,H);
			vgb		= _mm_mul_sd(qq,VV);
			fijC	= _mm_mul_sd(qq,FF);
			fijC	= _mm_mul_sd(fijC,gbscale);
			
			dvdatmp = _mm_mul_sd(fijC,r);
			dvdatmp	= _mm_add_sd(vgb,dvdatmp);
			dvdatmp = _mm_mul_sd(dvdatmp,neg);
			dvdatmp = _mm_mul_sd(dvdatmp,half);
			dvdasum	= _mm_add_sd(dvdasum,dvdatmp);
			
			xmm1	= _mm_mul_sd(dvdatmp,isaj);
			xmm1	= _mm_mul_sd(xmm1,isaj);
			dvdaj	= _mm_add_sd(dvdaj,xmm1);
			
			/* store dvda */
			_mm_storel_pd(dvda+jnr1,dvdaj);
			
			vctot	= _mm_add_sd(vctot,vcoul);
			vgbtot  = _mm_add_sd(vgbtot,vgb);
						
			/* VdW interaction */
			rinvsix = _mm_mul_sd(rinvsq,rinvsq);
			rinvsix = _mm_mul_sd(rinvsix,rinvsq);
			
			Vvdw6   = _mm_mul_sd(c6,rinvsix);
			Vvdw12  = _mm_mul_sd(c12,rinvsix);
			Vvdw12  = _mm_mul_sd(Vvdw12,rinvsix);
			Vvdwtmp = _mm_sub_sd(Vvdw12,Vvdw6);
			Vvdwtot = _mm_add_sd(Vvdwtot,Vvdwtmp);
			
			xmm1    = _mm_mul_sd(twelwe,Vvdw12);
			xmm2    = _mm_mul_sd(six,Vvdw6);
			xmm1    = _mm_sub_sd(xmm1,xmm2);
			xmm1    = _mm_mul_sd(xmm1,rinvsq);
			
			/* Scalar force */
			fscal   = _mm_sub_sd(fijC,fscal);
			fscal   = _mm_mul_sd(fscal,rinv);
			fscal   = _mm_sub_sd(xmm1,fscal);
			
			/* calculate partial force terms */
			t1		= _mm_mul_sd(fscal,dx);
			t2		= _mm_mul_sd(fscal,dy);
			t3		= _mm_mul_sd(fscal,dz);
			
			/* update the i force */
			fix		= _mm_add_sd(fix,t1);
			fiy		= _mm_add_sd(fiy,t2);
			fiz		= _mm_add_sd(fiz,t3);
			
			/* accumulate forces from memory */
			xmm5	= _mm_load_sd(faction+j13);   /* fx */
			xmm6    = _mm_load_sd(faction+j13+1); /* fy */
			xmm7    = _mm_load_sd(faction+j13+2); /* fz */
			
			/* subtract partial forces */
			xmm5	= _mm_sub_sd(xmm5,t1);
			xmm6	= _mm_sub_sd(xmm6,t2);
			xmm7	= _mm_sub_sd(xmm7,t3);
			
			/* store forces */
			_mm_store_sd(faction+j13,xmm5);
			_mm_store_sd(faction+j13+1,xmm6);
			_mm_store_sd(faction+j13+2,xmm7);
		}
		
		/* fix/fiy/fiz now contain four partial terms, that all should be
		 * added to the i particle forces
		 */
		t1		 = _mm_unpacklo_pd(t1,fix);
		t2		 = _mm_unpacklo_pd(t2,fiy);
		t3		 = _mm_unpacklo_pd(t3,fiz);
		
		fix		 = _mm_add_pd(fix,t1);
		fiy		 = _mm_add_pd(fiy,t2);
		fiz		 = _mm_add_pd(fiz,t3);
		
		fix      = _mm_shuffle_pd(fix,fix,_MM_SHUFFLE2(1,1));
		fiy      = _mm_shuffle_pd(fiy,fiy,_MM_SHUFFLE2(1,1));
		fiz      = _mm_shuffle_pd(fiz,fiz,_MM_SHUFFLE2(1,1));
		
		/* Load i forces from memory */
		xmm1     = _mm_load_sd(faction+ii3);
		xmm2     = _mm_load_sd(faction+ii3+1);
		xmm3     = _mm_load_sd(faction+ii3+2);
		
		/* Add to i force */
		fix      = _mm_add_sd(fix,xmm1);
		fiy      = _mm_add_sd(fiy,xmm2);
		fiz      = _mm_add_sd(fiz,xmm3);
		
		/* store i forces to memory */
		_mm_store_sd(faction+ii3,fix);
		_mm_store_sd(faction+ii3+1,fiy);
		_mm_store_sd(faction+ii3+2,fiz);
		
		/* Load i shift forces from memory */
		xmm1     = _mm_load_sd(fshift+ii3);
		xmm2     = _mm_load_sd(fshift+ii3+1);
		xmm3     = _mm_load_sd(fshift+ii3+2);
		
		/* Add to i force */
		fix      = _mm_add_sd(fix,xmm1);
		fiy      = _mm_add_sd(fiy,xmm2);
		fiz      = _mm_add_sd(fiz,xmm3);
		
		/* store i shift forces to memory */
		_mm_store_sd(fshift+ii3,fix);
		_mm_store_sd(fshift+ii3+1,fiy);
		_mm_store_sd(fshift+ii3+2,fiz);
		
		/* now do dvda */
		dvdatmp  = _mm_unpacklo_pd(dvdatmp,dvdasum);
		dvdasum  = _mm_add_pd(dvdasum,dvdatmp);
		_mm_storeh_pd(&dva,dvdasum);
		dvda[ii] = dvda[ii] + dva*isai_d*isai_d;
		
		ggid	 = gid[n];
		
		/* Coulomb potential */
		vcoul	 = _mm_unpacklo_pd(vcoul,vctot);
		vctot	 = _mm_add_pd(vctot,vcoul);
		_mm_storeh_pd(&vct,vctot);
		Vc[ggid] = Vc[ggid] + vct;
		
		/* VdW potential */
		Vvdwtmp	 = _mm_unpacklo_pd(Vvdwtmp,Vvdwtot);
		Vvdwtot	 = _mm_add_pd(Vvdwtot,Vvdwtmp);
		_mm_storeh_pd(&vdwt,Vvdwtot);
		Vvdw[ggid] = Vvdw[ggid] + vdwt;
		
		/* GB potential */
		vgb  	 = _mm_unpacklo_pd(vgb,vgbtot);
		vgbtot	 = _mm_add_pd(vgbtot,vgb);
		_mm_storeh_pd(&vgbt,vgbtot);
		gpol[ggid] = gpol[ggid] + vgbt;
	}

	*outeriter   = nri;            
    *inneriter   = nj1; 	
	
	
	
	
	
}
